"""
Header codec classes for different types of FASTQ headers
Each codec handles one header type with methods to:
- Test if a FASTQ header matches the format
- Parse the header into a tags dictionary (compatible with the tags used in SAM files)
- Encode a tags dictionary back into a FASTQ header
"""

import re
from typing import Dict, Optional, List, Type


class HeaderCodec:
    """Base class for header codecs."""
    _codecs: List[Type['HeaderCodec']] = []
    priority: int = 100  # Lower = higher priority (checked first)
    
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        HeaderCodec._codecs.append(cls)
        # Sort by priority (lower = higher priority)
        HeaderCodec._codecs.sort(key=lambda x: x.priority)
    
    @classmethod
    def matches(cls, header: str) -> bool:
        """Test if this codec can handle the given header.
        
        Args:
            header: Raw header string including '@' prefix
            
        Returns:
            True if this codec can parse the header
        """
        raise NotImplementedError
    
    @classmethod
    def parse(cls, header: str, indexFileParser=None, indexFileAlias=None) -> Dict[str, str]:
        """Parse header into tags dict.
        
        Args:
            header: Raw header string (with or without '@')
            indexFileParser: Optional parser for index correction
            indexFileAlias: Optional alias for index file
            
        Returns:
            Dictionary of tag_name -> value
        """
        raise NotImplementedError
    
    @classmethod
    def encode(cls, tags: Dict[str, str]) -> str:
        """Encode tags back to header string.
        
        Args:
            tags: Dictionary of tag_name -> value
            
        Returns:
            Header string (without '@' prefix)
        """
        raise NotImplementedError
    
    @classmethod
    def get_all_codecs(cls) -> List[Type['HeaderCodec']]:
        """Get all registered codecs sorted by priority."""
        return cls._codecs.copy()


class Illumina11Codec(HeaderCodec):
    """Standard 11-element Illumina header format with space separator.
    Format: @instrument:run:flowcell:lane:tile:x:y read:filter:control:index
    """
    
    priority = 10
    
    # Pre-compiled pattern - only used for matching during detection
    PATTERN = re.compile(
        r'^@?([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):[^:]+\s+[^:]+:[^:]+:[^:]+:[^:]+$'
    )
    
    @classmethod
    def matches(cls, header: str) -> bool:
        return cls.PATTERN.match(header) is not None
    
    @classmethod
    def parse(cls, header: str, indexFileParser=None, indexFileAlias=None) -> Dict[str, str]:
        # Remove @ prefix if present
        if header[0] == '@':
            header = header[1:]
        
        # Split on space to separate the two parts
        space_idx = header.find(' ')
        if space_idx == -1:
            # No space found - malformed header
            raise ValueError(f"Header does not match Illumina11 format: {header}")
        
        part1 = header[:space_idx]
        part2 = header[space_idx + 1:]
        
        # Split each part by colon
        p1_fields = part1.split(':')
        p2_fields = part2.split(':')
        
        if len(p1_fields) != 7 or len(p2_fields) != 4:
            raise ValueError(f"Header does not match Illumina11 format: {header}")
        
        tags = {
            'Is': p1_fields[0],
            'RN': p1_fields[1], 
            'Fc': p1_fields[2],
            'La': p1_fields[3],
            'Ti': p1_fields[4],
            'CX': p1_fields[5],
            'CY': p1_fields[6],
            'RP': p2_fields[0],
            'Fi': p2_fields[1],
            'CN': p2_fields[2],
            'aa': p2_fields[3]
        }
        
        # If index file parser provided, do index correction
        if indexFileParser is not None and indexFileAlias is not None:
            index_seq = p2_fields[3]
            try:
                int(index_seq)  # Check if numeric: if so, no lookup needed
            except ValueError:
                result = indexFileParser.getIndexCorrectedBarcodeAndHammingDistance(
                    alias=indexFileAlias, barcode=index_seq)
                index_id, corrected, _ = result
                if corrected is not None:
                    tags['aA'] = corrected
                    tags['aI'] = index_id
                else:
                    # The index is not present in the index file
                    pass

        return tags
    
    @classmethod
    def encode(cls, tags: Dict[str, str]) -> str:
        # Format: instrument:run:flowcell:lane:tile:x:y read:filter:control:index
        return (
            f"{tags.get('Is', '')}:{tags.get('RN', '')}:{tags.get('Fc', '')}:"
            f"{tags.get('La', '')}:{tags.get('Ti', '')}:{tags.get('CX', '')}:"
            f"{tags.get('CY', '')} "
            f"{tags.get('RP', '')}:{tags.get('Fi', '')}:{tags.get('CN', '')}:"
            f"{tags.get('aa', 'N')}"
        )


class Illumina10Codec(HeaderCodec):
    """Illumina header without index (10 elements): @instrument:run:flowcell:lane:tile:x:y:read:filter:control"""
    
    priority = 20
    
    # Pre-compiled pattern for 10 elements (no index)
    PATTERN = re.compile(
        r'^@?([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+)$'
    )
    
    @classmethod
    def matches(cls, header: str) -> bool:
        # Must match 10-element pattern but NOT be 11-element
        if not cls.PATTERN.match(header):
            return False
        # Check it's not actually an 11-element header by looking at content
        # 11-element headers have one more colon-separated field
        content = header[1:] if header.startswith('@') else header
        parts = content.split(':')
        return len(parts) == 10
    
    @classmethod
    def parse(cls, header: str, indexFileParser=None, indexFileAlias=None) -> Dict[str, str]:
        match = cls.PATTERN.match(header)
        if not match:
            raise ValueError(f"Header does not match Illumina10 format: {header}")
        
        groups = match.groups()
        return {
            'Is': groups[0],
            'RN': groups[1],
            'Fc': groups[2],
            'La': groups[3],
            'Ti': groups[4],
            'CX': groups[5],
            'CY': groups[6],
            'RP': groups[7],
            'Fi': groups[8],
            'CN': groups[9],
            'aa': 'N'
        }
    
    @classmethod
    def encode(cls, tags: Dict[str, str]) -> str:
        
        return (
            f"{tags.get('Is', '')}:{tags.get('RN', '')}:{tags.get('Fc', '')}:"
            f"{tags.get('La', '')}:{tags.get('Ti', '')}:{tags.get('CX', '')}:"
            f"{tags.get('CY', '')}:{tags.get('RP', '')}:{tags.get('Fi', '')}:"
            f"{tags.get('CN', '')}"
        )


class Illumina7Codec(HeaderCodec):
    """Header with 7 elements: @instrument:run:flowcell:lane:tile:x:y"""
    
    priority = 30
    
    # Pre-compiled pattern for 7 elements
    PATTERN = re.compile(
        r'^@?([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+):([^:]+)$'
    )
    
    @classmethod
    def matches(cls, header: str) -> bool:
        # Must match 7-element pattern exactly
        if not cls.PATTERN.match(header):
            return False
        content = header[1:] if header.startswith('@') else header
        parts = content.split(':')
        return len(parts) == 7
    
    @classmethod
    def parse(cls, header: str, indexFileParser=None, indexFileAlias=None) -> Dict[str, str]:
        match = cls.PATTERN.match(header)
        if not match:
            raise ValueError(f"Header does not match Illumina7 format: {header}")
        
        groups = match.groups()
        return {
            'Is': groups[0],
            'RN': groups[1],
            'Fc': groups[2],
            'La': groups[3],
            'Ti': groups[4],
            'CX': groups[5],
            'CY': groups[6],
            'RP': '1',
            'Fi': '0',
            'CN': '0',
            'aa': 'N'
        }
    
    @classmethod
    def encode(cls, tags: Dict[str, str]) -> str:
        return (
            f"{tags.get('Is', '')}:{tags.get('RN', '')}:{tags.get('Fc', '')}:"
            f"{tags.get('La', '')}:{tags.get('Ti', '')}:{tags.get('CX', '')}:"
            f"{tags.get('CY', '')}"
        )


class ScmoHeaderCodec(HeaderCodec):
    """SCMO header format: key:value pairs separated by semicolons"""
    
    priority = 40
    
    @classmethod
    def matches(cls, header: str) -> bool:
        # There are : in the header and at least one semicolon
        if ';' not in header or ':' not in header:
            return False
        return True
        
    
    @classmethod
    def parse(cls, header: str, indexFileParser=None, indexFileAlias=None) -> Dict[str, str]:
        content = header[1:] if header.startswith('@') else header
        content = content.strip()
        tags = {}
        
        for kv in content.split(';'):
            if ':' in kv:
                key, value = kv.split(':', 1)
                # Strip leading '@' from values (backwards compatibility)
                value = value.strip()
                if value.startswith('@'):
                    value = value[1:]
                tags[key.strip()] = value
        
        return tags
    
    @classmethod
    def encode(cls, tags: Dict[str, str]) -> str:
        return ';'.join(f"{k}:{v}" for k, v in tags.items())


class ThreePrimeDecoderCodec(HeaderCodec):
    """3-DEC format: @Cluster_s_lane_tile_readpair"""
    
    priority = 50
    
    @classmethod
    def matches(cls, header: str) -> bool:
        return header.startswith('@Cluster') and header.count('_') == 4
    
    @classmethod
    def parse(cls, header: str, indexFileParser=None, indexFileAlias=None) -> Dict[str, str]:
        # Format: @Cluster_s_1_1101_2
        parts = header.split('_')
        if len(parts) != 5:
            raise ValueError(f"Invalid 3DEC format: {header}")
        
        if parts[1] != 's':
            raise ValueError(f"Invalid 3DEC format (expected 's' at position 2): {header}")
        
        return {
            'Is': 'UNK',
            'RN': 'UNK',
            'Fc': 'UNK',
            'La': parts[2],
            'Ti': parts[3],
            'CX': '-1',
            'CY': '-1',
            'RP': parts[4],
            'Fi': '0',
            'CN': '0',
            'aa': 'N'
        }
    
    @classmethod
    def encode(cls, tags: Dict[str, str]) -> str:
        return f"Cluster_s_{tags.get('La', '')}_{tags.get('Ti', '')}_{tags.get('RP', '')}"


class FallbackCodec(HeaderCodec):
    """Fallback codec for non defined headers
    Stores the entire original header in 'oh' tag.
    """
    
    priority = 999  # Lowest priority, checked last
    
    @classmethod
    def matches(cls, header: str) -> bool:
        return True  # Always matches as fallback
    
    @classmethod
    def parse(cls, header: str, indexFileParser=None, indexFileAlias=None) -> Dict[str, str]:
        """Parse unknown format, store original header."""
        content = header[1:] if header.startswith('@') else header
        # Store the first part before any whitespace as 'oh'
        oh_value = content.replace(';', '').split()[0]
        return {'oh': oh_value}
    
    @classmethod
    def encode(cls, tags: Dict[str, str]) -> str:
        return tags.get('oh', 'unknown')


class HeaderCodecRegistry:
    """Registry that auto-detects and caches the appropriate read header codec."""
    
    def __init__(self):
        self._codec_class: Optional[Type[HeaderCodec]] = None
        self._detected: bool = False
        self._parse_func = None # Cache the parse function to avoid attribute lookups
    
    def detect(self, sample_header: str) -> Type[HeaderCodec]:
        """Detect the appropriate codec from a sample header.
        
        Args:
            sample_header: header to use for format detection
            
        Returns:
            The detected codec class
        """
        for codec_class in HeaderCodec.get_all_codecs():
            if codec_class.matches(sample_header):
                self._codec_class = codec_class
                self._detected = True
                self._parse_func = codec_class.parse # Cache the parse method to avoid attribute lookup overhead
                return codec_class
        
        # This should never happen as long as we have codecs registered
        raise ValueError(f"No codec found for header: {sample_header}")
    
    def parse(self, header: str, indexFileParser=None, indexFileAlias=None) -> Dict[str, str]:
        """Parse a header using the detected codec.
        
        Args:
            header: Header string to parse
            indexFileParser: Optional parser for index correction
            indexFileAlias: Optional alias for index file
            
        Returns:
            Dictionary of tag_name -> value
        """
        # Use cached parse function for speed (avoids attribute lookup, which takes a bit more time)
        return self._parse_func(header, indexFileParser, indexFileAlias)
    
    def encode(self, tags: Dict[str, str]) -> str:
        """Encode tags back to header.
        
        Args:
            tags: Dictionary of tag_name -> value
            
        Returns:
            Header string (without '@' prefix)
        """
        if not self._detected:
            raise RuntimeError("Codec not detected. Call detect() first.")
        return self._codec_class.encode(tags)
    
    @property
    def codec_class(self) -> Optional[Type[HeaderCodec]]:
        return self._codec_class
    
    @property
    def codec_name(self) -> Optional[str]:
        return self._codec_class.__name__ if self._codec_class else None
    
    @property
    def is_detected(self) -> bool:
        return self._detected
