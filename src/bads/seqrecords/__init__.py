from dataclasses import dataclass, field
import re
from bads.sequences import DNASequence, RNASequence, ProteinSequence


@dataclass
class SeqRecord:
    seq: DNASequence | RNASequence | ProteinSequence
    header: str = ""

    @property
    def accession(self) -> str | None:
        """Parses the header for an accession. If successful, returns a string
        otherwise returns None."""
        accession_re = r"[A-Z]{2}_\d+\.\d+"
        match = re.search(accession_re, self.header)
        return match.group()
