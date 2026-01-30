from pathlib import Path
import pytest
from bads.sequences import DNASequence
from bads.seqrecords import SeqRecord
from bads.io import read_fasta


@pytest.fixture
def fasta_fp():
    return Path(__file__).parent / "data" / "poliovirus.fasta"


def test_read_fasta(fasta_fp):
    genome = read_fasta(fasta_fp)[0]
    header = "NC_002058.3 Poliovirus, complete genome"

    assert isinstance(genome, SeqRecord)
    assert isinstance(genome.seq, DNASequence)
    assert isinstance(genome.header, str)
    assert genome.header == header
    assert genome.accession == "NC_002058.3"
