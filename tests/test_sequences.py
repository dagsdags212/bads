import pytest
from bads.sequences import DNASequence


@pytest.fixture
def valid_dna_sequence():
    return "ACGT"


@pytest.fixture
def invalid_dna_sequence():
    return "AUCXT"


def test_valid_dna_sequence(valid_dna_sequence):
    dna = DNASequence(valid_dna_sequence)
    assert len(dna) == 4


def test_invalid_dna_sequence(invalid_dna_sequence):
    with pytest.raises(ValueError):
        DNASequence(invalid_dna_sequence)


@pytest.mark.parametrize(
    "sequence,start,end,expected",
    [
        ("ACGT", 0, 1, "A"),
        ("ACGT", 0, 2, "AC"),
        ("ACGT", 0, 3, "ACG"),
    ],
)
def test_dna_indexing(sequence, start, end, expected):
    s = DNASequence(sequence)
    assert s[start:end] == expected


@pytest.mark.parametrize(
    "sequence,gc_content",
    [
        ("GGCC", 1.0),
        ("ACGT", 0.50),
        ("GATT", 0.25),
    ],
)
def test_gc_content(sequence, gc_content):
    s = DNASequence(sequence)
    assert s.gc_content() == gc_content
