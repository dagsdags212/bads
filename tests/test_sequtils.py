import pytest
from bads.sequtils import search_first_occurrence, search_all_occurrences, BoyerMoore


# Pattern matching algorithms
@pytest.fixture
def seq() -> str:
    return "ATAGAATAGATAATAGTC"

def test_search_first_occurrence(seq):
    assert search_first_occurrence(seq, "ATAG") == 0
    assert search_first_occurrence(seq, "TAGA") == 1
    assert search_first_occurrence(seq, "ATA") == 0
    assert search_first_occurrence(seq, "TC") == 16
    assert search_first_occurrence(seq, "XYZ") == -1

def test_search_all_occurrences(seq):
    assert search_all_occurrences(seq, "ATAG") == [0, 5, 12]
    assert search_all_occurrences(seq, "TAGA") == [1, 6]
    assert search_all_occurrences(seq, "ATA") == [0, 5, 9, 12]
    assert search_all_occurrences(seq, "TC") == [16]
    assert search_all_occurrences(seq, "XYZ") == []

# Boyer-Moore
def test_boyer_moore():
    bm = BoyerMoore("ACTG", "ACCA")
    match = bm.search("ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC")
    assert match == [5, 13, 23, 37]
