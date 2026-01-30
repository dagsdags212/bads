import pytest
from bads.sequences import DNASequence, ProteinSequence
from bads.sequences.sequences import RNASequence, TranslationError


@pytest.fixture
def valid_dna_sequence():
    return "ACGT"


@pytest.fixture
def valid_rna_sequence():
    return "ACGU"


# ===== BaseSequence Tests =====


def test_valid_dna_sequence(valid_dna_sequence):
    dna = DNASequence(valid_dna_sequence)
    assert len(dna) == 4


def test_valid_rna_sequence(valid_rna_sequence):
    rna = RNASequence(valid_rna_sequence)
    assert len(rna) == 4


def test_invalid_dna_sequence(valid_rna_sequence):
    with pytest.raises(ValueError):
        DNASequence(valid_rna_sequence)


def test_invalid_rna_sequence(valid_dna_sequence):
    with pytest.raises(ValueError):
        RNASequence(valid_dna_sequence)


def test_dna_repr():
    seq = DNASequence("ACGTACGTACGTACGT")
    repr_str = repr(seq)
    assert "DNASequence" in repr_str
    assert "length=16" in repr_str
    assert "ACGTACGTA...CGT" in repr_str


def test_dna_repr_short():
    seq = DNASequence("ACGT")
    repr_str = repr(seq)
    assert "DNASequence" in repr_str
    assert "length=4" in repr_str


def test_dna_str():
    seq = DNASequence("ACGT")
    assert str(seq) == "ACGT"


def test_dna_hash():
    seq1 = DNASequence("ACGT")
    seq2 = DNASequence("ACGT")
    seq3 = DNASequence("TGCA")
    assert hash(seq1) == hash(seq2)
    assert hash(seq1) != hash(seq3)


def test_dna_iter():
    seq = DNASequence("ACGT")
    result = "".join([base for base in seq])
    assert result == "ACGT"


def test_dna_contains():
    seq = DNASequence("ACGTACGT")
    assert "ACG" in seq
    assert "CGT" in seq
    assert "XYZ" not in seq


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


def test_dna_getitem_single():
    seq = DNASequence("ACGT")
    assert seq[0] == "A"
    assert seq[1] == "C"
    assert seq[-1] == "T"


def test_frequency():
    seq = DNASequence("AACGT")
    freq = seq.frequency()
    assert freq == {"A": 2, "C": 1, "G": 1, "T": 1}


def test_reverse():
    seq = DNASequence("ACGT")
    rev_seq = seq.reverse()
    assert isinstance(rev_seq, DNASequence)
    assert rev_seq.seq == "TGCA"


# ===== DNASequence Tests =====


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


def test_dna_complement():
    seq = DNASequence("ACGT")
    comp = seq.complement()
    assert isinstance(comp, DNASequence)
    assert comp.seq == "TGCA"


def test_dna_complement_all_bases():
    seq = DNASequence("AAATTTGGGCCC")
    comp = seq.complement()
    assert comp.seq == "TTTAAACCCGGG"


def test_dna_transcribe():
    seq = DNASequence("ACGT")
    rna = seq.transcribe()
    assert isinstance(rna, RNASequence)
    assert rna.seq == "ACGU"


def test_dna_transcribe_longer():
    seq = DNASequence("ATGATGTAA")
    rna = seq.transcribe()
    assert rna.seq == "AUGAUGUAA"


def test_dna_translate():
    # ATG = M (methionine), TAA = _ (stop)
    # Stop codons are not included in the protein sequence
    seq = DNASequence("ATGTAA")
    protein = seq.translate()
    assert isinstance(protein, ProteinSequence)
    assert protein.seq == "M"


def test_dna_translate_to_stop():
    # ATG = M (methionine), TAA = _ (stop)
    # With to_stop=True, translation stops at first stop codon (not included)
    seq = DNASequence("ATGTAA")
    protein = seq.translate(to_stop=True)
    assert isinstance(protein, ProteinSequence)
    assert protein.seq == "M"


def test_dna_translate_longer():
    # ATG=M, GCT=A, TAA=stop
    # By default, stop codons are skipped (not included)
    seq = DNASequence("ATGGCTTAA")
    protein = seq.translate()
    assert isinstance(protein, ProteinSequence)
    assert protein.seq == "MA"


def test_dna_translate_longer_to_stop():
    # ATG=M, GCT=A, TAA=stop, GGG=G (after stop)
    # With to_stop=True, stops at TAA and doesn't translate GGG
    seq = DNASequence("ATGGCTTAAGGG")
    protein = seq.translate(to_stop=True)
    assert protein.seq == "MA"


def test_dna_translate_multiple_stops():
    # ATG=M, GCT=A, TAA=stop, TTT=F, TAG=stop
    # Without to_stop, skips both stop codons
    seq = DNASequence("ATGGCTTAATTTTAG")
    protein = seq.translate()
    assert protein.seq == "MAF"


def test_dna_translate_multiple_stops_to_stop():
    # ATG=M, GCT=A, TAA=stop, TTT=F, TAG=stop
    # With to_stop=True, stops at first stop (TAA)
    seq = DNASequence("ATGGCTTAATTTTAG")
    protein = seq.translate(to_stop=True)
    assert protein.seq == "MA"


def test_dna_translate_cds():
    # Valid CDS: starts with ATG (M), ends with TAA (stop)
    seq = DNASequence("ATGGCTTAA")
    protein = seq.translate(cds=True)
    assert isinstance(protein, ProteinSequence)
    assert protein.seq == "MA"  # Stop codon excluded in CDS mode


def test_dna_translate_cds_no_start_codon():
    # Invalid CDS: doesn't start with ATG
    seq = DNASequence("GCTTAA")
    with pytest.raises(TranslationError, match="not a start codon"):
        seq.translate(cds=True)


def test_dna_translate_cds_no_stop_codon():
    # Invalid CDS: no stop codon
    seq = DNASequence("ATGGCT")
    with pytest.raises(TranslationError, match="does not contain a stop codon"):
        seq.translate(cds=True)


def test_dna_translate_with_offset():
    # From position 1: TGG=W, CTT=L
    seq = DNASequence("ATGGCTT")
    protein = seq.translate(initial_position=1)
    assert protein.seq == "WL"


def test_dna_translate_no_start_codon():
    # GCT=A, TAA=stop (no ATG start codon) - but this is fine without cds=True
    seq = DNASequence("GCTTAA")
    protein = seq.translate()
    assert protein.seq == "A"  # Stop codon not included


def test_dna_reading_frames():
    # Short sequence for testing
    seq = DNASequence("ATGGCTTAA")
    frames = seq.reading_frames()
    assert len(frames) == 6  # 3 forward + 3 reverse complement
    # Each frame returns a ProteinSequence
    assert all(isinstance(frame, ProteinSequence) for frame in frames)
    # First frame should find "MA" (stop codon not included)
    assert frames[0].seq == "MA"


def test_dna_codon_usage():
    # Test with leucine: CTT, CTC, CTA, CTG, TTA, TTG all code for L
    seq = DNASequence("CTTCTCCTACTTTTG")
    usage = seq.codon_usage("L")
    assert "CTT" in usage
    assert "CTC" in usage
    assert "CTA" in usage
    # Should be normalized frequencies (sum to 1.0)
    assert abs(sum(usage.values()) - 1.0) < 0.001


def test_dna_codon_usage_no_matches():
    seq = DNASequence("ATGATGATG")
    usage = seq.codon_usage("L")
    assert usage == {}


# ===== RNASequence Tests =====


def test_rna_complement():
    seq = RNASequence("ACGU")
    comp = seq.complement()
    assert isinstance(comp, RNASequence)
    assert comp.seq == "UGCA"


def test_rna_complement_all_bases():
    seq = RNASequence("AAAUUUGGGCCC")
    comp = seq.complement()
    assert comp.seq == "UUUAAACCCGGG"


def test_rna_translate():
    # AUG = M (methionine), UAA = _ (stop)
    # Stop codons are not included in the protein sequence
    seq = RNASequence("AUGUAA")
    protein = seq.translate()
    assert isinstance(protein, ProteinSequence)
    assert protein.seq == "M"


def test_rna_translate_to_stop():
    # AUG = M (methionine), UAA = _ (stop)
    # With to_stop=True, translation stops at first stop codon (not included)
    seq = RNASequence("AUGUAA")
    protein = seq.translate(to_stop=True)
    assert isinstance(protein, ProteinSequence)
    assert protein.seq == "M"


def test_rna_translate_longer():
    # AUG=M, GCU=A, UAA=stop
    # By default, stop codons are skipped (not included)
    seq = RNASequence("AUGGCUUAA")
    protein = seq.translate()
    assert isinstance(protein, ProteinSequence)
    assert protein.seq == "MA"


def test_rna_translate_longer_to_stop():
    # AUG=M, GCU=A, UAA=stop, GGG=G (after stop)
    # With to_stop=True, stops at UAA and doesn't translate GGG
    seq = RNASequence("AUGGCUUAAGGG")
    protein = seq.translate(to_stop=True)
    assert protein.seq == "MA"


def test_rna_translate_multiple_stops():
    # AUG=M, GCU=A, UAA=stop, UUU=F, UAG=stop
    # Without to_stop, skips both stop codons
    seq = RNASequence("AUGGCUUAAUUUUAG")
    protein = seq.translate()
    assert protein.seq == "MAF"


def test_rna_translate_multiple_stops_to_stop():
    # AUG=M, GCU=A, UAA=stop, UUU=F, UAG=stop
    # With to_stop=True, stops at first stop (UAA)
    seq = RNASequence("AUGGCUUAAUUUUAG")
    protein = seq.translate(to_stop=True)
    assert protein.seq == "MA"


def test_rna_translate_cds():
    # Valid CDS: starts with AUG (M), ends with UAA (stop)
    seq = RNASequence("AUGGCUUAA")
    protein = seq.translate(cds=True)
    assert isinstance(protein, ProteinSequence)
    assert protein.seq == "MA"  # Stop codon excluded in CDS mode


def test_rna_translate_cds_no_start_codon():
    # Invalid CDS: doesn't start with AUG
    seq = RNASequence("GCUUAA")
    with pytest.raises(TranslationError, match="not a start codon"):
        seq.translate(cds=True)


def test_rna_translate_cds_no_stop_codon():
    # Invalid CDS: no stop codon
    seq = RNASequence("AUGGCU")
    with pytest.raises(TranslationError, match="does not contain a stop codon"):
        seq.translate(cds=True)


def test_rna_codon_usage():
    # Test with leucine: CUU, CUC, CUA, CUG, UUA all code for L
    seq = RNASequence("CUUCUCCUACUGUUA")
    usage = seq.codon_usage("L")
    assert "CUU" in usage
    assert "CUC" in usage
    assert "CUA" in usage
    assert "CUG" in usage
    assert "UUA" in usage
    # Should be normalized frequencies (5 codons, each 0.2)
    assert abs(sum(usage.values()) - 1.0) < 0.001
    assert len(usage) == 5


# ===== ProteinSequence Tests =====


def test_valid_protein_sequence():
    protein = ProteinSequence("ACDEFGHIKLMNPQRSTVWY")
    assert len(protein) == 20


def test_invalid_protein_sequence():
    with pytest.raises(ValueError):
        ProteinSequence("ACGTX")  # X is not a standard amino acid


def test_protein_reverse():
    seq = ProteinSequence("MKLI")
    rev = seq.reverse()
    assert isinstance(rev, ProteinSequence)
    assert rev.seq == "ILKM"


def test_protein_frequency():
    seq = ProteinSequence("AACDE")
    freq = seq.frequency()
    assert freq == {"A": 2, "C": 1, "D": 1, "E": 1}


# ===== find_proteins Tests =====


def test_find_proteins_single():
    # ATG=M, GCT=A, TAA=stop -> protein "MA"
    seq = DNASequence("ATGGCTTAA")
    proteins = seq.find_proteins()
    assert len(proteins) == 1
    assert proteins[0].seq == "MA"


def test_find_proteins_nested_starts():
    # ATG=M, ATG=M, TAA=stop -> proteins "MM" and "M"
    seq = DNASequence("ATGATGTAA")
    proteins = seq.find_proteins()
    assert len(proteins) == 2
    assert proteins[0].seq == "MM"
    assert proteins[1].seq == "M"


def test_find_proteins_multiple_separate():
    # ATG=M, TAA=stop, ATG=M, GCT=A, TAA=stop -> proteins "MA" and "M"
    seq = DNASequence("ATGTAAATGGCTTAA")
    proteins = seq.find_proteins()
    assert len(proteins) == 2
    # Sorted by length descending
    assert proteins[0].seq == "MA"
    assert proteins[1].seq == "M"


def test_find_proteins_min_length():
    # ATG=M, ATG=M, TAA=stop -> proteins "MM" and "M"
    # With min_length=2, only "MM" should be returned
    seq = DNASequence("ATGATGTAA")
    proteins = seq.find_proteins(min_length=2)
    assert len(proteins) == 1
    assert proteins[0].seq == "MM"


def test_find_proteins_no_start_codon():
    # GCT=A, TAA=stop -> no proteins (no M start)
    seq = DNASequence("GCTTAA")
    proteins = seq.find_proteins()
    assert len(proteins) == 0


def test_find_proteins_no_stop_codon():
    # ATG=M, GCT=A -> no stop codon, proteins not captured
    seq = DNASequence("ATGGCT")
    proteins = seq.find_proteins()
    assert len(proteins) == 0


def test_find_proteins_sorted_by_length():
    # Create sequence with proteins of varying lengths
    # ATG=M, GCT=A, GCT=A, ATG=M, TAA=stop -> "MAAM" and "M"
    seq = DNASequence("ATGGCTGCTATGTAA")
    proteins = seq.find_proteins()
    # Should be sorted longest first
    lengths = [len(p) for p in proteins]
    assert lengths == sorted(lengths, reverse=True)


def test_find_proteins_rna():
    # AUG=M, GCU=A, UAA=stop -> protein "MA"
    seq = RNASequence("AUGGCUUAA")
    proteins = seq.find_proteins()
    assert len(proteins) == 1
    assert proteins[0].seq == "MA"
