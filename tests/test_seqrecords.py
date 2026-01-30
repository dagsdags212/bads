import pytest
from bads.sequences import DNASequence, RNASequence, ProteinSequence
from bads.seqrecords import SeqRecord


# ===== SeqRecord Initialization Tests =====


def test_seqrecord_with_dna():
    dna = DNASequence("ATGGCTTAA")
    record = SeqRecord(seq=dna, header=">NC_000001.1 Example gene")
    assert record.seq == dna
    assert record.header == ">NC_000001.1 Example gene"


def test_seqrecord_with_rna():
    rna = RNASequence("AUGGCUUAA")
    record = SeqRecord(seq=rna, header=">NM_001234.2 Example transcript")
    assert record.seq == rna
    assert record.header == ">NM_001234.2 Example transcript"


def test_seqrecord_with_protein():
    protein = ProteinSequence("MAPLK")
    record = SeqRecord(seq=protein, header=">NP_001234.1 Example protein")
    assert record.seq == protein
    assert record.header == ">NP_001234.1 Example protein"


def test_seqrecord_default_header():
    dna = DNASequence("ACGT")
    record = SeqRecord(seq=dna)
    assert record.header == ""


# ===== Accession Parsing Tests =====


def test_accession_ncbi_format():
    dna = DNASequence("ATGGCTTAA")
    record = SeqRecord(seq=dna, header=">NC_000001.11 Homo sapiens chromosome 1")
    assert record.accession == "NC_000001.11"


def test_accession_nm_format():
    rna = RNASequence("AUGGCUUAA")
    record = SeqRecord(seq=rna, header=">NM_001234.5 Some transcript")
    assert record.accession == "NM_001234.5"


def test_accession_np_format():
    protein = ProteinSequence("MAPLK")
    record = SeqRecord(seq=protein, header=">NP_005678.2 Some protein")
    assert record.accession == "NP_005678.2"


def test_accession_xm_format():
    dna = DNASequence("ACGT")
    record = SeqRecord(seq=dna, header=">XM_012345.1 Predicted transcript")
    assert record.accession == "XM_012345.1"


def test_accession_in_middle_of_header():
    dna = DNASequence("ACGT")
    record = SeqRecord(seq=dna, header="gene_name NC_000001.1 description")
    assert record.accession == "NC_000001.1"


def test_accession_no_match():
    dna = DNASequence("ACGT")
    record = SeqRecord(seq=dna, header=">invalid_header_no_accession")
    with pytest.raises(AttributeError):
        _ = record.accession


def test_accession_empty_header():
    dna = DNASequence("ACGT")
    record = SeqRecord(seq=dna, header="")
    with pytest.raises(AttributeError):
        _ = record.accession
