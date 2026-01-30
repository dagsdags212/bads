"""
This module provides functions for reading and writing FASTA files.
"""

from pathlib import Path
from bads.sequences import DNASequence, RNASequence, ProteinSequence
from bads.alphabet import (
    DNA_AMBIGUOUS_ALPHABET,
    RNA_AMBIGUOUS_ALPHABET,
    PROTEIN_ALPHABET,
)
from bads.seqrecords import SeqRecord


def read_fasta(fp: Path) -> list[SeqRecord]:
    """Reads a FASTA file and returns a dictionary of sequences."""
    sequences = {}
    with open(fp, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                sequences[name] = ""
            else:
                sequences[name] += line

    records = []
    for header, sequence in sequences.items():
        if set(sequence) <= DNA_AMBIGUOUS_ALPHABET:
            try:
                seq_obj = DNASequence(sequence)
            except ValueError as e:
                print("Cannot parse as DNASequence")
        elif set(sequence) <= RNA_AMBIGUOUS_ALPHABET:
            try:
                seq_obj = RNASequence(sequence)
            except ValueError as e:
                print("Cannot parse as RNASequence")
        elif set(sequence) <= PROTEIN_ALPHABET:
            try:
                seq_obj = ProteinSequence(sequence)
            except ValueError as e:
                print("Cannot parse as ProteinSequence")

        records.append(SeqRecord(header=header, seq=seq_obj))

    return records
