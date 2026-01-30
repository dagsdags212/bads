from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Iterator
from typing import Self

from bads.sequences.maps import DNA_COMPLEMENT, DNA_TO_AA, RNA_COMPLEMENT, RNA_TO_AA


class BaseSequence(ABC):
    def __init__(self, seq: str) -> None:
        self.seq = self.validate_sequence(seq, self.__class__.alphabet())

    def __repr__(self) -> str:
        seq_repr = self.seq[:9] + "..." + self.seq[-3:]
        return f"{self.__class__.__name__}(length={len(self.seq)}, seq={seq_repr})"

    def __len__(self):
        return len(self.seq)

    def __hash__(self):
        return hash(self.seq)

    def __iter__(self) -> Iterator:
        return iter(self.seq)

    def __getitem__(self, key: int | slice) -> str:
        if isinstance(key, slice):
            s = [self.seq[i] for i in range(*key.indices(len(self.seq)))]
            return "".join(s)
        return self.seq[key]

    def __str__(self) -> str:
        return self.seq

    def __contains__(self, s: str) -> bool:
        return s in self.seq

    @classmethod
    @abstractmethod
    def alphabet(cls) -> str: ...

    def validate_sequence(self, seq: str, alphabet: str) -> str:
        valid_len = sum([seq.count(base) for base in alphabet])
        if valid_len != len(seq):
            invalid_chars = set(seq).difference(alphabet)
            raise ValueError(
                f"Invalid characters for {self.__class__.__name__}: {invalid_chars}"
            )
        return seq

    def reverse(self) -> Self:
        return type(self)(seq=self.seq[::-1])

    def frequency(self) -> dict[str, int]:
        """Returns counts of unique characters in a given sequence."""
        counts = {}
        for s in self.seq.upper():
            if s in counts:
                counts[s] += 1
            else:
                counts[s] = 1
        return counts


class TranslatableSequence(BaseSequence):
    @classmethod
    def codon_map(cls) -> dict[str, str]:
        raise NotImplementedError("`codom_map` must be implemented for ", cls.__name__)

    @classmethod
    def _translate_codon(cls, codon: str) -> str:
        """Translates a codon into an amino acid using an internal
        dictionary with the standard genetic code."""
        assert len(codon) == 3, "codon must be a three-character string"
        codon_map = cls.codon_map()
        if codon not in codon_map:
            raise ValueError("Invalid codon")

        return codon_map[codon]

    def translate(self, initial_position: int = 0) -> ProteinSequence:
        """Translates a DNA sequence into an amino acid sequence."""
        seq_aa = ""
        for pos in range(initial_position, len(self.seq) - 2, 3):
            codon = self.seq[pos : pos + 3]
            amino_acid = self._translate_codon(codon)
            if amino_acid == "_":
                break
            seq_aa += amino_acid

        return ProteinSequence(seq=seq_aa)

    def codon_usage(self, amino_acid: str) -> dict[str, int]:
        """Provides the frequency of each codon encoding a given
        amino acid in a DNA sequence."""
        codon_counts = {}
        total = 0
        for i in range(0, len(self.seq) - 2, 3):
            codon = self.seq[i : i + 3]
            if self._translate_codon(codon) == amino_acid:
                if codon in codon_counts:
                    codon_counts[codon] += 1
                else:
                    codon_counts[codon] = 1
                total += 1

        if total > 0:
            for codon in codon_counts:
                codon_counts[codon] /= total

        return codon_counts


class DNASequence(TranslatableSequence):
    type = "DNA"

    def __init__(self, seq: str) -> None:
        super().__init__(seq)

    @classmethod
    def alphabet(cls) -> str:
        return "ACGT"

    @classmethod
    def complement_map(cls) -> dict[str, str]:
        return DNA_COMPLEMENT

    @classmethod
    def codon_map(cls) -> dict[str, str]:
        return DNA_TO_AA

    def complement(self) -> DNASequence:
        """Compute the complementary strand of a given DNA sequence."""
        comp_map = self.complement_map()
        comp = ""
        for base in self.seq:
            comp += comp_map[base]
        return DNASequence(seq=comp)

    def gc_content(self, precision: int = 2) -> float:
        """Returns percentage of G and C nucleotides in a DNA sequence."""
        gc = (self.seq.count("G") + self.seq.count("C")) / len(self.seq)
        return round(gc, precision)

    def transcribe(self) -> RNASequence:
        """Computes the transcribed sequence of a DNA sequence."""
        rna_seq = self.seq.replace("T", "U")
        return RNASequence(seq=rna_seq)

    def reading_frames(self) -> list[ProteinSequence]:
        """Computes the six reading frames of a DNA sequence,
        including reverse complements"""
        rfs = []
        rc = self.reverse().complement()
        for i in range(3):
            rfs.append(self.translate(initial_position=i))
            rfs.append(rc.translate(initial_position=i))
        return rfs


class RNASequence(TranslatableSequence):
    type = "RNA"

    def __init__(self, seq: str) -> None:
        super().__init__(seq)

    @classmethod
    def alphabet(cls) -> str:
        return "ACGU"

    @classmethod
    def complement_map(cls) -> dict[str, str]:
        return RNA_COMPLEMENT

    @classmethod
    def codon_map(cls) -> dict[str, str]:
        return RNA_TO_AA

    def complement(self) -> Self:
        """Compute the complementary strand of a given DNA sequence."""
        comp_map = self.complement_map()
        comp = ""
        for base in self.seq:
            comp += comp_map[base]
        return type(self)(seq=comp)


class ProteinSequence(BaseSequence):
    type = "PROTEIN"

    def __init__(self, seq: str) -> None:
        super().__init__(seq)

    @classmethod
    def alphabet(cls) -> str:
        return "ACDEFGHIKLMNPQRSTVWY"


# class BioSequence:
#     def __init__(self, seq: str, seq_type: str):
