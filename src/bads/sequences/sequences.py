from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Iterator
from typing import Self

from bads.alphabet import (
    DNA_AMBIGUOUS_ALPHABET,
    DNA_UNAMBIGUOUS_ALPHABET,
    RNA_AMBIGUOUS_ALPHABET,
    RNA_UNAMBIGUOUS_ALPHABET,
    PROTEIN_ALPHABET,
)
from bads.sequences.maps import DNA_COMPLEMENT, DNA_TO_AA, RNA_COMPLEMENT, RNA_TO_AA
from bads.exceptions import InvalidCharacterError, TranslationError


class BaseSequence(ABC):
    """Abstract base class for biological sequences.

    This class provides common functionality for all sequence types including
    DNA, RNA, and protein sequences. It handles sequence validation, iteration,
    slicing, and basic sequence operations.

    Attributes:
        seq (str): The sequence string.
    """

    def __init__(self, seq: str) -> None:
        """Initialize a biological sequence.

        Args:
            seq (str): The sequence string to validate and store.

        Raises:
            ValueError: If sequence contains invalid characters for this sequence type.
        """
        self.seq = self.validate_sequence(seq, self.__class__.alphabet())

    def __repr__(self) -> str:
        """Return a detailed string representation of the sequence.

        Returns:
            str: A string showing the class name, length, and abbreviated sequence.
                For sequences < 10 characters, shows the full sequence.
                For longer sequences, shows first 9 and last 3 characters.
        """
        if len(self.seq) < 10:
            seq_repr = self.seq[:9]
        else:
            seq_repr = self.seq[:9] + "..." + self.seq[-3:]
        return f"{self.__class__.__name__}(length={len(self.seq)}, seq={seq_repr})"

    def __len__(self) -> int:
        """Return the length of the sequence.

        Returns:
            int: The number of characters in the sequence.
        """
        return len(self.seq)

    def __hash__(self) -> int:
        """Return a hash of the sequence.

        Returns:
            int: Hash value of the sequence string.
        """
        return hash(self.seq)

    def __iter__(self) -> Iterator:
        """Return an iterator over the sequence characters.

        Returns:
            Iterator: Iterator yielding each character in the sequence.
        """
        return iter(self.seq)

    def __getitem__(self, key: int | slice) -> str:
        """Get a character or subsequence by index or slice.

        Args:
            key (int | slice): Index or slice object.

        Returns:
            str: Single character if key is int, or subsequence if key is slice.
        """
        if isinstance(key, slice):
            s = [self.seq[i] for i in range(*key.indices(len(self.seq)))]
            return "".join(s)
        return self.seq[key]

    def __str__(self) -> str:
        """Return the sequence as a string.

        Returns:
            str: The sequence string.
        """
        return self.seq

    def __contains__(self, s: str) -> bool:
        """Check if a substring is present in the sequence.

        Args:
            s (str): The substring to search for.

        Returns:
            bool: True if substring is found, False otherwise.
        """
        return s in self.seq

    def __add__(self, other: "Sequence") -> "Sequence":
        """Concatenate two sequences of the same types.

        Args:
            other (Sequence): The sequence to concatenate.

        Returns:
            Sequence: The concatenated sequence.
        """
        assert type(self) == type(other), f"Cannot concatenate sequence of different types ({type(self)} + {type(other)})"
        return self.__class__(self.seq + other.seq)

    def __radd__(self, other: "Sequence") -> "Sequence":
        """Concatenate two sequences of the same types.

        Args:
            other (Sequence): The sequence to concatenate.

        Returns:
            Sequence: The concatenated sequence.
        """
        assert type(self) == type(other), f"Cannot concatenate sequence of different types ({type(other)} + {type(self)})"
        return self.__class__(other.seq + self.seq)

    @classmethod
    @abstractmethod
    def alphabet(cls) -> set:
        """Return the valid alphabet for this sequence type.

        Returns:
            set: a collection containing all unique and valid characters for this sequence type.
        """
        ...

    def validate_sequence(self, seq: str, alphabet: set[str]) -> str:
        """Validate that all characters in the sequence are in the alphabet.

        Args:
            seq (str): The sequence to validate.
            alphabet (set): Set of valid characters.

        Returns:
            str: The validated sequence.

        Raises:
            ValueError: If sequence contains characters not in the alphabet.
        """
        if set(seq) <= alphabet:
            return seq
        invalid_chars = set(seq).difference(alphabet)
        raise InvalidCharacterError(
            f"Invalid characters for {self.__class__.__name__}: {invalid_chars}"
        )

    def reverse(self) -> Self:
        """Return the reversed sequence.

        Returns:
            Self: A new sequence object with the sequence reversed.

        Example:
            >>> seq = DNASequence("ACGT")
            >>> seq.reverse().seq
            'TGCA'
        """
        return type(self)(seq=self.seq[::-1])

    def count(self, substr: str) -> int:
        """Count the number of occurrences of a substring in the sequence.

        Args:
            substr (str): The substring to count.

        Returns:
            int: The number of occurrences of the substring.
        """
        return self.seq.count(substr)

    def find(self, substr: str) -> int:
        """Find the index of the first occurrence of a substring in the sequence.

        Args:
            substr (str): The substring to find.

        Returns:
            int: The index of the first occurrence of the substring.
        """
        return self.seq.find(substr)

    def frequency(self) -> dict[str, int]:
        """Calculate the frequency of each character in the sequence.

        Returns:
            dict[str, int]: Dictionary mapping each unique character to its count.

        Example:
            >>> seq = DNASequence("AACGT")
            >>> seq.frequency()
            {'A': 2, 'C': 1, 'G': 1, 'T': 1}
        """
        counts = {}
        for s in self.seq.upper():
            if s in counts:
                counts[s] += 1
            else:
                counts[s] = 1
        return counts


class TranslatableSequence(BaseSequence):
    """Abstract base class for sequences that can be translated to proteins.

    This class extends BaseSequence with methods for translating nucleotide
    sequences (DNA or RNA) into protein sequences using the genetic code.
    """

    @classmethod
    def codon_map(cls) -> dict[str, str]:
        """Return the codon to amino acid mapping for this sequence type.

        Returns:
            dict[str, str]: Dictionary mapping three-letter codons to single-letter
                amino acid codes. Stop codons are represented as "_".

        Raises:
            NotImplementedError: If not implemented by subclass.
        """
        raise NotImplementedError("`codon_map` must be implemented for ", cls.__name__)

    @classmethod
    def _translate_codon(cls, codon: str) -> str:
        """Translate a single codon into its corresponding amino acid.

        Uses the standard genetic code mapping defined in the codon_map.

        Args:
            codon (str): A three-letter codon sequence.

        Returns:
            str: Single-letter amino acid code, or "_" for stop codons.

        Raises:
            AssertionError: If codon is not exactly 3 characters.
            ValueError: If codon is not found in the codon map.

        Example:
            >>> DNASequence._translate_codon("ATG")
            'M'
            >>> DNASequence._translate_codon("TAA")
            '_'
        """
        assert len(codon) == 3, "codon must be a three-character string"
        codon_map = cls.codon_map()
        if codon not in codon_map:
            raise ValueError("Invalid codon")

        return codon_map[codon]

    def translate(
        self, initial_position: int = 0, to_stop: bool = False, cds: bool = False
    ) -> ProteinSequence:
        """Translate the nucleotide sequence into a protein sequence.

        Translates the sequence using the genetic code. Stop codons are never
        included in the resulting protein sequence.

        Args:
            initial_position (int): Position (0-indexed) to start translation.
                Use this to specify reading frames. Default is 0.
            to_stop (bool): If True, stop translation at the first stop codon
                encountered. If False, skip stop codons and continue translating
                until the end of the sequence. Default is False.
            cds (bool): If True, validate the sequence as a coding sequence (CDS).
                A valid CDS must start with a start codon (M) and contain at least
                one stop codon. Default is False.

        Returns:
            ProteinSequence: The translated protein sequence. Stop codons are
                never included in the output.

        Raises:
            TranslationError: If cds=True and the sequence doesn't start with
                a start codon or doesn't contain a stop codon.

        Example:
            >>> dna = DNASequence("ATGGCTTAA")
            >>> dna.translate().seq
            'MA'
            >>> dna.translate(to_stop=True).seq
            'MA'
            >>> dna.translate(cds=True).seq
            'MA'
            >>> dna.translate(initial_position=1).seq
            'WL'
        """
        # Translate to full amino acid sequence
        aa_seq = ""
        for pos in range(initial_position, len(self.seq) - 2, 3):
            codon = self.seq[pos : pos + 3]
            amino_acid = self._translate_codon(codon)
            if amino_acid == "_":
                if to_stop or cds:
                    break
                # Skip stop codons in normal translation
                continue
            aa_seq += amino_acid

        if cds:
            if not aa_seq or aa_seq[0] != "M":
                raise TranslationError(
                    f"First triplet {self.seq[initial_position : initial_position + 3]} is not a start codon"
                )
            # In CDS mode, we already stopped at the first stop codon
            # Just verify we actually found one
            remaining_pos = initial_position + (len(aa_seq) + 1) * 3
            if remaining_pos > len(self.seq):
                raise TranslationError(
                    f"Amino acid sequence does not contain a stop codon"
                )

        return ProteinSequence(seq=aa_seq)

    def find_proteins(self, min_length: int = 0):
        # Translate sequence including stop codons
        seq_aa = ""
        for pos in range(0, len(self.seq) - 2, 3):
            codon = self.seq[pos : pos + 3]
            seq_aa += self._translate_codon(codon)

        proteins = []
        current_prot = []
        for aa in seq_aa:
            if aa == "_":
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                if aa == "M":
                    current_prot.append("")  # Start with empty string
                for i in range(len(current_prot)):
                    current_prot[i] += aa  # Add aa to all proteins, including new ones

        filtered_proteins = []
        for prot in proteins:
            if len(prot) >= min_length:
                filtered_proteins.append(ProteinSequence(seq=prot))
        return sorted(filtered_proteins, key=lambda p: len(p), reverse=True)

    def codon_usage(self, amino_acid: str) -> dict[str, int]:
        """Calculate the codon usage frequencies for a specific amino acid.

        Analyzes the sequence in triplets (starting from position 0) and calculates
        the normalized frequency of each codon that encodes the specified amino acid.

        Args:
            amino_acid (str): Single-letter amino acid code to analyze.

        Returns:
            dict[str, int]: Dictionary mapping codons to their normalized frequencies
                (values sum to 1.0). Empty dict if amino acid is not found.

        Example:
            >>> dna = DNASequence("CTTCTCCTACTG")
            >>> dna.codon_usage("L")
            {'CTT': 0.25, 'CTC': 0.25, 'CTA': 0.25, 'CTG': 0.25}
        """
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
    """Represents a DNA sequence.

    DNA sequences use the alphabet A, C, G, T. This class provides methods
    for DNA-specific operations including complementation, transcription to RNA,
    translation to protein, and reading frame analysis.

    Attributes:
        type (str): Sequence type identifier ("DNA").
        seq (str): The DNA sequence string.

    Example:
        >>> dna = DNASequence("ATGGCTTAA")
        >>> dna.complement().seq
        'TACCGAATT'
        >>> dna.transcribe().seq
        'AUGGCUUAA'
        >>> dna.translate().seq
        'MA'
    """

    type = "DNA"

    def __init__(self, seq: str) -> None:
        """Initialize a DNA sequence.

        Args:
            seq (str): DNA sequence string containing only A, C, G, T.

        Raises:
            ValueError: If sequence contains invalid characters.
        """
        super().__init__(seq)

    @classmethod
    def alphabet(cls, ambiguous: bool = True) -> set[str]:
        """Return the DNA alphabet.

        Returns:
            set[str]
        """
        if ambiguous:
            return DNA_AMBIGUOUS_ALPHABET
        return DNA_UNAMBIGUOUS_ALPHABET

    @classmethod
    def complement_map(cls) -> dict[str, str]:
        """Return the DNA base pairing map.

        Returns:
            dict[str, str]: Dictionary mapping each base to its complement.
                A↔T, G↔C
        """
        return DNA_COMPLEMENT

    @classmethod
    def codon_map(cls) -> dict[str, str]:
        """Return the DNA codon to amino acid mapping.

        Returns:
            dict[str, str]: Dictionary mapping DNA codons to amino acids
                using the standard genetic code.
        """
        return DNA_TO_AA

    def complement(self) -> DNASequence:
        """Compute the complementary DNA strand.

        Generates the complement by applying Watson-Crick base pairing rules:
        A↔T, G↔C.

        Returns:
            DNASequence: The complementary DNA sequence.

        Example:
            >>> dna = DNASequence("ACGT")
            >>> dna.complement().seq
            'TGCA'
        """
        comp_map = self.complement_map()
        comp = ""
        for base in self.seq:
            comp += comp_map[base]
        return DNASequence(seq=comp)

    def gc_content(self, precision: int = 2) -> float:
        """Calculate the GC content of the DNA sequence.

        GC content is the percentage of bases that are either guanine (G)
        or cytosine (C). This is an important metric in molecular biology
        as it affects DNA stability and melting temperature.

        Args:
            precision (int): Number of decimal places to round to. Default is 2.

        Returns:
            float: GC content as a value between 0 and 1, rounded to the
                specified precision.

        Example:
            >>> dna = DNASequence("GGCC")
            >>> dna.gc_content()
            1.0
            >>> dna = DNASequence("ACGT")
            >>> dna.gc_content()
            0.5
        """
        gc = (self.seq.count("G") + self.seq.count("C")) / len(self.seq)
        return round(gc, precision)

    def transcribe(self) -> RNASequence:
        """Transcribe the DNA sequence to RNA.

        Performs transcription by replacing all thymine (T) bases with
        uracil (U) bases, converting DNA to RNA.

        Returns:
            RNASequence: The transcribed RNA sequence.

        Example:
            >>> dna = DNASequence("ATGGCT")
            >>> dna.transcribe().seq
            'AUGGCU'
        """
        rna_seq = self.seq.replace("T", "U")
        return RNASequence(seq=rna_seq)

    def reading_frames(self) -> list[ProteinSequence]:
        """Compute all six reading frames of the DNA sequence.

        Generates all possible reading frames: three on the forward strand
        (starting at positions 0, 1, 2) and three on the reverse complement
        strand (also starting at positions 0, 1, 2).

        Returns:
            list[ProteinSequence]: List of 6 protein sequences, one for each
                reading frame. Order is: forward frame 0, reverse frame 0,
                forward frame 1, reverse frame 1, forward frame 2, reverse frame 2.

        Example:
            >>> dna = DNASequence("ATGGCTTAA")
            >>> frames = dna.reading_frames()
            >>> len(frames)
            6
            >>> frames[0].seq  # Forward frame 0
            'MA'
        """
        rfs = []
        rc = self.reverse().complement()
        for i in range(3):
            rfs.append(self.translate(initial_position=i))
            rfs.append(rc.translate(initial_position=i))
        return rfs


class RNASequence(TranslatableSequence):
    """Represents an RNA sequence.

    RNA sequences use the alphabet A, C, G, U. This class provides methods
    for RNA-specific operations including complementation and translation to protein.

    Attributes:
        type (str): Sequence type identifier ("RNA").
        seq (str): The RNA sequence string.

    Example:
        >>> rna = RNASequence("AUGGCUUAA")
        >>> rna.complement().seq
        'UACCGAAUU'
        >>> rna.translate().seq
        'MA'
    """

    type = "RNA"

    def __init__(self, seq: str) -> None:
        """Initialize an RNA sequence.

        Args:
            seq (str): RNA sequence string containing only A, C, G, U.

        Raises:
            ValueError: If sequence contains invalid characters.
        """
        super().__init__(seq)

    @classmethod
    def alphabet(cls, ambiguous: bool = True) -> set[str]:
        """Return the RNA alphabet.

        Returns:
            str: "ACGU"
        """
        if ambiguous:
            return RNA_AMBIGUOUS_ALPHABET
        return RNA_UNAMBIGUOUS_ALPHABET

    @classmethod
    def complement_map(cls) -> dict[str, str]:
        """Return the RNA base pairing map.

        Returns:
            dict[str, str]: Dictionary mapping each base to its complement.
                A↔U, G↔C
        """
        return RNA_COMPLEMENT

    @classmethod
    def codon_map(cls) -> dict[str, str]:
        """Return the RNA codon to amino acid mapping.

        Returns:
            dict[str, str]: Dictionary mapping RNA codons to amino acids
                using the standard genetic code.
        """
        return RNA_TO_AA

    def complement(self) -> Self:
        """Compute the complementary RNA strand.

        Generates the complement by applying Watson-Crick base pairing rules:
        A↔U, G↔C.

        Returns:
            RNASequence: The complementary RNA sequence.

        Example:
            >>> rna = RNASequence("ACGU")
            >>> rna.complement().seq
            'UGCA'
        """
        comp_map = self.complement_map()
        comp = ""
        for base in self.seq:
            comp += comp_map[base]
        return type(self)(seq=comp)


class ProteinSequence(BaseSequence):
    """Represents a protein (amino acid) sequence.

    Protein sequences use the standard 20 amino acid single-letter codes.
    This class provides basic sequence operations inherited from BaseSequence.

    Attributes:
        type (str): Sequence type identifier ("PROTEIN").
        seq (str): The protein sequence string.

    Example:
        >>> protein = ProteinSequence("MAPLK")
        >>> len(protein)
        5
        >>> protein.frequency()
        {'M': 1, 'A': 1, 'P': 1, 'L': 1, 'K': 1}
    """

    type = "PROTEIN"

    def __init__(self, seq: str) -> None:
        """Initialize a protein sequence.

        Args:
            seq (str): Protein sequence using single-letter amino acid codes.
                Valid characters: A, C, D, E, F, G, H, I, K, L, M, N, P, Q,
                R, S, T, V, W, Y.

        Raises:
            ValueError: If sequence contains invalid characters.
        """
        super().__init__(seq)

    @classmethod
    def alphabet(cls) -> set[str]:
        """Return the protein alphabet.

        Returns:
            str: "ACDEFGHIKLMNPQRSTVWY" (20 standard amino acids)
        """
        return PROTEIN_ALPHABET


# class BioSequence:
#     def __init__(self, seq: str, seq_type: str):
