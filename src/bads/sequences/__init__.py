from __future__ import annotations
from abc import ABC, abstractmethod
from collections.abc import Iterator


class BaseSequence(ABC):
    def __init__(self, seq: str) -> None:
        self.seq = self.validate_sequence(seq)

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

    @abstractmethod
    def validate_sequence(self, seq: str) -> str:
        raise NotImplementedError

    def frequency(self) -> dict[str, int]:
        """Returns counts of unique characters in a given sequence."""
        counts = {}
        for s in self.seq.upper():
            if s in counts:
                counts[s] += 1
            else:
                counts[s] = 1
        return counts


class DNASequence(BaseSequence):
    type = "DNA"

    def __init__(self, seq: str) -> None:
        super().__init__(seq)

    def validate_sequence(self, seq: str) -> str:
        """Checks if DNA sequence is valid. Returns True if sequence is valid,
        or False otherwise."""
        s = seq.upper()
        valid_len = s.count("A") + s.count("T") + s.count("G") + s.count("C")
        if valid_len != len(s):
            raise ValueError("Invalid DNA string")
        return s

    def gc_content(self, precision: int = 2) -> float:
        """Returns percentage of G and C nucleotides in a DNA sequence."""
        gc = (self.seq.count("G") + self.seq.count("C")) / len(self.seq)
        return round(gc, precision)
