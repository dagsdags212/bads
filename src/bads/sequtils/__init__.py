def search_first_occurrence(seq: str, pattern: str) -> int:
    """Search for the first occurrence of a pattern in a sequence.

    Args:
        seq (str): The sequence to search in.
        pattern (str): The pattern to search for.

    Returns:
        int: The index of the first occurrence of the pattern in the sequence.
    """
    found = False
    i = 0
    while i <= len(seq) - len(pattern) and not found:
        j = 0
        while j < len(pattern) and pattern[j] == seq[i+j]:
            j += 1
        if j == len(pattern):
            found = True
        else:
            i += 1
    if found:
        return i
    return -1

def search_all_occurrences(seq: str, pattern: str) -> list[int]:
    """Search for all occurrences of a pattern in a sequence.

    Args:
        seq (str): The sequence to search in.
        pattern (str): The pattern to search for.

    Returns:
        list[int]: A list of indices where the pattern is found in the sequence.
    """
    res = []
    for i in range(len(seq) - len(pattern) + 1):
        j = 0
        while j < len(pattern) and pattern[j] == seq[i+j]:
            j += 1
        if j == len(pattern):
            res.append(i)
    return res


class BoyerMoore:

    def __init__(self, alphabet: str, pattern: str):
        self.alphabet = alphabet
        self.pattern = pattern
        self.preprocess()

    def preprocess(self) -> None:
        self.process_bcr()
        self.process_gsr()

    def process_bcr(self) -> None:
        """Preprocessing step for the Bad Character Rule."""
        self.occ = {}
        for symbol in self.alphabet:
            self.occ[symbol] = -1
        for j in range(len(self.pattern)):
            c = self.pattern[j]
            self.occ[c] = j

    def process_gsr(self) -> None:
        """Preprocessing step for the Good Suffix Rule."""
        m = len(self.pattern)
        self.f = [0] * (m + 1)
        self.s = [m] * (m + 1)

        # Case 1: Suffix exists elsewhere in pattern
        i, j = m, m + 1
        self.f[i] = j
        while i > 0:
            while j <= m and self.pattern[i-1] != self.pattern[j-1]:
                if self.s[j] == m:
                    self.s[j] = j-i
                j = self.f[j]
            i -= 1
            j -= 1
            self.f[i] = j

        # Case 2: A portion of the suffix is also a prefix
        j = self.f[0]
        for i in range(m + 1):
            if self.s[i] == m:
                self.s[i] = j
            if i == j:
                j = self.f[j]

    def search(self, text: str) -> list[int]:
        n, m = len(text), len(self.pattern)
        if m == 0:
            return []

        res = []
        shift = 0
        while shift <= (n - m):
            j = m - 1
            while j >= 0 and self.pattern[j] == text[shift + j]:
                j -= 1
            if j < 0:
                res.append(shift)
                shift += self.s[0]
            else:
                bad_char_val = self.occ.get(text[shift + j], -1)
                bc_shift = j - bad_char_val
                gs_shift = self.s[j + 1]
                shift += max(bc_shift, gs_shift)
        return res
