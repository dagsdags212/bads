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
