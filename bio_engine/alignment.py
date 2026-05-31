"""Needleman-Wunsch global alignment algorithm for DNA sequences."""

from typing import Tuple


def needleman_wunsch(
    seq1: str,
    seq2: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> Tuple[str, str, float]:
    """Align two DNA sequences globally and return the alignment and final score.

    Args:
        seq1: First DNA sequence (typically wild-type/reference).
        seq2: Second DNA sequence (typically mutated/query).
        match: Score for match.
        mismatch: Penalty for mismatch.
        gap: Penalty for gap.

    Returns:
        A tuple of (aligned_seq1, aligned_seq2, score).
    """
    m, n = len(seq1), len(seq2)
    # Initialize DP matrix with dimensions (m+1) x (n+1)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    # Fill base cases (first row and column)
    for i in range(1, m + 1):
        score_matrix[i][0] = i * gap
    for j in range(1, n + 1):
        score_matrix[0][j] = j * gap

    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = score_matrix[i - 1][j - 1] + match_score
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(diag, up, left)

    # Traceback to reconstruct alignment
    aligned1 = []
    aligned2 = []
    i, j = m, n

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            match_score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            if score_matrix[i][j] == score_matrix[i - 1][j - 1] + match_score:
                aligned1.append(seq1[i - 1])
                aligned2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0:
            if score_matrix[i][j] == score_matrix[i - 1][j] + gap:
                aligned1.append(seq1[i - 1])
                aligned2.append("-")
                i -= 1
                continue
        if j > 0:
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1

    # Reverse alignments since we traced backwards
    aligned_seq1 = "".join(reversed(aligned1))
    aligned_seq2 = "".join(reversed(aligned2))
    score = score_matrix[m][n]

    return aligned_seq1, aligned_seq2, float(score)


def calculate_identity_and_distance(
    aligned_seq1: str, aligned_seq2: str
) -> Tuple[float, float]:
    """Calculate sequence identity and evolutionary distance from aligned sequences.

    Sequence identity = matches / alignment_length
    Distance = 1.0 - identity

    Args:
        aligned_seq1: Aligned reference sequence (with gaps).
        aligned_seq2: Aligned query sequence (with gaps).

    Returns:
        A tuple of (identity_percentage, distance_score).
    """
    if not aligned_seq1 or not aligned_seq2 or len(aligned_seq1) != len(aligned_seq2):
        return 0.0, 1.0

    alignment_len = len(aligned_seq1)
    matches = sum(
        1
        for a, b in zip(aligned_seq1, aligned_seq2)
        if a == b and a != "-" and b != "-"
    )

    identity = matches / alignment_len
    distance = 1.0 - identity

    return identity, distance
