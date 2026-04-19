from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple


@dataclass
class AlignmentResult:
    score: int
    aligned_seq1: str
    aligned_seq2: str
    start1: int
    end1: int
    start2: int
    end2: int

    @property
    def matches(self) -> int:
        return sum(
            1 for a, b in zip(self.aligned_seq1, self.aligned_seq2)
            if a == b and a != "-"
        )

    @property
    def mismatches(self) -> int:
        return sum(
            1 for a, b in zip(self.aligned_seq1, self.aligned_seq2)
            if a != "-" and b != "-" and a != b
        )

    @property
    def gaps(self) -> int:
        return sum(
            1 for a, b in zip(self.aligned_seq1, self.aligned_seq2)
            if a == "-" or b == "-"
        )


def _score_pair(
    a: str,
    b: str,
    match_score: int,
    mismatch_score: int,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> int:
    if substitution_matrix is not None:
        if (a, b) in substitution_matrix:
            return substitution_matrix[(a, b)]
        if (b, a) in substitution_matrix:
            return substitution_matrix[(b, a)]
        raise KeyError(f"Pair ({a}, {b}) not found in substitution matrix.")
    return match_score if a == b else mismatch_score


def needleman_wunsch(
    seq1: str,
    seq2: str,
    *,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -1,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> AlignmentResult:
    """
    Global pairwise alignment using the Needleman-Wunsch algorithm
    with a linear gap penalty.

    Parameters
    ----------
    seq1, seq2 : str
        Input sequences.
    match_score : int
        Score for a match when no substitution matrix is supplied.
    mismatch_score : int
        Score for a mismatch when no substitution matrix is supplied.
    gap_penalty : int
        Linear penalty for inserting a gap.
    substitution_matrix : dict[(char, char), int] | None
        Optional substitution matrix. If supplied, it overrides
        match_score and mismatch_score.

    Returns
    -------
    AlignmentResult
        Score and reconstructed global alignment.
    """
    if gap_penalty > 0:
        raise ValueError("gap_penalty should usually be zero or negative for alignment.")
    if any(ch == "-" for ch in seq1) or any(ch == "-" for ch in seq2):
        raise ValueError("Input sequences must not already contain gap characters.")

    n = len(seq1)
    m = len(seq2)

    # DP score matrix
    score = [[0] * (m + 1) for _ in range(n + 1)]

    # Traceback matrix:
    # 'D' = diagonal, 'U' = up, 'L' = left, 'E' = entry/start
    trace = [[""] * (m + 1) for _ in range(n + 1)]
    trace[0][0] = "E"

    # Initialization
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap_penalty
        trace[i][0] = "U"

    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap_penalty
        trace[0][j] = "L"

    # Fill
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + _score_pair(
                seq1[i - 1],
                seq2[j - 1],
                match_score,
                mismatch_score,
                substitution_matrix,
            )
            up = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty

            best = max(diag, up, left)
            score[i][j] = best

            # Deterministic tie-breaking:
            # prefer diagonal, then up, then left
            if best == diag:
                trace[i][j] = "D"
            elif best == up:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

    # Traceback from bottom-right
    aligned1 = []
    aligned2 = []
    i, j = n, m

    while i > 0 or j > 0:
        direction = trace[i][j]

        if direction == "D":
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == "U":
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
        elif direction == "L":
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1
        else:
            raise RuntimeError(f"Invalid traceback direction at ({i}, {j}): {direction}")

    aligned1.reverse()
    aligned2.reverse()

    return AlignmentResult(
        score=score[n][m],
        aligned_seq1="".join(aligned1),
        aligned_seq2="".join(aligned2),
        start1=0,
        end1=n,
        start2=0,
        end2=m,
    )