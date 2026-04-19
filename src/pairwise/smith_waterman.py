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


def smith_waterman(
    seq1: str,
    seq2: str,
    *,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -1,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> AlignmentResult:
    """
    Local pairwise alignment using the Smith-Waterman algorithm
    with a linear gap penalty.
    """
    if gap_penalty > 0:
        raise ValueError("gap_penalty should usually be zero or negative for alignment.")
    if any(ch == "-" for ch in seq1) or any(ch == "-" for ch in seq2):
        raise ValueError("Input sequences must not already contain gap characters.")

    n = len(seq1)
    m = len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_pos = (0, 0)

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

            best = max(0, diag, up, left)
            score[i][j] = best

            if best == 0:
                trace[i][j] = "E"
            elif best == diag:
                trace[i][j] = "D"
            elif best == up:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

            if best > max_score:
                max_score = best
                max_pos = (i, j)

    aligned1 = []
    aligned2 = []

    i, j = max_pos
    end1 = i
    end2 = j

    while i > 0 and j > 0 and score[i][j] > 0:
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
            break

    start1 = i
    start2 = j

    aligned1.reverse()
    aligned2.reverse()

    return AlignmentResult(
        score=max_score,
        aligned_seq1="".join(aligned1),
        aligned_seq2="".join(aligned2),
        start1=start1,
        end1=end1,
        start2=start2,
        end2=end2,
    )