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


def _nw_score_row(
    seq1: str,
    seq2: str,
    match_score: int,
    mismatch_score: int,
    gap_penalty: int,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> list[int]:
    """
    Return only the final DP row for global alignment of seq1 vs seq2.
    Uses O(len(seq2)) space.
    """
    prev = [j * gap_penalty for j in range(len(seq2) + 1)]

    for i in range(1, len(seq1) + 1):
        curr = [i * gap_penalty] + [0] * len(seq2)
        for j in range(1, len(seq2) + 1):
            diag = prev[j - 1] + _score_pair(
                seq1[i - 1],
                seq2[j - 1],
                match_score,
                mismatch_score,
                substitution_matrix,
            )
            up = prev[j] + gap_penalty
            left = curr[j - 1] + gap_penalty
            curr[j] = max(diag, up, left)
        prev = curr

    return prev


def _needleman_wunsch_small(
    seq1: str,
    seq2: str,
    match_score: int,
    mismatch_score: int,
    gap_penalty: int,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> tuple[str, str]:
    """
    Standard Needleman-Wunsch used as base case in Hirschberg recursion.
    """
    n = len(seq1)
    m = len(seq2)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]
    trace[0][0] = "E"

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap_penalty
        trace[i][0] = "U"

    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap_penalty
        trace[0][j] = "L"

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

            if best == diag:
                trace[i][j] = "D"
            elif best == up:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

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
    return "".join(aligned1), "".join(aligned2)


def _hirschberg_recursive(
    seq1: str,
    seq2: str,
    match_score: int,
    mismatch_score: int,
    gap_penalty: int,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> tuple[str, str]:
    """
    Hirschberg divide-and-conquer global alignment.
    """
    n = len(seq1)
    m = len(seq2)

    if n == 0:
        return "-" * m, seq2
    if m == 0:
        return seq1, "-" * n
    if n == 1 or m == 1:
        return _needleman_wunsch_small(
            seq1,
            seq2,
            match_score,
            mismatch_score,
            gap_penalty,
            substitution_matrix,
        )

    mid = n // 2

    score_l = _nw_score_row(
        seq1[:mid],
        seq2,
        match_score,
        mismatch_score,
        gap_penalty,
        substitution_matrix,
    )
    score_r = _nw_score_row(
        seq1[mid:][::-1],
        seq2[::-1],
        match_score,
        mismatch_score,
        gap_penalty,
        substitution_matrix,
    )

    split_j = max(
        range(m + 1),
        key=lambda j: score_l[j] + score_r[m - j],
    )

    left_a1, left_a2 = _hirschberg_recursive(
        seq1[:mid],
        seq2[:split_j],
        match_score,
        mismatch_score,
        gap_penalty,
        substitution_matrix,
    )
    right_a1, right_a2 = _hirschberg_recursive(
        seq1[mid:],
        seq2[split_j:],
        match_score,
        mismatch_score,
        gap_penalty,
        substitution_matrix,
    )

    return left_a1 + right_a1, left_a2 + right_a2


def _score_alignment(
    aln1: str,
    aln2: str,
    match_score: int,
    mismatch_score: int,
    gap_penalty: int,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> int:
    total = 0
    for a, b in zip(aln1, aln2):
        if a == "-" or b == "-":
            total += gap_penalty
        else:
            total += _score_pair(a, b, match_score, mismatch_score, substitution_matrix)
    return total


def hirschberg(
    seq1: str,
    seq2: str,
    *,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_penalty: int = -1,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> AlignmentResult:
    """
    Global alignment using Hirschberg's linear-space divide-and-conquer algorithm.
    """
    if gap_penalty > 0:
        raise ValueError("gap_penalty should usually be zero or negative for alignment.")
    if any(ch == "-" for ch in seq1) or any(ch == "-" for ch in seq2):
        raise ValueError("Input sequences must not already contain gap characters.")

    aln1, aln2 = _hirschberg_recursive(
        seq1,
        seq2,
        match_score,
        mismatch_score,
        gap_penalty,
        substitution_matrix,
    )

    score = _score_alignment(
        aln1,
        aln2,
        match_score,
        mismatch_score,
        gap_penalty,
        substitution_matrix,
    )

    return AlignmentResult(
        score=score,
        aligned_seq1=aln1,
        aligned_seq2=aln2,
        start1=0,
        end1=len(seq1),
        start2=0,
        end2=len(seq2),
    )
