from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple


NEG_INF = -10**12


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


def gotoh(
    seq1: str,
    seq2: str,
    *,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_open: int = -2,
    gap_extend: int = -1,
    substitution_matrix: Optional[Dict[Tuple[str, str], int]] = None,
) -> AlignmentResult:
    """
    Global pairwise alignment using Gotoh's affine-gap dynamic programming.

    States:
    M[i][j] = best score ending in match/mismatch
    X[i][j] = best score ending with a gap in seq2 (vertical move)
    Y[i][j] = best score ending with a gap in seq1 (horizontal move)
    """
    if gap_open > 0 or gap_extend > 0:
        raise ValueError("Gap penalties should usually be zero or negative.")
    if any(ch == "-" for ch in seq1) or any(ch == "-" for ch in seq2):
        raise ValueError("Input sequences must not already contain gap characters.")

    n = len(seq1)
    m = len(seq2)

    M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

    trace_M = [[None] * (m + 1) for _ in range(n + 1)]
    trace_X = [[None] * (m + 1) for _ in range(n + 1)]
    trace_Y = [[None] * (m + 1) for _ in range(n + 1)]

    M[0][0] = 0

    for i in range(1, n + 1):
        X[i][0] = gap_open + (i - 1) * gap_extend
        trace_X[i][0] = "X" if i > 1 else "M"

    for j in range(1, m + 1):
        Y[0][j] = gap_open + (j - 1) * gap_extend
        trace_Y[0][j] = "Y" if j > 1 else "M"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = _score_pair(
                seq1[i - 1],
                seq2[j - 1],
                match_score,
                mismatch_score,
                substitution_matrix,
            )

            # M state
            candidates_M = [
                (M[i - 1][j - 1] + s, "M"),
                (X[i - 1][j - 1] + s, "X"),
                (Y[i - 1][j - 1] + s, "Y"),
            ]
            M[i][j], trace_M[i][j] = max(candidates_M, key=lambda x: x[0])

            # X state: gap in seq2
            candidates_X = [
                (M[i - 1][j] + gap_open, "M"),
                (X[i - 1][j] + gap_extend, "X"),
            ]
            X[i][j], trace_X[i][j] = max(candidates_X, key=lambda x: x[0])

            # Y state: gap in seq1
            candidates_Y = [
                (M[i][j - 1] + gap_open, "M"),
                (Y[i][j - 1] + gap_extend, "Y"),
            ]
            Y[i][j], trace_Y[i][j] = max(candidates_Y, key=lambda x: x[0])

    final_candidates = [
        (M[n][m], "M"),
        (X[n][m], "X"),
        (Y[n][m], "Y"),
    ]
    final_score, state = max(final_candidates, key=lambda x: x[0])

    aligned1 = []
    aligned2 = []
    i, j = n, m

    while i > 0 or j > 0:
        if state == "M":
            prev_state = trace_M[i][j]
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1
            state = prev_state
        elif state == "X":
            prev_state = trace_X[i][j]
            aligned1.append(seq1[i - 1])
            aligned2.append("-")
            i -= 1
            state = prev_state
        elif state == "Y":
            prev_state = trace_Y[i][j]
            aligned1.append("-")
            aligned2.append(seq2[j - 1])
            j -= 1
            state = prev_state
        else:
            raise RuntimeError("Invalid traceback state.")

    aligned1.reverse()
    aligned2.reverse()

    return AlignmentResult(
        score=final_score,
        aligned_seq1="".join(aligned1),
        aligned_seq2="".join(aligned2),
        start1=0,
        end1=n,
        start2=0,
        end2=m,
    )