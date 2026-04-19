import pytest

from src.pairwise.needleman_wunsch import needleman_wunsch


def test_identical_sequences_align_perfectly():
    result = needleman_wunsch("ACGT", "ACGT", match_score=2, mismatch_score=-1, gap_penalty=-2)

    assert result.score == 8
    assert result.aligned_seq1 == "ACGT"
    assert result.aligned_seq2 == "ACGT"
    assert result.matches == 4
    assert result.mismatches == 0
    assert result.gaps == 0


def test_empty_sequence_aligns_with_gaps():
    result = needleman_wunsch("", "ACG", match_score=1, mismatch_score=-1, gap_penalty=-2)

    assert result.score == -6
    assert result.aligned_seq1 == "---"
    assert result.aligned_seq2 == "ACG"


def test_alignment_preserves_original_sequences():
    seq1 = "GATTACA"
    seq2 = "GCATGCU"

    result = needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1)

    assert result.aligned_seq1.replace("-", "") == seq1
    assert result.aligned_seq2.replace("-", "") == seq2
    assert len(result.aligned_seq1) == len(result.aligned_seq2)


def test_custom_substitution_matrix_is_used():
    matrix = {
        ("A", "A"): 3,
        ("C", "C"): 3,
        ("A", "C"): -2,
    }

    result = needleman_wunsch(
        "A",
        "C",
        gap_penalty=-4,
        substitution_matrix=matrix,
    )

    assert result.score == -2
    assert result.aligned_seq1 == "A"
    assert result.aligned_seq2 == "C"


def test_rejects_input_sequences_with_existing_gap_characters():
    try:
        needleman_wunsch("A-C", "AC")
        assert False, "Expected ValueError for pre-gapped input."
    except ValueError as exc:
        assert "must not already contain gap characters" in str(exc)


def test_rejects_positive_gap_penalty():
    with pytest.raises(ValueError, match="gap_penalty"):
        needleman_wunsch("AC", "AC", gap_penalty=1)


def test_substitution_matrix_supports_reverse_pair_lookup():
    matrix = {
        ("A", "C"): -2,
    }

    result = needleman_wunsch("C", "A", gap_penalty=-4, substitution_matrix=matrix)

    assert result.score == -2
    assert result.aligned_seq1 == "C"
    assert result.aligned_seq2 == "A"


def test_missing_substitution_matrix_entry_raises_keyerror():
    matrix = {
        ("A", "A"): 3,
    }

    with pytest.raises(KeyError, match="Pair \\(A, C\\) not found"):
        needleman_wunsch("A", "C", substitution_matrix=matrix)


def test_alignment_statistics_track_gaps_and_matches():
    result = needleman_wunsch("AC", "A", match_score=2, mismatch_score=-1, gap_penalty=-2)

    assert result.score == 0
    assert result.aligned_seq1 == "AC"
    assert result.aligned_seq2 == "A-"
    assert result.matches == 1
    assert result.mismatches == 0
    assert result.gaps == 1
