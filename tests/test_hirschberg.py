from src.pairwise.hirschberg import hirschberg
from src.pairwise.needleman_wunsch import needleman_wunsch


def test_identical_sequences_align_perfectly():
    result = hirschberg("ACGT", "ACGT", match_score=2, mismatch_score=-1, gap_penalty=-2)

    assert result.score == 8
    assert result.aligned_seq1 == "ACGT"
    assert result.aligned_seq2 == "ACGT"
    assert result.matches == 4
    assert result.mismatches == 0
    assert result.gaps == 0


def test_empty_sequence_aligns_with_gaps():
    result = hirschberg("", "ACG", match_score=1, mismatch_score=-1, gap_penalty=-2)

    assert result.score == -6
    assert result.aligned_seq1 == "---"
    assert result.aligned_seq2 == "ACG"


def test_hirschberg_matches_needleman_wunsch_score():
    seq1 = "GATTACA"
    seq2 = "GCATGCU"

    h = hirschberg(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1)
    nw = needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1)

    assert h.score == nw.score
    assert h.aligned_seq1.replace("-", "") == seq1
    assert h.aligned_seq2.replace("-", "") == seq2
    assert len(h.aligned_seq1) == len(h.aligned_seq2)


def test_custom_substitution_matrix_is_used():
    matrix = {
        ("A", "A"): 3,
        ("C", "C"): 3,
        ("A", "C"): -2,
    }

    result = hirschberg("A", "C", gap_penalty=-4, substitution_matrix=matrix)

    assert result.score == -2
    assert result.aligned_seq1 == "A"
    assert result.aligned_seq2 == "C"


def test_rejects_input_sequences_with_existing_gap_characters():
    try:
        hirschberg("A-C", "AC")
        assert False, "Expected ValueError for pre-gapped input."
    except ValueError as exc:
        assert "must not already contain gap characters" in str(exc)