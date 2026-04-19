from src.pairwise.smith_waterman import smith_waterman


def test_identical_sequences_align_perfectly():
    result = smith_waterman("ACGT", "ACGT", match_score=2, mismatch_score=-1, gap_penalty=-2)

    assert result.score == 8
    assert result.aligned_seq1 == "ACGT"
    assert result.aligned_seq2 == "ACGT"
    assert result.matches == 4
    assert result.mismatches == 0
    assert result.gaps == 0


def test_local_alignment_finds_internal_match():
    result = smith_waterman("TTACGTAA", "GGACGTCC", match_score=2, mismatch_score=-1, gap_penalty=-2)

    assert result.score == 8
    assert result.aligned_seq1 == "ACGT"
    assert result.aligned_seq2 == "ACGT"
    assert result.start1 == 2
    assert result.end1 == 6
    assert result.start2 == 2
    assert result.end2 == 6


def test_no_similarity_returns_zero_score():
    result = smith_waterman("AAAA", "TTTT", match_score=2, mismatch_score=-1, gap_penalty=-2)

    assert result.score == 0
    assert result.aligned_seq1 == ""
    assert result.aligned_seq2 == ""


def test_rejects_input_sequences_with_existing_gap_characters():
    try:
        smith_waterman("A-C", "AC")
        assert False, "Expected ValueError for pre-gapped input."
    except ValueError as exc:
        assert "must not already contain gap characters" in str(exc)


def test_custom_substitution_matrix_is_used():
    matrix = {
        ("A", "A"): 3,
        ("C", "C"): 3,
        ("A", "C"): -2,
    }

    result = smith_waterman(
        "TAAC",
        "GAAC",
        gap_penalty=-4,
        substitution_matrix=matrix,
    )

    assert result.score == 6
    assert result.aligned_seq1 == "AC"
    assert result.aligned_seq2 == "AC"