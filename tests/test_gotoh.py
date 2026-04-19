from src.pairwise.gotoh import gotoh


def test_identical_sequences_align_perfectly():
    result = gotoh("ACGT", "ACGT", match_score=2, mismatch_score=-1, gap_open=-3, gap_extend=-1)

    assert result.score == 8
    assert result.aligned_seq1 == "ACGT"
    assert result.aligned_seq2 == "ACGT"
    assert result.matches == 4
    assert result.mismatches == 0
    assert result.gaps == 0


def test_empty_sequence_aligns_with_one_affine_gap():
    result = gotoh("", "ACG", match_score=1, mismatch_score=-1, gap_open=-2, gap_extend=-1)

    # one gap of length 3 = -2 + 2*(-1) = -4
    assert result.score == -4
    assert result.aligned_seq1 == "---"
    assert result.aligned_seq2 == "ACG"


def test_affine_prefers_one_long_gap():
    result = gotoh("AAAAAA", "AAA", match_score=1, mismatch_score=-1, gap_open=-2, gap_extend=-1)

    assert result.matches == 3
    assert result.gaps == 3
    # should be 3 matches + one gap length 3 => 3 + (-2) + 2*(-1) = -1
    assert result.score == -1


def test_rejects_input_sequences_with_existing_gap_characters():
    try:
        gotoh("A-C", "AC")
        assert False, "Expected ValueError for pre-gapped input."
    except ValueError as exc:
        assert "must not already contain gap characters" in str(exc)