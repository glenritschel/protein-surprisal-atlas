import pytest
from src.protein_atlas.score_windows import get_windows, score_sequence_windowed

def test_get_windows():
    seq = "A" * 1500
    windows = get_windows(seq, max_len=900, overlap=100)

    assert len(windows) == 2
    assert windows[0][:2] == (0, 900)
    assert windows[1][:2] == (800, 1500)

    seq_short = "A" * 500
    windows_short = get_windows(seq_short, max_len=900, overlap=100)
    assert len(windows_short) == 1
    assert windows_short[0][:2] == (0, 500)

def test_score_sequence_windowed_merge():
    # A fake scoring function
    def fake_score_fn(scorer, sequence, **kwargs):
        L = len(sequence)
        # Mock surprisal bits uniquely identifying distance from start
        return {
            "residue_surprisals": list(range(L)),
            "residue_probabilities": [0.5] * L,
            "scoring_method": "mock"
        }

    seq = "A" * 1500

    res = score_sequence_windowed(
        scorer=None,
        sequence=seq,
        score_fn=fake_score_fn,
        max_len=900,
        overlap=100,
        merge_rule="central"
    )

    assert res["windowed"] is True
    assert res["number_of_windows"] == 2

    surp = res["residue_surprisals"]
    assert len(surp) == 1500
    assert None not in surp

    # In window 0 (0 to 900), pos 850 is distance 49 from right edge.
    # In window 1 (800 to 1500), global pos 850 is local pos 50, distance 50 from left edge.
    # So for global pos 850, window 1 is more central (dist 50 > dist 49).
    # In window 1, local pos 50 has mock surprisal 50.
    assert surp[850] == 50
    assert res["residue_window_ids"][850] == 1
