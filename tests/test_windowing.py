import pytest
from src.protein_atlas.score_windows import get_windows

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
