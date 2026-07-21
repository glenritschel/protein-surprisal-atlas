import pytest
import math
import numpy as np
from src.protein_atlas.baselines import uniform_baseline, order0_baseline, HUMAN_BACKGROUND

def test_uniform_baseline():
    L = 10
    expected = L * math.log2(20)
    assert np.isclose(uniform_baseline(L), expected)

def test_order0_baseline():
    seq = "A"
    expected = -math.log2(HUMAN_BACKGROUND['A'])
    assert np.isclose(order0_baseline(seq), expected)

    seq = "AC"
    expected = -math.log2(HUMAN_BACKGROUND['A']) - math.log2(HUMAN_BACKGROUND['C'])
    assert np.isclose(order0_baseline(seq), expected)

def test_surprisal_ratio():
    seq = "A"
    L = len(seq)
    S = 2.0
    uni_b = uniform_baseline(L)
    sr_uni = S / uni_b
    assert np.isclose(sr_uni, 2.0 / math.log2(20))
