import math
import numpy as np
from src.protein_atlas.annotate_lowcomplexity import compute_entropy, compute_low_complexity_fraction

def test_entropy_computation():
    # Homopolymer has 0 entropy
    assert compute_entropy("AAAAAAAAAAAAAAAAAAAA") == 0.0

    # Max entropy for 20 distinct amino acids is log2(20) ~ 4.32
    seq_20_distinct = "ACDEFGHIKLMNPQRSTVWY"
    expected = math.log2(20)
    assert math.isclose(compute_entropy(seq_20_distinct), expected)

def test_low_complexity_fraction():
    # Sequence with a homopolymer run in the middle
    # 10 diverse + 20 A's + 10 diverse
    seq = "ACDEFGHIKL" + "A" * 20 + "MNPQRSTVWY"

    fraction = compute_low_complexity_fraction(seq, window_size=20, entropy_threshold=0.1)

    # The only window with H < 0.1 is the pure 20 A's. So exactly 20 residues are flagged.
    # Total length = 40.
    assert math.isclose(fraction, 20 / 40)

    fraction_0 = compute_low_complexity_fraction(seq, window_size=20, entropy_threshold=0.01)
    assert math.isclose(fraction_0, 0.5)

def test_short_sequence():
    assert compute_low_complexity_fraction("ACD", window_size=20) == 0.0
