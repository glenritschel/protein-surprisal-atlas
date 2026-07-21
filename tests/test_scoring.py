import pytest
import math
import numpy as np

def test_surprisal_arithmetic():
    probs = [0.5, 0.25, 0.125]
    surprisals = [-math.log2(p) for p in probs]
    assert np.allclose(surprisals, [1.0, 2.0, 3.0])
    assert sum(surprisals) == 6.0
