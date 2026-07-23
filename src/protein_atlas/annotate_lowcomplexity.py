import math
import numpy as np

def compute_entropy(window_seq: str) -> float:
    counts = {}
    for char in window_seq:
        counts[char] = counts.get(char, 0) + 1

    length = len(window_seq)
    entropy = 0.0
    for count in counts.values():
        p = count / length
        if p > 0:
            entropy -= p * math.log2(p)
    return entropy

def compute_low_complexity_fraction(sequence: str, window_size: int = 20, entropy_threshold: float = 3.0) -> float:
    """
    Computes the fraction of residues that are covered by at least one
    window of `window_size` which has an amino acid Shannon entropy
    below `entropy_threshold`.
    """
    if len(sequence) < window_size:
        return 0.0

    is_low_complexity = np.zeros(len(sequence), dtype=bool)

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        entropy = compute_entropy(window)
        if entropy < entropy_threshold:
            is_low_complexity[i:i + window_size] = True

    return float(np.mean(is_low_complexity))
