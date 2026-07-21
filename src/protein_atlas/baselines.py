import math
from typing import Dict, List

VALID_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

HUMAN_BACKGROUND = {
    'A': 0.070, 'C': 0.023, 'D': 0.047, 'E': 0.071, 'F': 0.036,
    'G': 0.066, 'H': 0.026, 'I': 0.043, 'K': 0.057, 'L': 0.099,
    'M': 0.021, 'N': 0.036, 'P': 0.063, 'Q': 0.047, 'R': 0.056,
    'S': 0.083, 'T': 0.053, 'V': 0.060, 'W': 0.012, 'Y': 0.026
}

sum_freq = sum(HUMAN_BACKGROUND.values())
for k in HUMAN_BACKGROUND:
    HUMAN_BACKGROUND[k] /= sum_freq

def uniform_baseline(length: int) -> float:
    return length * math.log2(20)

def order0_baseline(sequence: str, background: Dict[str, float] = HUMAN_BACKGROUND) -> float:
    surprisal = 0.0
    for aa in sequence:
        if aa in background:
            surprisal -= math.log2(background[aa])
        else:
            surprisal -= math.log2(1.0 / 20.0)
    return surprisal
