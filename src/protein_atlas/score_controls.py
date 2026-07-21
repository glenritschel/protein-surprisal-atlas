import random
import numpy as np
from src.protein_atlas.baselines import HUMAN_BACKGROUND

def shuffle_sequence(seq: str, seed: int = None) -> str:
    if seed is not None:
        random.seed(seed)
    seq_list = list(seq)
    random.shuffle(seq_list)
    return "".join(seq_list)

def generate_random_uniform(length: int, seed: int = None) -> str:
    if seed is not None:
        np.random.seed(seed)
    aas = list(HUMAN_BACKGROUND.keys())
    return "".join(np.random.choice(aas, size=length))

def generate_random_background(length: int, seed: int = None) -> str:
    if seed is not None:
        np.random.seed(seed)
    aas = list(HUMAN_BACKGROUND.keys())
    probs = list(HUMAN_BACKGROUND.values())
    return "".join(np.random.choice(aas, p=probs, size=length))

def reverse_sequence(seq: str) -> str:
    return seq[::-1]
