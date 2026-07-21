from typing import Dict, Any, Callable
from src.protein_atlas.esm_model import ESM2Scorer

def get_windows(sequence: str, max_len: int = 900, overlap: int = 100):
    L = len(sequence)
    if L <= max_len:
        return [(0, L, sequence)]

    windows = []
    step = max_len - overlap
    for i in range(0, L, step):
        start = i
        end = min(i + max_len, L)
        windows.append((start, end, sequence[start:end]))
        if end == L:
            break

    return windows

def score_sequence_windowed(scorer: ESM2Scorer, sequence: str, score_fn: Callable, max_len: int = 900, overlap: int = 100, merge_rule: str = "central", **kwargs) -> Dict[str, Any]:
    L = len(sequence)
    windows = get_windows(sequence, max_len, overlap)

    if len(windows) == 1:
        res = score_fn(scorer, sequence, **kwargs)
        res["windowed"] = False
        res["number_of_windows"] = 1
        return res

    merged_surprisals = [None] * L
    merged_probs = [None] * L
    max_dist_to_edge = [-1] * L

    sum_surprisals = [0.0] * L
    count_surprisals = [0] * L

    method = None

    for (start, end, subseq) in windows:
        res = score_fn(scorer, subseq, **kwargs)
        method = res["scoring_method"]
        sub_surp = res["residue_surprisals"]
        sub_prob = res["residue_probabilities"]

        w_len = end - start

        for i in range(w_len):
            global_pos = start + i
            dist = min(i, w_len - 1 - i)

            if merge_rule == "central":
                if dist > max_dist_to_edge[global_pos]:
                    max_dist_to_edge[global_pos] = dist
                    merged_surprisals[global_pos] = sub_surp[i]
                    merged_probs[global_pos] = sub_prob[i]
            elif merge_rule == "average":
                sum_surprisals[global_pos] += sub_surp[i]
                count_surprisals[global_pos] += 1

                if merged_probs[global_pos] is None:
                    merged_probs[global_pos] = sub_prob[i]
                else:
                    merged_probs[global_pos] = (merged_probs[global_pos] * (count_surprisals[global_pos]-1) + sub_prob[i]) / count_surprisals[global_pos]

    if merge_rule == "average":
        for i in range(L):
            merged_surprisals[i] = sum_surprisals[i] / count_surprisals[i]

    total_surprisal = sum(merged_surprisals)

    return {
        "total_surprisal_bits": total_surprisal,
        "bits_per_residue": total_surprisal / L,
        "residue_surprisals": merged_surprisals,
        "residue_probabilities": merged_probs,
        "scoring_method": method,
        "windowed": True,
        "number_of_windows": len(windows)
    }
