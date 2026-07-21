import torch
import math
import random
from typing import Dict, Any
from src.protein_atlas.esm_model import ESM2Scorer

def score_sequence_approximate(scorer: ESM2Scorer, sequence: str, passes: int = 7) -> Dict[str, Any]:
    L = len(sequence)
    token_ids = scorer.tokenizer(sequence, add_special_tokens=True, return_tensors="pt")["input_ids"][0]

    surprisals = [None] * L
    probabilities = [None] * L

    groups = []
    for i in range(passes):
        group = [pos for pos in range(L) if pos % passes == i]
        if group:
            groups.append(group)

    with torch.no_grad():
        for group in groups:
            variant = token_ids.clone()

            for pos in group:
                variant[pos+1] = scorer.mask_token_id

            variant_t = variant.unsqueeze(0).to(scorer.device)
            outputs = scorer.model(variant_t)
            logits = outputs.logits[0]

            log_probs = torch.nn.functional.log_softmax(logits, dim=-1)
            probs = torch.exp(log_probs)

            for pos in group:
                target_token = token_ids[pos+1].item()
                pos_log_prob = log_probs[pos+1, target_token].item()
                pos_prob = probs[pos+1, target_token].item()

                surprisal = - (pos_log_prob / math.log(2))

                surprisals[pos] = surprisal
                probabilities[pos] = pos_prob

    total_surprisal = sum(surprisals)

    return {
        "total_surprisal_bits": total_surprisal,
        "bits_per_residue": total_surprisal / L,
        "residue_surprisals": surprisals,
        "residue_probabilities": probabilities,
        "scoring_method": "sampled_mask"
    }

def score_sequence_naive(scorer: ESM2Scorer, sequence: str) -> Dict[str, Any]:
    L = len(sequence)
    token_ids = scorer.tokenizer(sequence, add_special_tokens=True, return_tensors="pt")["input_ids"][0]

    surprisals = []
    probabilities = []

    with torch.no_grad():
        inputs = token_ids.unsqueeze(0).to(scorer.device)
        outputs = scorer.model(inputs)
        logits = outputs.logits[0]

        log_probs = torch.nn.functional.log_softmax(logits, dim=-1)
        probs = torch.exp(log_probs)

        for i in range(L):
            target_token = token_ids[i+1].item()
            pos_log_prob = log_probs[i+1, target_token].item()
            pos_prob = probs[i+1, target_token].item()

            surprisal = - (pos_log_prob / math.log(2))
            surprisals.append(surprisal)
            probabilities.append(pos_prob)

    total_surprisal = sum(surprisals)

    return {
        "total_surprisal_bits": total_surprisal,
        "bits_per_residue": total_surprisal / L,
        "residue_surprisals": surprisals,
        "residue_probabilities": probabilities,
        "scoring_method": "naive_onepass"
    }
