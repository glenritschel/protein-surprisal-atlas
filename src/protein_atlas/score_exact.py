import torch
import math
from typing import Dict, Any
from src.protein_atlas.esm_model import ESM2Scorer

def score_sequence_exact(scorer: ESM2Scorer, sequence: str, batch_size: int = 32) -> Dict[str, Any]:
    L = len(sequence)
    token_ids = scorer.tokenizer(sequence, add_special_tokens=True, return_tensors="pt")["input_ids"][0]

    masked_variants = []
    target_ids = []

    for i in range(L):
        variant = token_ids.clone()
        target_ids.append(variant[i+1].item())
        variant[i+1] = scorer.mask_token_id
        masked_variants.append(variant)

    masked_variants = torch.stack(masked_variants)

    surprisals = []
    probabilities = []

    with torch.no_grad():
        for i in range(0, L, batch_size):
            batch = masked_variants[i:i+batch_size].to(scorer.device)
            outputs = scorer.model(batch)
            logits = outputs.logits

            log_probs = torch.nn.functional.log_softmax(logits, dim=-1)
            probs = torch.exp(log_probs)

            for j in range(batch.size(0)):
                idx = i + j
                target_token = target_ids[idx]
                pos_log_prob = log_probs[j, idx+1, target_token].item()
                pos_prob = probs[j, idx+1, target_token].item()

                surprisal = - (pos_log_prob / math.log(2))

                surprisals.append(surprisal)
                probabilities.append(pos_prob)

    total_surprisal = sum(surprisals)

    return {
        "total_surprisal_bits": total_surprisal,
        "bits_per_residue": total_surprisal / L,
        "residue_surprisals": surprisals,
        "residue_probabilities": probabilities,
        "scoring_method": "exact"
    }
