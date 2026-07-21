import torch
import math
from typing import List, Dict, Any
import pandas as pd

def score_variants(scorer, sequence: str, variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Scores missense variants using the masked-marginal log-likelihood ratio.
    variant: dict with 'position' (1-based), 'reference', 'alternate'
    Delta S = -log2 P(alt | context) + log2 P(ref | context)
    Positive Delta S: alt is less probable.
    """
    L = len(sequence)
    token_ids = scorer.tokenizer(sequence, add_special_tokens=True, return_tensors="pt")["input_ids"][0]

    results = []

    with torch.no_grad():
        for var in variants:
            pos_1based = var['position']
            ref_aa = var['reference']
            alt_aa = var['alternate']

            if pos_1based < 1 or pos_1based > L:
                continue

            idx = pos_1based - 1
            if sequence[idx] != ref_aa:
                continue # Ref mismatch

            variant = token_ids.clone()
            variant[idx+1] = scorer.mask_token_id

            inputs = variant.unsqueeze(0).to(scorer.device)
            outputs = scorer.model(inputs)
            logits = outputs.logits[0, idx+1]

            log_probs = torch.nn.functional.log_softmax(logits, dim=-1)

            ref_token = scorer.get_token_id(ref_aa)
            alt_token = scorer.get_token_id(alt_aa)

            ref_log_prob = log_probs[ref_token].item()
            alt_log_prob = log_probs[alt_token].item()

            ref_surprisal = - (ref_log_prob / math.log(2))
            alt_surprisal = - (alt_log_prob / math.log(2))

            delta_s = alt_surprisal - ref_surprisal

            results.append({
                "position": pos_1based,
                "reference_residue": ref_aa,
                "alternate_residue": alt_aa,
                "ref_surprisal_bits": ref_surprisal,
                "alt_surprisal_bits": alt_surprisal,
                "delta_s": delta_s
            })

    return results

def score_variants_batch(scorer, sequence: str, variants: List[Dict[str, Any]], batch_size: int = 32) -> List[Dict[str, Any]]:
    """
    More efficient batch version of score_variants.
    """
    L = len(sequence)
    token_ids = scorer.tokenizer(sequence, add_special_tokens=True, return_tensors="pt")["input_ids"][0]

    valid_variants = []
    masked_seqs = []

    for var in variants:
        pos_1based = var['position']
        ref_aa = var['reference']
        if pos_1based < 1 or pos_1based > L:
            continue
        idx = pos_1based - 1
        if sequence[idx] != ref_aa:
            continue

        variant_toks = token_ids.clone()
        variant_toks[idx+1] = scorer.mask_token_id
        masked_seqs.append(variant_toks)
        valid_variants.append(var)

    if not masked_seqs:
        return []

    masked_seqs = torch.stack(masked_seqs)

    results = []

    with torch.no_grad():
        for i in range(0, len(masked_seqs), batch_size):
            batch = masked_seqs[i:i+batch_size].to(scorer.device)
            outputs = scorer.model(batch)

            for j in range(batch.size(0)):
                idx = i + j
                var = valid_variants[idx]
                pos_1based = var['position']
                seq_idx = pos_1based - 1

                logits = outputs.logits[j, seq_idx+1]
                log_probs = torch.nn.functional.log_softmax(logits, dim=-1)

                ref_token = scorer.get_token_id(var['reference'])
                alt_token = scorer.get_token_id(var['alternate'])

                ref_log_prob = log_probs[ref_token].item()
                alt_log_prob = log_probs[alt_token].item()

                ref_surprisal = - (ref_log_prob / math.log(2))
                alt_surprisal = - (alt_log_prob / math.log(2))

                delta_s = alt_surprisal - ref_surprisal

                results.append({
                    "position": pos_1based,
                    "reference_residue": var['reference'],
                    "alternate_residue": var['alternate'],
                    "ref_surprisal_bits": ref_surprisal,
                    "alt_surprisal_bits": alt_surprisal,
                    "delta_s": delta_s
                })

    return results
