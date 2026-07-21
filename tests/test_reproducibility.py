import pytest
import pandas as pd
from src.protein_atlas.select_pilot import select_pilot_set
from src.protein_atlas.score_approximate import score_sequence_approximate

def test_deterministic_pilot():
    df = pd.DataFrame({
        "uniprot_id": [f"P{i}" for i in range(100)],
        "gene_symbol": [f"G{i}" for i in range(100)],
        "sequence_length": list(range(100, 200))
    })

    pilot1 = select_pilot_set(df, target_size=10, seed=42)
    pilot2 = select_pilot_set(df, target_size=10, seed=42)

    assert pilot1.equals(pilot2)

    pilot3 = select_pilot_set(df, target_size=10, seed=43)
    assert not pilot1.equals(pilot3)

class MockModel:
    def __call__(self, inputs):
        class Output:
            def __init__(self, logits):
                self.logits = logits
        import torch
        # Model gets a batch of sequences. `outputs.logits[0]` means we expect the logits
        # to have shape (batch_size, seq_len, vocab_size).
        # Here inputs.shape is (1, L)
        B, L = inputs.shape
        # Create logits for the single sequence
        logits = torch.arange(L * 150).float().view(B, L, 150) % 10.0
        # The scorer takes outputs.logits[0], so we need outputs.logits to have shape (B, L, 150)
        return Output(logits) # Notice I passed logits directly here, not [logits], as output.logits should be the tensor

class MockScorer:
    def __init__(self):
        self.device = "cpu"
        self.mask_token_id = 99
        self.model = MockModel()

    def tokenizer(self, seq, add_special_tokens=True, return_tensors="pt"):
        import torch
        # A tiny tokenizer stub: sequence -> ascii values, with pseudo special tokens at ends
        ids = [1] + [ord(c) for c in seq] + [2]
        return {"input_ids": [torch.tensor(ids)]}

def test_score_sequence_approximate_reproducibility():
    scorer = MockScorer()
    seq = "ACDEFGHIKLMNPQRSTVWY"

    res1 = score_sequence_approximate(scorer, seq, passes=3)
    res2 = score_sequence_approximate(scorer, seq, passes=3)

    assert res1["total_surprisal_bits"] == res2["total_surprisal_bits"]
    assert res1["residue_surprisals"] == res2["residue_surprisals"]
    assert res1["residue_probabilities"] == res2["residue_probabilities"]
