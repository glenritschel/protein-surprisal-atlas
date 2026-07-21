import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

class ESM2Scorer:
    def __init__(self, model_name: str = "facebook/esm2_t12_35M_UR50D", device: str = "auto", precision: str = "float32"):
        self.model_name = model_name

        if device == "auto":
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self.device = torch.device(device)

        print(f"Loading {model_name} on {self.device}...")
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModelForMaskedLM.from_pretrained(model_name)

        if precision == "float16" and self.device.type == "cuda":
            self.model = self.model.half()

        self.model.to(self.device)
        self.model.eval()

        self.mask_token_id = self.tokenizer.mask_token_id
        self.vocab = self.tokenizer.get_vocab()

    def get_token_id(self, aa: str) -> int:
        return self.vocab.get(aa, self.tokenizer.unk_token_id)
