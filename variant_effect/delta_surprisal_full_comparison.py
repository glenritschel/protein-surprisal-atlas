# SELF-CONTAINED: ESM-2 650M delta-surprisal vs AlphaMissense vs BLOSUM62, 8-gene ClinVar set.
# Disconnect-proof: depends on no prior cell. Just re-run this one cell. GPU required.
import torch, os, re, urllib.request, functools, requests, numpy as np, pandas as pd
from transformers import AutoTokenizer, AutoModelForMaskedLM
from sklearn.metrics import roc_auc_score
from scipy.stats import rankdata
from Bio.Align import substitution_matrices
B = substitution_matrices.load("BLOSUM62")

MODEL = "facebook/esm2_t33_650M_UR50D"
tok = AutoTokenizer.from_pretrained(MODEL)
model = AutoModelForMaskedLM.from_pretrained(MODEL).eval()
DEV = "cuda" if torch.cuda.is_available() else "cpu"; model.to(DEV)
print("loaded", MODEL, "on", DEV)
assert DEV == "cuda", "No GPU — Runtime > Change runtime type > T4 GPU, then re-run."

W = 1022


@torch.no_grad()
def dsurp(seq, p, wt, mut):
    assert seq[p-1] == wt
    half = W // 2; start = max(0, (p-1) - half); end = min(len(seq), start + W); start = max(0, end - W)
    sub = seq[start:end]; loc = (p-1) - start + 1; assert sub[loc-1] == wt
    enc = tok(sub, return_tensors="pt").to(DEV); ids = enc["input_ids"]
    m = ids.clone(); m[0, loc] = tok.mask_token_id
    lp = torch.log_softmax(model(input_ids=m, attention_mask=enc["attention_mask"]).logits[0, loc], -1)
    return float((lp[tok.convert_tokens_to_ids(wt)] - lp[tok.convert_tokens_to_ids(mut)]).item())


@functools.lru_cache(maxsize=512)
def useq(acc):
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30); r.raise_for_status()
    return "".join(r.text.splitlines()[1:])


def clinvar(UP, cap=200):
    url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"; p = "/content/clinvar.txt.gz"
    if not os.path.exists(p):
        print("downloading ClinVar ~250MB ..."); urllib.request.urlretrieve(url, p)
    cols = ["Type", "Name", "GeneSymbol", "ClinicalSignificance"]; genes = set(UP); hits = []
    for ch in pd.read_csv(p, sep="\t", usecols=cols, dtype=str, chunksize=200000, compression="gzip"):
        s = ch[(ch.Type == "single nucleotide variant") & (ch.GeneSymbol.isin(genes))]
        if len(s): hits.append(s)
    df = pd.concat(hits, ignore_index=True)
    t2o = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G','His':'H',
           'Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
    pat = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})"); rows = []
    for _, r in df.iterrows():
        m = pat.search(str(r.Name)); sig = str(r.ClinicalSignificance).lower()
        if not m: continue
        w, pos, mt = m.group(1), int(m.group(2)), m.group(3)
        if w not in t2o or mt not in t2o: continue
        if "pathogenic" in sig and "conflict" not in sig: lab = 1
        elif "benign" in sig and "conflict" not in sig:   lab = 0
        else: continue
        rows.append(dict(gene=r.GeneSymbol, uniprot=UP[r.GeneSymbol], pos=pos, wt=t2o[w], mut=t2o[mt], label=lab))
    out = pd.DataFrame(rows).drop_duplicates(["uniprot", "pos", "wt", "mut"])
    return out.groupby(["gene", "label"], group_keys=False).apply(lambda g: g.sample(min(len(g), cap), random_state=0))


UP = {"TP53": "P04637", "PTEN": "P60484", "BRCA1": "P38398", "BRCA2": "P51587",
      "MLH1": "P40692", "MSH2": "P43246", "VHL": "P40337", "RET": "P07949"}
V = clinvar(UP); SEQ = {a: useq(a) for a in V.uniprot.unique()}
rows = []
for _, v in V.iterrows():
    try:
        d = dsurp(SEQ[v.uniprot], int(v.pos), v.wt, v.mut)
    except Exception:
        continue
    rows.append(dict(gene=v.gene, label=v.label, uniprot=v.uniprot, pos=int(v.pos),
                     wt=v.wt, mut=v.mut, ds=d, blosum=-float(B[v.wt, v.mut])))
M = pd.DataFrame(rows); M["key"] = M.wt + M.pos.astype(str) + M.mut
print(f"scored {len(M)} variants")

url = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"; ap = "/content/alphamissense_aa.tsv.gz"
if not os.path.exists(ap):
    print("downloading AlphaMissense ~1GB ..."); urllib.request.urlretrieve(url, ap)
want = set(M.uniprot.unique()); hits = []
for ch in pd.read_csv(ap, sep="\t", comment="#", chunksize=1_000_000, compression="gzip",
                      usecols=["uniprot_id", "protein_variant", "am_pathogenicity"],
                      dtype={"uniprot_id": str, "protein_variant": str, "am_pathogenicity": float}):
    s = ch[ch.uniprot_id.isin(want)]
    if len(s): hits.append(s)
AM = pd.concat(hits, ignore_index=True); amap = dict(zip(AM.uniprot_id + "_" + AM.protein_variant, AM.am_pathogenicity))
M["am"] = [amap.get(f"{u}_{k}", np.nan) for u, k in zip(M.uniprot, M.key)]
print(f"AlphaMissense coverage: {M.am.notna().mean():.1%}")


def ci(y, s, n=1000):
    rng = np.random.default_rng(0); y = np.asarray(y); s = np.asarray(s); idx = np.arange(len(y)); out = []
    for _ in range(n):
        b = rng.choice(idx, len(idx), replace=True)
        if len(set(y[b])) == 2:
            out.append(roc_auc_score(y[b], s[b]))
    return np.percentile(out, [2.5, 97.5])


def report(sub, name, cols):
    if sub.label.nunique() < 2:
        print(f"\n{name}: one class only (n={len(sub)}) — skipped"); return
    y = sub.label.values
    print(f"\n{name}  (n={len(sub)}, {int((y==1).sum())} path / {int((y==0).sum())} benign)")
    for col, lab in cols:
        d = sub.dropna(subset=[col])
        a = roc_auc_score(d.label.values, d[col].values); lo, hi = ci(d.label.values, d[col].values)
        print(f"    {lab:<22} AUROC={a:.3f}  95%CI[{lo:.3f},{hi:.3f}]")


cols = [("ds", "delta-surprisal-650M"), ("am", "AlphaMissense"), ("blosum", "-BLOSUM62")]
report(M, "POOLED", cols)
for g, sub in sorted(M.groupby("gene"), key=lambda kv: -len(kv[1])):
    report(sub, g, cols)

Mm = M.dropna(subset=["am"]); y = Mm.label.values
ens = (rankdata(Mm.ds) + rankdata(Mm.am)) / 2
print(f"\nENSEMBLE (rank-avg delta+AM)  AUROC={roc_auc_score(y, ens):.3f}   "
      f"[delta {roc_auc_score(y, Mm.ds.values):.3f}, AM {roc_auc_score(y, Mm.am.values):.3f}]")
