import os
import hashlib
import pandas as pd
import requests

DEPMAP_RELEASE_URL_EFFECT = "https://ndownloader.figshare.com/files/51064667"
DEPMAP_RELEASE_URL_ESSENTIAL = "https://ndownloader.figshare.com/files/51064916"
DEPMAP_DOI = "10.6084/m9.figshare.27993248"
DEPMAP_LABEL = "DepMap Public 24Q4"

DEPMAP_EFFECT_PATH = "data/external/CRISPRGeneEffect.csv"
DEPMAP_ESSENTIAL_PATH = "data/external/CRISPRInferredCommonEssentials.csv"

def fetch_depmap_data():
    os.makedirs(os.path.dirname(DEPMAP_EFFECT_PATH), exist_ok=True)

    # Download CRISPRGeneEffect
    if not os.path.exists(DEPMAP_EFFECT_PATH):
        print(f"Downloading {DEPMAP_EFFECT_PATH} from {DEPMAP_RELEASE_URL_EFFECT}...")
        response = requests.get(DEPMAP_RELEASE_URL_EFFECT, stream=True)
        response.raise_for_status()
        with open(DEPMAP_EFFECT_PATH, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print("Download complete.")

    # Download CRISPRInferredCommonEssentials
    if not os.path.exists(DEPMAP_ESSENTIAL_PATH):
        print(f"Downloading {DEPMAP_ESSENTIAL_PATH} from {DEPMAP_RELEASE_URL_ESSENTIAL}...")
        response = requests.get(DEPMAP_RELEASE_URL_ESSENTIAL)
        response.raise_for_status()
        with open(DEPMAP_ESSENTIAL_PATH, 'wb') as f:
            f.write(response.content)
        print("Download complete.")

    sha256_effect = hashlib.sha256()
    with open(DEPMAP_EFFECT_PATH, 'rb') as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256_effect.update(block)

    sha256_essential = hashlib.sha256()
    with open(DEPMAP_ESSENTIAL_PATH, 'rb') as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256_essential.update(block)

    return (
        DEPMAP_EFFECT_PATH, sha256_effect.hexdigest(),
        DEPMAP_ESSENTIAL_PATH, sha256_essential.hexdigest(),
        DEPMAP_LABEL, DEPMAP_DOI
    )

def parse_depmap_label(label: str) -> str:
    """Parses a DepMap gene label like 'TP53 (7157)' to just 'TP53'."""
    return label.split(' ')[0]

def process_depmap_data(effect_path, essential_path):
    # Read common essentials
    # Usually a single column 'Essentials' or 'Gene' depending on format. We will see.
    # Actually, from DepMap format, it's typically a list of genes. Let's read it.
    essentials_df = pd.read_csv(essential_path)
    if 'Essentials' in essentials_df.columns:
        essential_col = 'Essentials'
    else:
        essential_col = essentials_df.columns[0] # assume first column if 'Essentials' not there

    # Process genes strings to just symbol
    common_essentials = set(parse_depmap_label(g) for g in essentials_df[essential_col].dropna())

    # Calculate per-gene mean CRISPR effect memory-efficiently
    # CRISPRGeneEffect.csv is CellLines x Genes, first col usually 'ModelID'
    chunk_size = 500
    sum_effects = None
    count_effects = None

    for chunk in pd.read_csv(effect_path, chunksize=chunk_size):
        if 'ModelID' in chunk.columns:
            chunk = chunk.drop(columns=['ModelID'])
        elif 'Unnamed: 0' in chunk.columns:
            chunk = chunk.drop(columns=['Unnamed: 0'])

        if sum_effects is None:
            sum_effects = chunk.sum()
            count_effects = chunk.notna().sum()
        else:
            sum_effects += chunk.sum()
            count_effects += chunk.notna().sum()

    mean_effects = sum_effects / count_effects

    results = []
    for raw_label, mean_effect in mean_effects.items():
        sym = parse_depmap_label(raw_label)
        results.append({
            'gene_symbol': sym,
            'depmap_gene_effect': mean_effect,
            'depmap_common_essential': sym in common_essentials
        })

    df = pd.DataFrame(results)

    # In case of duplicate gene symbols after stripping (rare, but possible), take mean
    df = df.groupby('gene_symbol').agg({
        'depmap_gene_effect': 'mean',
        'depmap_common_essential': 'any'
    }).reset_index()

    return df

def join_depmap_to_pilot(pilot_df, depmap_df):
    merged = pilot_df.copy()

    merged['depmap_gene_effect'] = None
    merged['depmap_common_essential'] = None

    # Pilot gene_symbol can contain space-separated synonyms like "SMIM26 LINC00493"
    merged['first_gene_symbol'] = merged['gene_symbol'].astype(str).str.split(' ').str[0]

    # create mapping
    effect_map = dict(zip(depmap_df['gene_symbol'], depmap_df['depmap_gene_effect']))
    essential_map = dict(zip(depmap_df['gene_symbol'], depmap_df['depmap_common_essential']))

    merged['depmap_gene_effect'] = merged['first_gene_symbol'].map(effect_map)
    merged['depmap_common_essential'] = merged['first_gene_symbol'].map(essential_map)

    merged.drop(columns=['first_gene_symbol'], errors='ignore', inplace=True)
    return merged
