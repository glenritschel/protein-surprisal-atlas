import os
import hashlib
import pandas as pd
import requests

GNOMAD_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
GNOMAD_FILE_PATH = "data/external/gnomad.v4.1.constraint_metrics.tsv"

def fetch_gnomad_metrics():
    os.makedirs(os.path.dirname(GNOMAD_FILE_PATH), exist_ok=True)
    if not os.path.exists(GNOMAD_FILE_PATH):
        print(f"Downloading gnomAD metrics from {GNOMAD_URL}...")
        response = requests.get(GNOMAD_URL)
        response.raise_for_status()
        with open(GNOMAD_FILE_PATH, 'wb') as f:
            f.write(response.content)
        print("Download complete.")

    sha256 = hashlib.sha256()
    with open(GNOMAD_FILE_PATH, 'rb') as f:
        for block in iter(lambda: f.read(4096), b""):
            sha256.update(block)

    return GNOMAD_FILE_PATH, sha256.hexdigest(), "v4.1", GNOMAD_URL

def process_gnomad_metrics(filepath):
    df = pd.read_csv(filepath, sep='\t')

    # Based on the user's feedback, we need to inspect the actual column names
    # But since we've already done that previously and know they are lof.oe_ci.upper and lof.pLI in v4.1, we use them.
    cols_to_keep = ['gene', 'gene_id', 'transcript', 'mane_select', 'canonical', 'lof.oe_ci.upper', 'lof.pLI']

    # Rename for consistency
    df = df[cols_to_keep].copy()
    df.rename(columns={'lof.oe_ci.upper': 'gnomad_loeuf', 'lof.pLI': 'gnomad_pli'}, inplace=True)

    # Priority: MANE = 1, Canonical = 2, other = 3
    df['priority'] = 3
    df.loc[df['canonical'] == True, 'priority'] = 2
    df.loc[df['mane_select'] == True, 'priority'] = 1

    # Filter to only ENSG ids
    df_ens = df[df['gene_id'].astype(str).str.startswith('ENSG')].copy()

    # Gene-level table: one row per gene id / gene symbol
    df_ens_gene = df_ens.sort_values('priority').drop_duplicates(subset=['gene_id'], keep='first')
    df_sym_gene = df_ens.sort_values('priority').drop_duplicates(subset=['gene'], keep='first')

    # Transcript-level table (fallback)
    df_transcript = df_ens.drop_duplicates(subset=['transcript'], keep='first')

    return df_ens_gene, df_sym_gene, df_transcript

def extract_base_ensts(ensembl_str):
    if pd.isna(ensembl_str):
        return []
    # format: "ENST00000275493.7 [P00533-1];ENST00000342916.7 [P00533-4];"
    tokens = [t.strip() for t in str(ensembl_str).split(';') if t.strip()]
    base_ids = []
    for token in tokens:
        # remove bracketed isoform info e.g. [P00533-1]
        t = token.split('[')[0].strip()
        # remove version suffix e.g. .7
        t = t.split('.')[0].strip()
        if t:
            base_ids.append(t)
    return base_ids

def join_gnomad_to_pilot(pilot_df, gnomad_ens_gene, gnomad_sym_gene, gnomad_transcript):
    merged = pilot_df.copy()

    merged['gnomad_loeuf'] = None
    merged['gnomad_pli'] = None
    merged['gnomad_join_method'] = None

    # extract first gene symbol
    merged['first_gene_symbol'] = merged['gene_symbol'].astype(str).str.split(' ').str[0]

    for idx, row in merged.iterrows():
        matched = False

        # 1. Primary join: by gene_symbol to gene-level table
        sym = row['first_gene_symbol']
        if sym in gnomad_sym_gene['gene'].values:
            match_row = gnomad_sym_gene[gnomad_sym_gene['gene'] == sym].iloc[0]
            merged.at[idx, 'gnomad_loeuf'] = match_row['gnomad_loeuf']
            merged.at[idx, 'gnomad_pli'] = match_row['gnomad_pli']
            merged.at[idx, 'gnomad_join_method'] = 'gene_symbol'
            matched = True
            continue

        # 2. Fallback join: ENST transcript match
        ensts = extract_base_ensts(row.get('Ensembl', ''))
        for enst in ensts:
            if enst in gnomad_transcript['transcript'].values:
                match_row = gnomad_transcript[gnomad_transcript['transcript'] == enst].iloc[0]
                merged.at[idx, 'gnomad_loeuf'] = match_row['gnomad_loeuf']
                merged.at[idx, 'gnomad_pli'] = match_row['gnomad_pli']
                merged.at[idx, 'gnomad_join_method'] = 'transcript'
                matched = True
                break

    merged.drop(columns=['first_gene_symbol'], errors='ignore', inplace=True)
    return merged
