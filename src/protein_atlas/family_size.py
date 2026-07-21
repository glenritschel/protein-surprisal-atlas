import requests
import math
from typing import Optional, Dict, Any, Tuple
import pandas as pd

import time

def get_uniref_cluster_sizes(uniprot_id: str) -> Dict[str, Any]:
    uniref90_size = None
    uniref50_size = None
    uniref50_id = None

    def fetch_size(identity: str) -> Tuple[Optional[str], Optional[int]]:
        u_url = f"https://rest.uniprot.org/uniref/search?query=(member:{uniprot_id})AND(identity:{identity})&fields=id,count,name"

        retries = 3
        backoff = 1
        for attempt in range(retries):
            try:
                r = requests.get(u_url, timeout=10)
                if r.status_code == 200:
                    d = r.json()
                    results = d.get('results', [])
                    if results:
                        # Take the first matched cluster
                        cluster = results[0]
                        return cluster.get('id'), cluster.get('memberCount')
                    else:
                        # Empty results means not found
                        return None, None
                elif r.status_code == 429:
                    # Throttled, wait and retry
                    time.sleep(backoff)
                    backoff *= 2
                    continue
                else:
                    # Other HTTP error
                    return None, None
            except requests.RequestException:
                time.sleep(backoff)
                backoff *= 2
                continue
        return None, None

    _, uniref90_size = fetch_size("0.9")
    uniref50_id, uniref50_size = fetch_size("0.5")

    log_family_size = math.nan
    if uniref50_size is not None:
        log_family_size = math.log(uniref50_size + 1)

    return {
        "uniref50_cluster_id": uniref50_id if uniref50_id else f"Unknown_{uniprot_id}",
        "uniref90_cluster_size": uniref90_size,
        "uniref50_cluster_size": uniref50_size,
        "log_family_size": log_family_size
    }

def add_family_sizes(df, max_workers=4):
    import os
    if os.environ.get("FAST_MOCK_UNIREF") == "1":
        import numpy as np
        df['uniref50_cluster_id'] = "UR50_" + df['uniprot_id']
        df['uniref90_cluster_size'] = np.random.randint(1, 100, size=len(df))
        df['uniref50_cluster_size'] = np.random.randint(1, 1000, size=len(df))
        df['log_family_size'] = np.log(df['uniref50_cluster_size'] + 1)
        return df

    from concurrent.futures import ThreadPoolExecutor
    import tqdm

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for res in tqdm.tqdm(executor.map(get_uniref_cluster_sizes, df['uniprot_id']), total=len(df)):
            results.append(res)

    df_res = pd.DataFrame(results)
    for col in df_res.columns:
        df[col] = df_res[col].values

    missing_count = df['log_family_size'].isna().sum()
    if missing_count > 0:
        import logging
        logging.warning(f"Failed to retrieve family size for {missing_count} proteins.")

    return df
