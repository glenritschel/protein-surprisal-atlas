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
        u_url = f"https://rest.uniprot.org/uniref/search?query=(uniprot_id:{uniprot_id})AND(identity:{identity})&fields=id,count,name"

        retries = 3
        backoff = 1
        for attempt in range(retries):
            try:
                r = requests.get(u_url, timeout=10)
                if r.status_code == 200:
                    d = r.json()
                    results = d.get('results', [])
                    if results:
                        import logging
                        if not hasattr(fetch_size, "has_logged"):
                            logging.info(f"Raw JSON response for {uniprot_id} identity {identity}: {d}")
                            fetch_size.has_logged = True
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
    import json
    import threading
    if os.environ.get("FAST_MOCK_UNIREF") == "1":
        import numpy as np
        df['uniref50_cluster_id'] = "UR50_" + df['uniprot_id']
        df['uniref90_cluster_size'] = np.random.randint(1, 100, size=len(df))
        df['uniref50_cluster_size'] = np.random.randint(1, 1000, size=len(df))
        df['log_family_size'] = np.log(df['uniref50_cluster_size'] + 1)
        return df

    from concurrent.futures import ThreadPoolExecutor, as_completed
    import tqdm

    cache_file = "data/external/uniref_cache.json"
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    cache = {}
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r') as f:
                cache = json.load(f)
        except json.JSONDecodeError:
            cache = {}

    cache_lock = threading.Lock()

    def save_cache_safely():
        with cache_lock:
            with open(cache_file, 'w') as f:
                json.dump(cache, f)

    def fetch_and_cache(uid):
        if uid in cache:
            return uid, cache[uid]
        res = get_uniref_cluster_sizes(uid)
        with cache_lock:
            cache[uid] = res
        return uid, res

    to_fetch = df['uniprot_id'].tolist()
    results_map = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(fetch_and_cache, uid): uid for uid in to_fetch}
        for i, future in enumerate(tqdm.tqdm(as_completed(futures), total=len(to_fetch))):
            uid, res = future.result()
            results_map[uid] = res
            # Save cache periodically to survive interruption
            if i > 0 and i % 100 == 0:
                save_cache_safely()

    # Final save
    save_cache_safely()

    results = [results_map[uid] for uid in to_fetch]

    df_res = pd.DataFrame(results)
    for col in df_res.columns:
        df[col] = df_res[col].values

    missing_count = df['log_family_size'].isna().sum()
    if missing_count > 0:
        import logging
        logging.warning(f"Failed to retrieve family size for {missing_count} proteins.")

    return df
