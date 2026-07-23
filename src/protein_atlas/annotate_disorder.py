import requests
import json
import time
import math
import os
import pandas as pd
from typing import Dict, Any, Optional
import logging

MOBIDB_API_URL = "https://mobidb.org/api/download?acc={accession}&format=json"
CACHE_FILE = "data/external/mobidb_cache.json"

def fetch_disorder_fraction(uniprot_id: str) -> Optional[float]:
    url = MOBIDB_API_URL.format(accession=uniprot_id)
    retries = 3
    backoff = 1
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                d = r.json()
                if "prediction-disorder-th_50" in d:
                    return float(d["prediction-disorder-th_50"].get("content_fraction", math.nan))
                else:
                    return math.nan
            elif r.status_code == 429:
                time.sleep(backoff)
                backoff *= 2
                continue
            elif r.status_code == 404:
                return math.nan
            else:
                return math.nan
        except requests.RequestException:
            time.sleep(backoff)
            backoff *= 2
            continue
    return math.nan

def load_cache() -> Dict[str, Any]:
    if os.path.exists(CACHE_FILE):
        try:
            with open(CACHE_FILE, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError:
            return {}
    return {}

def save_cache(cache: Dict[str, Any]):
    os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
    with open(CACHE_FILE, 'w') as f:
        json.dump(cache, f)

def get_disorder_fractions(df: pd.DataFrame, max_workers: int = 4) -> pd.DataFrame:
    df = df.copy()
    cache = load_cache()

    results = []

    from concurrent.futures import ThreadPoolExecutor
    import tqdm

    to_fetch = []
    for uid in df['uniprot_id']:
        if uid not in cache:
            to_fetch.append(uid)

    if to_fetch:
        print(f"Fetching disorder fraction for {len(to_fetch)} proteins...")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            fetched_results = list(tqdm.tqdm(executor.map(fetch_disorder_fraction, to_fetch), total=len(to_fetch)))

        for uid, val in zip(to_fetch, fetched_results):
            cache[uid] = val
        save_cache(cache)

    df['disorder_fraction'] = df['uniprot_id'].map(cache)

    missing = df['disorder_fraction'].isna().sum()
    if missing > 0:
        logging.warning(f"Disorder fraction missing for {missing} proteins.")

    return df
