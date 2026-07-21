import requests
import pandas as pd
import io
import os
import hashlib
import datetime
from pathlib import Path

def download_human_proteome(output_dir: str):
    """
    Downloads the reviewed canonical human proteome from UniProtKB/Swiss-Prot.
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    url = "https://rest.uniprot.org/uniprotkb/stream"

    fields = (
        "accession,gene_names,protein_name,sequence,length,reviewed,"
        "organism_id,cc_subcellular_location,protein_families,ft_domain,"
        "ft_transmem,ft_signal,cc_disease,cc_catalytic_activity,"
        "cc_cofactor,go_id,xref_ensembl,xref_hgnc"
    )

    params = {
        "query": "organism_id:9606 AND reviewed:true",
        "format": "tsv",
        "fields": fields
    }

    print(f"Downloading from UniProt: {url}")
    response = requests.get(url, params=params, stream=True)
    response.raise_for_status()

    raw_file_path = Path(output_dir) / "uniprot_human_reviewed.tsv"
    with open(raw_file_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)

    checksum = hashlib.sha256()
    with open(raw_file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            checksum.update(chunk)

    checksum_hex = checksum.hexdigest()

    df = pd.read_csv(raw_file_path, sep="\t")

    metadata = {
        "download_date": datetime.datetime.now().isoformat(),
        "source_url": url,
        "query": params["query"],
        "checksum_sha256": checksum_hex
    }

    return df, metadata
