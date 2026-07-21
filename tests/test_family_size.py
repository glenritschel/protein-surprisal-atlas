import pytest
import os
import pandas as pd
import responses
from src.protein_atlas.family_size import get_uniref_cluster_sizes, add_family_sizes

@responses.activate
def test_get_uniref_cluster_sizes_live_mocked():
    # Mock UniRef90
    responses.add(
        responses.GET,
        "https://rest.uniprot.org/uniref/search?query=(uniprot_id:P12345)AND(identity:0.9)&fields=id,count,name",
        json={"results": [{"id": "UniRef90_P00000", "memberCount": 5}]},
        status=200
    )
    # Mock UniRef50
    responses.add(
        responses.GET,
        "https://rest.uniprot.org/uniref/search?query=(uniprot_id:P12345)AND(identity:0.5)&fields=id,count,name",
        json={"results": [{"id": "UniRef50_P00000", "memberCount": 150}]},
        status=200
    )

    res = get_uniref_cluster_sizes("P12345")
    assert res["uniref90_cluster_size"] == 5
    assert res["uniref50_cluster_size"] == 150
    assert res["uniref50_cluster_id"] == "UniRef50_P00000"

@responses.activate
def test_get_uniref_cluster_sizes_not_found():
    # Return empty results
    responses.add(
        responses.GET,
        "https://rest.uniprot.org/uniref/search?query=(member:P99999)AND(identity:0.9)&fields=id,count,name",
        json={"results": []},
        status=200
    )
    responses.add(
        responses.GET,
        "https://rest.uniprot.org/uniref/search?query=(member:P99999)AND(identity:0.5)&fields=id,count,name",
        json={"results": []},
        status=200
    )

    res = get_uniref_cluster_sizes("P99999")
    import math
    assert res["uniref90_cluster_size"] is None
    assert res["uniref50_cluster_size"] is None
    assert math.isnan(res["log_family_size"])

def test_add_family_sizes_fast_mock(monkeypatch):
    monkeypatch.setenv("FAST_MOCK_UNIREF", "1")
    df = pd.DataFrame({"uniprot_id": ["P1", "P2"]})
    res_df = add_family_sizes(df)
    assert "uniref50_cluster_size" in res_df.columns
    assert "log_family_size" in res_df.columns
