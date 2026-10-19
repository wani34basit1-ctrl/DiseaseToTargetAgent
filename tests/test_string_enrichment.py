import requests
from unittest.mock import Mock

import pytest

from services.target_enrichment import TargetEnrichmentService


def response(data, status_code=200):
    item = Mock()
    item.status_code = status_code
    item.raise_for_status.return_value = None
    item.json.return_value = data
    return item


def test_string_network_normalizes_nodes_edges_and_scores():
    service = TargetEnrichmentService()
    service.session.get = Mock(return_value=response([{
        "stringId_A": "9606.ENSP1",
        "stringId_B": "9606.ENSP2",
        "preferredName_A": "LDLR",
        "preferredName_B": "APOB",
        "score": 0.9,
        "nscore": 0.8,
    }]))

    result = service.get_string_network("LDLR")

    assert {node["preferred_name"] for node in result["nodes"]} == {"LDLR", "APOB"}
    assert result["edges"][0]["score_0_1000"] == 900
    assert result["edges"][0]["evidence_scores"]["nscore"] == 0.8
    assert service.session.get.call_args.args[0].endswith("/network")


def test_string_network_empty_gene_returns_empty_network_without_request():
    service = TargetEnrichmentService()
    service.session.get = Mock()

    assert service.get_string_network(" ") == {"nodes": [], "edges": []}
    service.session.get.assert_not_called()


def test_string_functional_enrichment_returns_go_and_kegg_terms():
    service = TargetEnrichmentService()
    service.session.get = Mock(return_value=response([
        {"category": "Process", "term": "GO:0001", "description": "Example process", "fdr": 0.01},
        {"category": "KEGG", "term": "hsa00010", "description": "Example pathway", "fdr": 0.02},
    ]))

    result = service.get_string_functional_enrichment(["LDLR", "APOB", "LDLR"])

    assert len(result) == 2
    assert result[0]["term_id"] == "GO:0001"
    assert result[1]["category"] == "KEGG"
    assert result[1]["fdr"] == 0.02
    assert service.session.get.call_args.kwargs["params"]["identifiers"] == "LDLR\rAPOB"


def test_string_functional_enrichment_propagates_http_errors():
    service = TargetEnrichmentService()
    item = response([])
    item.raise_for_status.side_effect = requests.HTTPError("HTTP 500")
    service.session.get = Mock(return_value=item)

    with pytest.raises(requests.HTTPError):
        service.get_string_functional_enrichment(["LDLR"])


def test_enrich_target_keeps_existing_results_if_string_fails():
    service = TargetEnrichmentService()
    service.get_uniprot = Mock(return_value=[{"accession": "P01130", "protein_name": "LDL receptor", "url": "https://example.org"}])
    service.get_pubmed_ids = Mock(return_value=[])
    service.get_reactome_pathways = Mock(return_value=[{"pathway_id": "R-HSA-1", "name": "Example", "url": "https://example.org"}])
    service.get_string_network = Mock(side_effect=requests.ConnectionError("offline"))
    service.get_string_functional_enrichment = Mock(return_value=[])

    target = Mock(symbol="LDLR")
    result = service.enrich_target(target, "familial hypercholesterolemia")

    assert result["uniprot_entries"][0]["accession"] == "P01130"
    assert result["pathways"][0]["pathway_id"] == "R-HSA-1"
    assert result["string_network"] == {"nodes": [], "edges": []}
    assert result["enrichment_errors"][0]["service"] == "STRING network"
