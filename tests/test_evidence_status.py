from unittest.mock import Mock

import requests

from services.target_enrichment import TargetEnrichmentService


def test_enrichment_marks_successful_empty_services_as_empty():
    service = TargetEnrichmentService()
    service.get_uniprot = Mock(return_value=[])
    service.get_pubmed_ids = Mock(return_value=[])
    service.get_reactome_pathways = Mock(return_value=[])
    service.get_string_network = Mock(
        return_value={"nodes": [], "edges": []}
    )
    service.get_string_functional_enrichment = Mock(return_value=[])

    target = Mock(symbol="LDLR")
    result = service.enrich_target(target, "familial hypercholesterolemia")

    statuses = result["evidence_status"]

    assert statuses["uniprot"]["status"] == "empty"
    assert statuses["publications"]["status"] == "empty"
    assert statuses["reactome_pathways"]["status"] == "empty"
    assert statuses["string_network"]["status"] == "empty"
    assert statuses["go_enrichment"]["status"] == "empty"
    assert statuses["kegg_enrichment"]["status"] == "empty"
    assert result["enrichment_errors"] == []


def test_enrichment_marks_failed_services_without_losing_other_results():
    service = TargetEnrichmentService()
    service.get_uniprot = Mock(
        side_effect=requests.ConnectionError("UniProt unavailable")
    )
    service.get_pubmed_ids = Mock(return_value=[{
        "pmid": "12345",
        "title": "Example publication",
        "pubdate": "2026",
        "url": "https://pubmed.ncbi.nlm.nih.gov/12345/",
    }])
    service.get_reactome_pathways = Mock(return_value=[])
    service.get_string_network = Mock(
        return_value={"nodes": [], "edges": []}
    )
    service.get_string_functional_enrichment = Mock(return_value=[])

    target = Mock(symbol="LDLR")
    result = service.enrich_target(target, "familial hypercholesterolemia")

    assert result["evidence_status"]["uniprot"]["status"] == "failed"
    assert result["evidence_status"]["publications"]["status"] == "retrieved"
    assert len(result["publications"]) == 1
    assert any(
        issue["service"] == "UniProt"
        for issue in result["enrichment_errors"]
    )


def test_string_network_status_includes_actual_node_and_edge_counts():
    service = TargetEnrichmentService()
    service.get_uniprot = Mock(return_value=[])
    service.get_pubmed_ids = Mock(return_value=[])
    service.get_reactome_pathways = Mock(return_value=[])
    service.get_string_network = Mock(return_value={
        "nodes": [
            {"preferred_name": "LDLR"},
            {"preferred_name": "APOB"},
        ],
        "edges": [{"protein_a": "LDLR", "protein_b": "APOB"}],
    })
    service.get_string_functional_enrichment = Mock(return_value=[])

    target = Mock(symbol="LDLR")
    result = service.enrich_target(target, "familial hypercholesterolemia")

    status = result["evidence_status"]["string_network"]
    assert status["status"] == "retrieved"
    assert status["node_count"] == 2
    assert status["record_count"] == 1
