import pytest
import requests
from unittest.mock import Mock

from services.target_enrichment import TargetEnrichmentService


def make_response(data=None, status_code=200):
    """Create a mocked HTTP response."""
    response = Mock()
    response.status_code = status_code
    response.raise_for_status.return_value = None
    response.json.return_value = data
    return response


def test_uniprot_http_error():
    service = TargetEnrichmentService()
    response = make_response()
    response.raise_for_status.side_effect = requests.HTTPError("HTTP 500")
    service.session.get = Mock(return_value=response)

    with pytest.raises(requests.HTTPError):
        service.get_uniprot("LDLR")


def test_uniprot_malformed_json():
    service = TargetEnrichmentService()
    response = make_response()
    response.json.side_effect = ValueError("Invalid JSON")
    service.session.get = Mock(return_value=response)

    with pytest.raises(ValueError):
        service.get_uniprot("LDLR")


def test_pubmed_connection_failure():
    service = TargetEnrichmentService()
    service.session.get = Mock(
        side_effect=requests.ConnectionError("Connection failed")
    )

    with pytest.raises(requests.ConnectionError):
        service.get_pubmed_ids("LDLR", "familial hypercholesterolemia")


def test_pubmed_empty_results():
    service = TargetEnrichmentService()
    response = make_response({"esearchresult": {"idlist": []}})
    service.session.get = Mock(return_value=response)

    result = service.get_pubmed_ids("LDLR", "familial hypercholesterolemia")

    assert result == []
    assert service.session.get.call_count == 1


def test_reactome_not_found():
    service = TargetEnrichmentService()
    response = make_response(status_code=404)
    service.session.get = Mock(return_value=response)

    result = service.get_reactome_pathways("INVALID_ACCESSION")

    assert result == []


def test_open_targets_graphql_error():
    service = TargetEnrichmentService()
    response = make_response({
        "errors": [{"message": "Invalid GraphQL query"}]
    })
    service.session.post = Mock(return_value=response)

    with pytest.raises(RuntimeError, match="GraphQL error"):
        service.get_known_drugs("ENSG00000130164")


def test_open_targets_missing_target():
    service = TargetEnrichmentService()
    response = make_response({"data": {"target": None}})
    service.session.post = Mock(return_value=response)

    result = service.get_known_drugs("ENSG00000130164")

    assert result == []


def test_uniprot_empty_results():
    service = TargetEnrichmentService()
    response = make_response({"results": []})
    service.session.get = Mock(return_value=response)

    result = service.get_uniprot("NONEXISTENT")

    assert result == []
