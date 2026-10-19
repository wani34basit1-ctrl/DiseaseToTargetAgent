from unittest.mock import Mock

from services.disease_treatments import get_disease_treatments
from services.target_enrichment import TargetEnrichmentService


def mock_response(payload):
    response = Mock()
    response.raise_for_status.return_value = None
    response.json.return_value = payload
    return response


def test_disease_treatments_have_disease_level_provenance():
    treatments = get_disease_treatments("familial hypercholesterolemia")

    assert len(treatments) == 5
    assert all(
        item["evidence_level"] == "disease-level"
        for item in treatments
    )
    assert all(
        item["evidence_source"] == "curated_disease_treatment_catalogue"
        for item in treatments
    )


def test_disease_treatment_alias_preserves_provenance():
    treatments = get_disease_treatments("FH")

    assert treatments
    assert all(item["evidence_level"] == "disease-level" for item in treatments)


def test_open_targets_empty_candidates_return_empty_list():
    service = TargetEnrichmentService()
    service.session.post = Mock(return_value=mock_response({
        "data": {
            "target": {
                "id": "ENSG00000130164",
                "drugAndClinicalCandidates": {
                    "count": 0,
                    "rows": [],
                },
            }
        }
    }))

    assert service.get_known_drugs("ENSG00000130164") == []


def test_open_targets_candidates_have_target_associated_provenance():
    service = TargetEnrichmentService()
    service.session.post = Mock(return_value=mock_response({
        "data": {
            "target": {
                "id": "ENSG00000130164",
                "drugAndClinicalCandidates": {
                    "count": 1,
                    "rows": [{
                        "maxClinicalStage": 3,
                        "drug": {
                            "id": "CHEMBL123",
                            "name": "Example drug",
                            "maximumClinicalStage": 3,
                            "mechanismsOfAction": {
                                "rows": [{
                                    "mechanismOfAction": "Example mechanism",
                                    "actionType": "INHIBITOR",
                                    "targetName": "LDLR",
                                }]
                            },
                        },
                        "diseases": [],
                        "clinicalReports": [],
                    }],
                },
            }
        }
    }))

    candidates = service.get_known_drugs("ENSG00000130164")

    assert len(candidates) == 1
    assert candidates[0]["evidence_level"] == "target-associated"
    assert candidates[0]["evidence_source"] == "Open Targets Platform"
    assert candidates[0]["drug"]["name"] == "Example drug"
