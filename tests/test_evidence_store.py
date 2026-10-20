from types import SimpleNamespace

from services.evidence_store import normalize_target_evidence


def sample_target():
    return SimpleNamespace(
        symbol="LDLR",
        name="Low-density lipoprotein receptor",
        id="ENSG00000130164",
        association_score=0.82,
    )


def test_normalizes_target_and_publication():
    enrichment = {
        "publications": [{
            "pmid": "123456",
            "title": "LDLR and familial hypercholesterolemia",
            "pubdate": "2024",
            "url": "https://pubmed.ncbi.nlm.nih.gov/123456/",
        }]
    }

    records = normalize_target_evidence(
        sample_target(), enrichment, "familial hypercholesterolemia"
    )

    assert len(records) == 2
    assert records[0]["source"] == "Open Targets"
    assert records[0]["metadata"]["association_score"] == 0.82

    paper = records[1]
    assert paper["source"] == "PubMed"
    assert paper["source_id"] == "123456"
    assert paper["url"] == "https://pubmed.ncbi.nlm.nih.gov/123456/"
    assert paper["disease_name"] == "familial hypercholesterolemia"


def test_normalizes_pathway_and_string_evidence():
    enrichment = {
        "pathways": [{
            "pathway_id": "R-HSA-123",
            "name": "Cholesterol metabolism",
        }],
        "string_network": {
            "nodes": [],
            "edges": [{
                "protein_a": "LDLR",
                "protein_a_id": "P01130",
                "protein_b": "APOB",
                "protein_b_id": "P04114",
                "score": 0.9,
                "score_0_1000": 900,
                "evidence_scores": {"dscore": 0.7},
            }],
        },
    }

    records = normalize_target_evidence(sample_target(), enrichment)

    assert len(records) == 3
    by_source = {record["source"]: record for record in records}
    assert set(by_source) == {"Open Targets", "Reactome", "STRING"}
    assert "does not by itself prove direct physical binding" in by_source["STRING"]["content"]


def test_records_have_stable_ids_and_required_fields():
    target = sample_target()
    enrichment = {
        "uniprot_entries": [{
            "accession": "P01130",
            "protein_name": "Low-density lipoprotein receptor",
        }]
    }

    first = normalize_target_evidence(target, enrichment)
    second = normalize_target_evidence(target, enrichment)

    assert first[0]["evidence_id"] == second[0]["evidence_id"]
    assert {
        "evidence_id", "source", "evidence_type", "title",
        "content", "url", "retrieved_at", "metadata",
    }.issubset(first[0])


def test_empty_input_returns_no_records():
    empty_target = SimpleNamespace(
        symbol=None, name=None, id=None, association_score=None
    )
    assert normalize_target_evidence(empty_target, {}) == []
    assert normalize_target_evidence(empty_target, None) == []
