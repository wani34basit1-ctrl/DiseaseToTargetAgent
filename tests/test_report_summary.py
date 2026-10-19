from types import SimpleNamespace

from main import build_target_summary, build_evidence_summary


def test_target_summary_identifies_top_ranked_target():
    first = SimpleNamespace(
        symbol="LDLR",
        name="low density lipoprotein receptor",
        gene_id="ENSG00000130164",
        association_score=0.9123,
    )
    second = SimpleNamespace(
        symbol="APOB",
        name="apolipoprotein B",
        gene_id="ENSG00000084674",
        association_score=0.85,
    )

    summary = build_target_summary(first, [first, second])

    assert summary["gene_symbol"] == "LDLR"
    assert summary["rank"] == 1
    assert summary["ranked_target_count"] == 2
    assert summary["association_score"] == 0.9123
    assert "not a probability" in summary["interpretation"]


def test_target_summary_handles_target_not_in_ranked_list():
    target = SimpleNamespace(
        symbol="LDLR",
        name="low density lipoprotein receptor",
        gene_id="ENSG00000130164",
        association_score=0.9,
    )

    summary = build_target_summary(target, [])

    assert summary["rank"] is None
    assert summary["ranked_target_count"] == 0


def test_evidence_summary_counts_retrieved_empty_and_failed():
    enrichment = {
        "evidence_status": {
            "uniprot": {
                "status": "retrieved",
                "source": "UniProt",
                "record_count": 1,
            },
            "publications": {
                "status": "failed",
                "source": "PubMed",
                "record_count": None,
            },
        }
    }

    summary = build_evidence_summary(
        enrichment,
        drugs=[],
        treatments=[{"drug": "Atorvastatin"}],
        drug_status={
            "status": "empty",
            "source": "Open Targets Platform",
            "record_count": 0,
        },
    )

    assert summary["status_counts"] == {
        "retrieved": 2,
        "empty": 1,
        "failed": 1,
    }
    assert summary["categories"]["known_drugs"]["status"] == "empty"
    assert summary["categories"]["disease_treatments"]["status"] == "retrieved"
    assert "category_interpretations" in summary
