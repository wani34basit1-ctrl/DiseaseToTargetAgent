from services.agentic_research import AgenticResearchService


class FakeExecutor:
    def __init__(self, results):
        self.results = results

    def run(self, question, model=None):
        return {
            "request": question,
            "planned_call_count": len(self.results),
            "completed_call_count": sum(
                item["status"] == "completed" for item in self.results
            ),
            "results": self.results,
        }


def test_literature_results_become_source_linked_evidence():
    execution = FakeExecutor([
        {
            "tool": "search_literature",
            "status": "completed",
            "result": {
                "gene_symbol": "LDLR",
                "disease_name": "familial hypercholesterolemia",
                "results": [{
                    "pmid": "123456",
                    "title": "LDLR and familial hypercholesterolemia",
                    "pubdate": "2024",
                    "url": "https://pubmed.ncbi.nlm.nih.gov/123456/",
                }],
            },
        }
    ])
    service = AgenticResearchService(
        executor=execution,
        answerer=lambda question, evidence: {
            "answer": "Evidence retrieved." if evidence else "Insufficient evidence.",
            "citations": [],
        },
    )

    result = service.run("Research LDLR in familial hypercholesterolemia")

    assert result["evidence_count"] == 1
    assert result["evidence"][0]["source"] == "PubMed"
    assert result["evidence"][0]["source_id"] == "123456"
    assert result["evidence"][0]["url"].endswith("/123456/")


def test_failed_tool_results_are_not_treated_as_evidence():
    execution = FakeExecutor([
        {
            "tool": "search_literature",
            "status": "failed",
            "error": "Database unavailable",
        }
    ])
    service = AgenticResearchService(
        executor=execution,
        answerer=lambda question, evidence: {
            "answer": "Insufficient evidence.",
            "evidence_count": len(evidence),
        },
    )

    result = service.run("Research LDLR")

    assert result["evidence"] == []
    assert result["evidence_count"] == 0
    assert result["research_answer"]["evidence_count"] == 0


def test_enrichment_evidence_preserves_original_record():
    record = {
        "evidence_id": "evidence-1",
        "source": "UniProt",
        "evidence_type": "protein",
        "gene_symbol": "LDLR",
        "disease_name": "familial hypercholesterolemia",
        "title": "LDLR protein",
        "content": "Protein evidence",
        "source_id": "P01130",
        "url": "https://www.uniprot.org/",
        "metadata": {"accession": "P01130"},
    }
    execution = FakeExecutor([
        {
            "tool": "enrich_target",
            "status": "completed",
            "result": {"evidence": [record]},
        }
    ])
    service = AgenticResearchService(
        executor=execution,
        answerer=lambda question, evidence: {"answer": "Reviewed."},
    )

    result = service.run("Review LDLR protein evidence")

    assert result["evidence"] == [record]
    assert result["evidence"][0]["source_id"] == "P01130"


def test_empty_evidence_is_passed_to_answerer():
    execution = FakeExecutor([])
    captured = {}
    def answerer(question, evidence):
        captured["evidence"] = evidence
        return {"answer": "Insufficient evidence."}

    service = AgenticResearchService(
        executor=execution,
        answerer=answerer,
    )
    result = service.run("Research a gene")

    assert captured["evidence"] == []
    assert result["research_answer"]["answer"] == "Insufficient evidence."


def test_duplicate_evidence_ids_are_deduplicated():
    record = {
        "evidence_id": "same-id",
        "source": "PubMed",
        "evidence_type": "publication",
        "title": "Paper",
        "content": "Evidence",
    }
    execution = FakeExecutor([
        {
            "tool": "enrich_target",
            "status": "completed",
            "result": {"evidence": [record]},
        },
        {
            "tool": "enrich_target",
            "status": "completed",
            "result": {"evidence": [record]},
        },
    ])
    service = AgenticResearchService(
        executor=execution,
        answerer=lambda question, evidence: {"answer": "Reviewed."},
    )

    result = service.run("Review evidence")

    assert result["evidence_count"] == 1


def test_rejects_empty_question():
    service = AgenticResearchService(
        executor=FakeExecutor([]),
        answerer=lambda question, evidence: {},
    )
    try:
        service.run("  ")
    except ValueError as exc:
        assert "research question" in str(exc).lower()
    else:
        raise AssertionError("Expected empty question to be rejected")


def test_failed_tool_evidence_is_excluded_even_if_payload_contains_records():
    record = {
        "evidence_id": "failed-evidence",
        "source": "PubMed",
        "evidence_type": "publication",
        "title": "Unverified paper",
        "content": "This must not reach the answerer.",
    }
    execution = FakeExecutor([
        {
            "tool": "enrich_target",
            "status": "failed",
            "error": "Tool failed",
            "result": {"evidence": [record]},
        }
    ])
    captured = {}

    def answerer(question, evidence):
        captured["evidence"] = evidence
        return {"answer": "Insufficient evidence."}

    service = AgenticResearchService(
        executor=execution,
        answerer=answerer,
    )

    result = service.run("Research LDLR")

    assert captured["evidence"] == []
    assert result["evidence"] == []
    assert result["evidence_count"] == 0


def test_rejects_malformed_executor_response():
    class MalformedExecutor:
        def run(self, question, model=None):
            return {"results": "not-a-list"}

    service = AgenticResearchService(
        executor=MalformedExecutor(),
        answerer=lambda question, evidence: {},
    )

    try:
        service.run("Research LDLR")
    except ValueError as exc:
        assert "results list" in str(exc)
    else:
        raise AssertionError("Expected malformed executor response to be rejected")


def test_reports_citation_validation_scope():
    service = AgenticResearchService(
        executor=FakeExecutor([]),
        answerer=lambda question, evidence: {
            "answer": "Evidence reviewed.",
            "citation_validation": {
                "validation_scope": "label_integrity_only",
            },
        },
    )

    result = service.run("Research LDLR")

    assert result["citation_validation_scope"] == "label_integrity_only"
    assert result["citation_support_verified"] is False


def test_rejects_malformed_individual_execution_result():
    class MalformedExecutor:
        def run(self, question, model=None):
            return {"results": ["not-a-result-dictionary"]}

    service = AgenticResearchService(
        executor=MalformedExecutor(),
        answerer=lambda question, evidence: {},
    )

    try:
        service.run("Research LDLR")
    except ValueError as exc:
        assert "malformed result" in str(exc)
    else:
        raise AssertionError("Expected malformed result to be rejected")
