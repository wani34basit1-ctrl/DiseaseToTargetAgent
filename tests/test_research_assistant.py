from unittest.mock import patch

import pytest

from services.research_assistant import answer_research_question


def evidence():
    return [
        {
            "evidence_id": "pub-ldlr-1",
            "source": "PubMed",
            "evidence_type": "publication",
            "gene_symbol": "LDLR",
            "disease_name": "familial hypercholesterolemia",
            "title": "LDLR and familial hypercholesterolemia",
            "content": "The publication discusses LDLR and familial hypercholesterolemia.",
            "source_id": "123456",
            "url": "https://pubmed.ncbi.nlm.nih.gov/123456/",
        }
    ]


def test_returns_grounded_answer_and_source_citation():
    with patch(
        "services.research_assistant.generate_response",
        return_value="The retrieved publication discusses LDLR and FH [E1].",
    ) as mocked:
        result = answer_research_question(
            "LDLR familial hypercholesterolemia",
            evidence(),
        )

    assert result["retrieved_count"] == 1
    assert result["citations"][0]["label"] == "E1"
    assert result["citations"][0]["url"] == evidence()[0]["url"]
    assert "[E1]" in result["answer"]
    mocked.assert_called_once()


def test_no_matching_evidence_does_not_call_model():
    with patch("services.research_assistant.generate_response") as mocked:
        result = answer_research_question(
            "unrelated astronomy topic",
            evidence(),
        )

    assert result["retrieved_count"] == 0
    assert result["citations"] == []
    assert "insufficient" in result["answer"].lower()
    mocked.assert_not_called()


def test_rejects_empty_question():
    with pytest.raises(ValueError):
        answer_research_question("", evidence())


def test_rejects_non_list_evidence():
    with pytest.raises(TypeError):
        answer_research_question("LDLR research", {})


def test_answer_includes_only_retrieved_source_metadata():
    with patch(
        "services.research_assistant.generate_response",
        return_value="The record reports the relevant topic [E1].",
    ):
        result = answer_research_question(
            "LDLR familial hypercholesterolemia",
            evidence(),
        )

    assert len(result["citations"]) == 1
    assert result["citations"][0]["source_id"] == "123456"


def test_integration_reports_valid_citation():
    with patch(
        "services.research_assistant.generate_response",
        return_value="The publication discusses LDLR and FH [E1].",
    ):
        result = answer_research_question(
            "LDLR familial hypercholesterolemia",
            evidence(),
        )

    validation = result["citation_validation"]
    assert validation["is_valid"] is True
    assert validation["valid_used_labels"] == ["E1"]
    assert validation["invalid_labels"] == []


def test_integration_flags_fabricated_citation():
    with patch(
        "services.research_assistant.generate_response",
        return_value="This claim is supported by [E99].",
    ):
        result = answer_research_question(
            "LDLR familial hypercholesterolemia",
            evidence(),
        )

    validation = result["citation_validation"]
    assert validation["is_valid"] is False
    assert validation["invalid_labels"] == ["E99"]


def test_integration_flags_answer_without_citations():
    with patch(
        "services.research_assistant.generate_response",
        return_value="The publication discusses LDLR and FH.",
    ):
        result = answer_research_question(
            "LDLR familial hypercholesterolemia",
            evidence(),
        )

    validation = result["citation_validation"]
    assert validation["has_uncited_answer"] is True
    assert validation["has_citations"] is False


def test_no_evidence_returns_validation_report_without_calling_model():
    with patch("services.research_assistant.generate_response") as mocked:
        result = answer_research_question(
            "unrelated astronomy topic",
            evidence(),
        )

    assert result["retrieved_count"] == 0
    assert result["citations"] == []
    assert "citation_validation" in result
    assert result["citation_validation"]["has_citations"] is False
    mocked.assert_not_called()
