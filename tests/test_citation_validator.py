import pytest

from services.citation_validator import validate_citations


def test_accepts_valid_citation_labels():
    result = validate_citations(
        "The evidence reports an association [E1].",
        [{"label": "E1", "source": "example"}],
    )
    assert result["is_valid"] is True
    assert result["valid_used_labels"] == ["E1"]
    assert result["invalid_labels"] == []


def test_flags_fabricated_citation_label():
    result = validate_citations(
        "This finding is supported by [E99].",
        [{"label": "E1"}],
    )
    assert result["is_valid"] is False
    assert result["invalid_labels"] == ["E99"]


def test_identifies_answer_without_citations():
    result = validate_citations(
        "The evidence suggests a possible association.",
        [{"label": "E1"}],
    )
    assert result["has_uncited_answer"] is True
    assert result["has_citations"] is False


def test_identifies_unused_evidence_records():
    result = validate_citations(
        "The finding is reported in [E1].",
        [{"label": "E1"}, {"label": "E2"}],
    )
    assert result["unused_labels"] == ["E2"]


def test_accepts_multiple_valid_labels():
    result = validate_citations(
        "Two records report relevant findings [E1] [E2].",
        [{"label": "E1"}, {"label": "E2"}],
    )
    assert result["is_valid"] is True
    assert result["valid_used_labels"] == ["E1", "E2"]


def test_empty_answer_has_no_uncited_claim():
    result = validate_citations("", [{"label": "E1"}])
    assert result["has_uncited_answer"] is False
    assert result["has_citations"] is False


def test_rejects_invalid_input_types():
    with pytest.raises(TypeError):
        validate_citations(None, [])

    with pytest.raises(TypeError):
        validate_citations("Answer [E1]", "E1")


def test_reports_validation_scope():
    result = validate_citations(
        "Finding [E1].",
        [{"label": "E1"}],
    )
    assert result["validation_scope"] == "label_integrity_only"
