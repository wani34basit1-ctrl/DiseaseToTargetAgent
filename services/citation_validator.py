"""Validate evidence citation labels in research assistant answers."""

import re


CITATION_PATTERN = re.compile(r"\[([Ee]\d+)\]")


def validate_citations(answer, citations):
    """Check cited labels against the evidence records supplied to the answer.

    This validates citation-label integrity, not whether a source supports
    the scientific claim attached to a citation.
    """
    if not isinstance(answer, str):
        raise TypeError("answer must be a string.")

    if not isinstance(citations, list):
        raise TypeError("citations must be a list.")

    valid_labels = set()

    for citation in citations:
        if not isinstance(citation, dict):
            raise TypeError("Each citation must be a dictionary.")

        label = citation.get("label")
        if isinstance(label, str) and re.fullmatch(r"E\d+", label):
            valid_labels.add(label)

    found_labels = {
        match.group(1).upper()
        for match in CITATION_PATTERN.finditer(answer)
    }

    invalid_labels = sorted(found_labels - valid_labels)
    valid_used_labels = sorted(found_labels & valid_labels)
    unused_labels = sorted(valid_labels - found_labels)

    return {
        "is_valid": not invalid_labels,
        "found_labels": sorted(found_labels),
        "valid_used_labels": valid_used_labels,
        "invalid_labels": invalid_labels,
        "unused_labels": unused_labels,
        "has_citations": bool(found_labels),
        "has_uncited_answer": bool(answer.strip()) and not bool(found_labels),
        "validation_scope": "label_integrity_only",
    }
