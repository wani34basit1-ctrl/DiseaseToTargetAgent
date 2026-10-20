"""Deterministic retrieval over normalized biological evidence."""

import re


def _tokens(text):
    """Tokenize text while preserving gene symbols and biological identifiers."""
    return re.findall(r"[a-zA-Z0-9]+(?:[-_.:/][a-zA-Z0-9]+)*", str(text).lower())


def _record_text(record):
    fields = (
        "title",
        "content",
        "source",
        "evidence_type",
        "gene_symbol",
        "disease_name",
        "source_id",
    )
    return " ".join(str(record.get(field) or "") for field in fields)


def retrieve_evidence(query, records, top_k=5, min_score=1):
    """Rank evidence records against a query and return relevant records.

    The score measures textual relevance only; it is not a measure of
    biological validity or clinical confidence.
    """
    if not isinstance(query, str) or not query.strip():
        return []

    if top_k < 1:
        return []

    query_tokens = set(_tokens(query))
    if not query_tokens:
        return []

    ranked = []

    for index, record in enumerate(records or []):
        if not isinstance(record, dict):
            continue

        title_tokens = set(_tokens(record.get("title") or ""))
        content_tokens = set(_tokens(record.get("content") or ""))
        gene_tokens = set(_tokens(record.get("gene_symbol") or ""))
        disease_tokens = set(_tokens(record.get("disease_name") or ""))
        identifier_tokens = set(_tokens(record.get("source_id") or ""))

        all_tokens = set(_tokens(_record_text(record)))
        overlap = query_tokens & all_tokens

        if not overlap:
            continue

        # Exact identifier matches receive extra weight, without implying truth.
        score = float(len(overlap))
        score += 4.0 * len(query_tokens & gene_tokens)
        score += 3.0 * len(query_tokens & disease_tokens)
        score += 2.0 * len(query_tokens & identifier_tokens)
        score += 1.5 * len(query_tokens & title_tokens)

        if score < min_score:
            continue

        result = dict(record)
        result["retrieval_score"] = round(score, 3)
        result["retrieval_rank"] = 0
        ranked.append((score, index, result))

    ranked.sort(key=lambda item: (-item[0], item[1]))

    results = []
    for rank, (_, _, record) in enumerate(ranked[:top_k], start=1):
        record["retrieval_rank"] = rank
        results.append(record)

    return results
