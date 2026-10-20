from services.evidence_retriever import retrieve_evidence


def sample_records():
    return [
        {
            "evidence_id": "pub-1",
            "source": "PubMed",
            "evidence_type": "publication",
            "gene_symbol": "LDLR",
            "disease_name": "familial hypercholesterolemia",
            "title": "LDLR variants in familial hypercholesterolemia",
            "content": "The abstract discusses LDL receptor function and FH.",
            "source_id": "123456",
            "url": "https://pubmed.ncbi.nlm.nih.gov/123456/",
        },
        {
            "evidence_id": "path-1",
            "source": "Reactome",
            "evidence_type": "pathway",
            "gene_symbol": "APOB",
            "disease_name": None,
            "title": "Lipoprotein metabolism",
            "content": "Pathway information about lipoproteins.",
            "source_id": "R-HSA-123",
            "url": "https://reactome.org/",
        },
    ]


def test_retrieves_relevant_gene_and_disease_evidence():
    results = retrieve_evidence(
        "LDLR familial hypercholesterolemia",
        sample_records(),
    )

    assert results
    assert results[0]["evidence_id"] == "pub-1"
    assert results[0]["retrieval_rank"] == 1
    assert results[0]["url"].startswith("https://pubmed.ncbi.nlm.nih.gov/")


def test_returns_empty_for_unmatched_query():
    assert retrieve_evidence("unrelated astronomy topic", sample_records()) == []


def test_respects_top_k():
    results = retrieve_evidence(
        "lipoprotein LDLR familial hypercholesterolemia metabolism",
        sample_records(),
        top_k=1,
    )

    assert len(results) <= 1


def test_retrieval_does_not_mutate_input_records():
    records = sample_records()
    results = retrieve_evidence("LDLR", records)

    assert results
    assert "retrieval_score" not in records[0]
    assert "retrieval_rank" in results[0]


def test_empty_query_returns_no_results():
    assert retrieve_evidence("", sample_records()) == []
    assert retrieve_evidence("   ", sample_records()) == []
