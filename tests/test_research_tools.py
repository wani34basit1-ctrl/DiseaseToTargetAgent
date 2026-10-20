import pytest

from agent.research_tools import ResearchToolError, ResearchTools


class FakeOpenTargets:
    def search_disease(self, disease_name):
        return [
            {"id": "EFO_0000249", "name": "Alzheimer disease"}
        ]


class FakeTargetAgent:
    def run(self, disease_name, top_n=10):
        from models.disease import Disease
        from models.target import Target

        return {
            "disease": Disease("Alzheimer disease", "EFO_0000249"),
            "total_associations": 1,
            "targets": [
                Target(
                    gene_id="ENSG00000141510",
                    symbol="TP53",
                    name="Example target",
                    association_score=0.8,
                )
            ][:top_n],
            "top_candidate": None,
        }


class FakeEnrichment:
    def get_pubmed_ids(self, gene_symbol, disease_name, limit=5):
        return [
            {
                "pmid": "12345",
                "title": "Example publication",
                "pubdate": "2024",
                "url": "https://pubmed.ncbi.nlm.nih.gov/12345/",
            }
        ][:limit]

    def enrich_target(self, target, disease_name):
        return {
            "uniprot_entries": [],
            "publications": [],
            "pathways": [],
            "string_network": {"nodes": [], "edges": []},
            "go_enrichment": [],
            "kegg_enrichment": [],
            "evidence_status": {
                "publications": {
                    "status": "empty",
                    "source": "PubMed",
                    "record_count": 0,
                }
            },
            "enrichment_errors": [],
        }


@pytest.fixture
def tools():
    return ResearchTools(
        open_targets=FakeOpenTargets(),
        enrichment_service=FakeEnrichment(),
        target_agent=FakeTargetAgent(),
    )


def test_registry_exposes_four_approved_tools(tools):
    names = {item["name"] for item in tools.available_tools()}
    assert names == {
        "search_disease",
        "find_disease_targets",
        "search_literature",
        "enrich_target",
    }


def test_unknown_tool_is_rejected(tools):
    with pytest.raises(ResearchToolError):
        tools.execute("run_arbitrary_code", {})


def test_non_dictionary_arguments_are_rejected(tools):
    with pytest.raises(ResearchToolError):
        tools.execute("search_disease", ["Alzheimer disease"])


def test_empty_disease_name_is_rejected(tools):
    with pytest.raises(ResearchToolError):
        tools.execute("search_disease", {"disease_name": "  "})


def test_search_disease_returns_source_metadata(tools):
    result = tools.execute(
        "search_disease",
        {"disease_name": "Alzheimer disease"},
    )
    assert result["results"][0]["disease_id"] == "EFO_0000249"
    assert result["results"][0]["source"] == "Open Targets"


def test_target_results_are_json_compatible(tools):
    result = tools.execute(
        "find_disease_targets",
        {"disease_name": "Alzheimer disease", "top_n": 1},
    )
    assert result["targets"][0]["symbol"] == "TP53"
    assert "not a probability of clinical efficacy" in result["selection_note"]


def test_invalid_target_limit_is_rejected(tools):
    with pytest.raises(ResearchToolError):
        tools.execute(
            "find_disease_targets",
            {"disease_name": "Alzheimer disease", "top_n": 1000},
        )


def test_literature_search_returns_publication_metadata(tools):
    result = tools.execute(
        "search_literature",
        {
            "gene_symbol": "TP53",
            "disease_name": "Alzheimer disease",
            "limit": 1,
        },
    )
    assert result["result_count"] == 1
    assert result["results"][0]["pmid"] == "12345"


def test_enrichment_preserves_evidence_status(tools):
    result = tools.execute(
        "enrich_target",
        {
            "gene_symbol": "TP53",
            "disease_name": "Alzheimer disease",
        },
    )
    assert result["status"] == "completed"
    assert "evidence_status" in result
    assert "enrichment_errors" in result


def test_enrichment_reports_target_not_found(tools):
    result = tools.execute(
        "enrich_target",
        {
            "gene_symbol": "BRCA1",
            "disease_name": "Alzheimer disease",
        },
    )
    assert result["status"] == "not_found_in_retrieved_candidates"
    assert result["evidence"] == []

def test_rejects_unexpected_tool_arguments(tools):
    with pytest.raises(ResearchToolError, match="Unexpected argument"):
        tools.execute(
            "search_disease",
            {
                "disease_name": "Alzheimer disease",
                "unexpected": True,
            },
        )


def test_rejects_missing_required_tool_arguments(tools):
    with pytest.raises(ResearchToolError, match="Missing required argument"):
        tools.execute("search_disease", {})



def test_missing_target_reports_resolution_scope_and_candidate_count(tools):
    result = tools.execute(
        "enrich_target",
        {
            "gene_symbol": "BRCA1",
            "disease_name": "Alzheimer disease",
        },
    )
    assert result["status"] == "not_found_in_retrieved_candidates"
    assert result["resolution_scope"] == "retrieved_candidates_only"
    assert result["candidate_count"] == 1
    assert "does not establish" in result["interpretation"]


def test_missing_target_does_not_claim_no_disease_association(tools):
    result = tools.execute(
        "enrich_target",
        {
            "gene_symbol": "BRCA1",
            "disease_name": "Alzheimer disease",
        },
    )
    assert "does not establish" in result["interpretation"].lower()
    assert result["evidence"] == []
