"""Controlled research tools for the biological research assistant."""

from agent.target_agent import TargetPrioritizationAgent
from services.open_targets import OpenTargetsService
from services.target_enrichment import TargetEnrichmentService
from services.evidence_store import normalize_target_evidence
import inspect


class ResearchToolError(ValueError):
    """Raised when a research tool request is invalid."""


class ResearchTools:
    """Expose an allowlist of validated biological research operations."""

    MAX_TARGETS = 20
    MAX_PUBLICATIONS = 10

    def __init__(
        self,
        open_targets=None,
        enrichment_service=None,
        target_agent=None,
    ):
        self.open_targets = open_targets or OpenTargetsService()
        self.enrichment_service = (
            enrichment_service or TargetEnrichmentService()
        )
        self.target_agent = target_agent or TargetPrioritizationAgent(
            service=self.open_targets
        )

        self._registry = {
            "search_disease": self.search_disease,
            "find_disease_targets": self.find_disease_targets,
            "search_literature": self.search_literature,
            "enrich_target": self.enrich_target,
        }

    @staticmethod
    def _validate_text(value, name, max_length=200):
        if not isinstance(value, str) or not value.strip():
            raise ResearchToolError(f"{name} must be a non-empty string.")
        value = value.strip()
        if len(value) > max_length:
            raise ResearchToolError(
                f"{name} must be no longer than {max_length} characters."
            )
        return value

    @staticmethod
    def _validate_limit(value, name, maximum):
        if isinstance(value, bool) or not isinstance(value, int):
            raise ResearchToolError(f"{name} must be an integer.")
        if not 1 <= value <= maximum:
            raise ResearchToolError(
                f"{name} must be between 1 and {maximum}."
            )
        return value

    def available_tools(self):
        """Return names and descriptions, not callable objects."""
        return [
            {
                "name": "search_disease",
                "description": "Search Open Targets for matching diseases.",
                "arguments": {"disease_name": "string"},
            },
            {
                "name": "find_disease_targets",
                "description": (
                    "Retrieve disease-associated targets ranked by "
                    "Open Targets association score."
                ),
                "arguments": {
                    "disease_name": "string",
                    "top_n": "integer from 1 to 20",
                },
            },
            {
                "name": "search_literature",
                "description": (
                    "Retrieve PubMed publication metadata for a gene "
                    "and disease."
                ),
                "arguments": {
                    "gene_symbol": "string",
                    "disease_name": "string",
                    "limit": "integer from 1 to 10",
                },
            },
            {
                "name": "enrich_target",
                "description": (
                    "Retrieve available protein, publication, pathway, "
                    "STRING, GO, and KEGG evidence for a gene and disease."
                ),
                "arguments": {
                    "gene_symbol": "string",
                    "disease_name": "string",
                    "limit": "integer from 1 to 10",
                },
            },
        ]

    def execute(self, tool_name, arguments):
        """Execute one allowlisted tool with validated arguments."""
        if not isinstance(tool_name, str) or tool_name not in self._registry:
            raise ResearchToolError("Unknown research tool.")

        if not isinstance(arguments, dict):
            raise ResearchToolError("Tool arguments must be a dictionary.")

        function = self._registry[tool_name]
        parameters = inspect.signature(function).parameters

        allowed = set(parameters)
        supplied = set(arguments)
        unexpected = sorted(supplied - allowed)
        missing = sorted(
            name
            for name, parameter in parameters.items()
            if parameter.default is inspect.Parameter.empty
            and name not in supplied
        )

        if unexpected:
            raise ResearchToolError(
                f"Unexpected argument(s) for {tool_name}: "
                + ", ".join(unexpected)
            )

        if missing:
            raise ResearchToolError(
                f"Missing required argument(s) for {tool_name}: "
                + ", ".join(missing)
            )

        return function(**arguments)

    def search_disease(self, disease_name):
        disease_name = self._validate_text(disease_name, "disease_name")
        hits = self.open_targets.search_disease(disease_name)

        results = []
        for hit in hits[:10]:
            disease_id = hit.get("id")
            url = "https://" + "platform.opentargets.org/"
            if disease_id:
                url += "disease/" + str(disease_id)
            results.append({
                "disease_id": disease_id,
                "name": hit.get("name"),
                "source": "Open Targets",
                "url": url,
            })

        return {
            "tool": "search_disease",
            "query": disease_name,
            "results": results,
        }

    def find_disease_targets(self, disease_name, top_n=10):
        disease_name = self._validate_text(disease_name, "disease_name")
        top_n = self._validate_limit(top_n, "top_n", self.MAX_TARGETS)

        result = self.target_agent.run(disease_name, top_n=top_n)
        disease = result["disease"]

        return {
            "tool": "find_disease_targets",
            "disease": {
                "name": disease.name,
                "disease_id": disease.disease_id,
            },
            "total_associations": result["total_associations"],
            "targets": [target.to_dict() for target in result["targets"]],
            "selection_note": (
                "Targets are ranked by Open Targets association score. "
                "This score is not a probability of clinical efficacy."
            ),
            "source": "Open Targets",
        }

    def search_literature(self, gene_symbol, disease_name, limit=5):
        gene_symbol = self._validate_text(
            gene_symbol, "gene_symbol", max_length=40
        )
        disease_name = self._validate_text(disease_name, "disease_name")
        limit = self._validate_limit(
            limit, "limit", self.MAX_PUBLICATIONS
        )

        publications = self.enrichment_service.get_pubmed_ids(
            gene_symbol, disease_name, limit=limit
        )

        return {
            "tool": "search_literature",
            "gene_symbol": gene_symbol,
            "disease_name": disease_name,
            "source": "PubMed",
            "results": publications[:limit],
            "result_count": len(publications[:limit]),
        }

    def enrich_target(self, gene_symbol, disease_name, limit=5):
        gene_symbol = self._validate_text(
            gene_symbol, "gene_symbol", max_length=40
        )
        disease_name = self._validate_text(disease_name, "disease_name")
        self._validate_limit(limit, "limit", self.MAX_PUBLICATIONS)

        # Resolve the target through Open Targets to obtain its stable ID
        # and association score before normalizing enrichment evidence.
        result = self.target_agent.run(disease_name, top_n=self.MAX_TARGETS)
        target = next(
            (
                item for item in result["targets"]
                if item.symbol.casefold() == gene_symbol.casefold()
            ),
            None,
        )

        if target is None:
            return {
                "tool": "enrich_target",
                "gene_symbol": gene_symbol,
                "disease_name": disease_name,
                "status": "not_found_in_retrieved_candidates",
                "resolution_scope": "retrieved_candidates_only",
                "candidate_count": len(result["targets"]),
                "interpretation": (
                    "The gene was not found among the retrieved disease-associated "
                    "target candidates. This does not establish that the gene has "
                    "no association with the disease."
                ),
                "evidence": [],
                "evidence_status": {},
                "enrichment_errors": [],
            }

        enrichment = self.enrichment_service.enrich_target(
            target, disease_name
        )
        evidence = normalize_target_evidence(
            target, enrichment, disease_name=disease_name
        )

        return {
            "tool": "enrich_target",
            "gene_symbol": target.symbol,
            "disease_name": disease_name,
            "target": target.to_dict(),
            "status": "completed",
            "evidence": evidence,
            "evidence_count": len(evidence),
            "evidence_status": enrichment.get("evidence_status", {}),
            "enrichment_errors": enrichment.get("enrichment_errors", []),
            "limitations": [
                "Target association scores are not probabilities of clinical efficacy.",
                "STRING associations do not by themselves prove direct physical binding.",
                "Missing or failed source results must not be interpreted as negative biological evidence.",
            ],
        }
