from models.disease import Disease
from models.target import Target
from services.open_targets import OpenTargetsService


class TargetPrioritizationAgent:
    """Retrieves, ranks, and selects the top disease-associated target."""

    def __init__(self, service=None):
        self.service = service or OpenTargetsService()

    def find_disease(self, disease_name):
        """Resolve a disease name to its database identifier."""

        if not disease_name.strip():
            raise ValueError("Disease name cannot be empty.")

        hits = self.service.search_disease(disease_name)

        if not hits:
            raise ValueError(f"No disease found for: {disease_name}")

        # Select the first matching disease record automatically.
        selected = hits[0]

        return Disease(
            name=selected["name"],
            disease_id=selected["id"],
        )

    def run(self, disease_name, top_n=10):
        """Retrieve candidates and select the highest-scoring target."""

        if top_n < 1:
            raise ValueError("top_n must be at least 1.")

        disease = self.find_disease(disease_name)

        result = self.service.get_disease_targets(
            disease_id=disease.disease_id,
            page_size=min(max(top_n, 20), 100),
        )

        targets = []

        for row in result["rows"]:
            target_data = row.get("target") or {}
            score = row.get("score")

            if score is None:
                continue

            gene_id = target_data.get("id")
            symbol = target_data.get("approvedSymbol")

            if not gene_id or not symbol:
                continue

            targets.append(
                Target(
                    gene_id=gene_id,
                    symbol=symbol,
                    name=target_data.get("approvedName") or "Unknown",
                    association_score=score,
                )
            )

        targets.sort(
            key=lambda target: target.association_score,
            reverse=True,
        )

        # Select the highest-scoring candidate automatically.
        top_candidate = targets[0] if targets else None

        return {
            "disease": disease,
            "total_associations": result["total_associations"],
            "targets": targets[:top_n],
            "top_candidate": top_candidate,
        }
