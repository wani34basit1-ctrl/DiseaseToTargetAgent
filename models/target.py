class Target:
    """Represents a potential therapeutic target."""

    def __init__(
        self,
        gene_id,
        symbol,
        name,
        association_score,
        evidence_source="Open Targets",
    ):
        self.gene_id = gene_id
        self.symbol = symbol
        self.name = name
        self.association_score = float(association_score)
        self.evidence_source = evidence_source

    def to_dict(self):
        return {
            "gene_id": self.gene_id,
            "symbol": self.symbol,
            "name": self.name,
            "association_score": self.association_score,
            "evidence_source": self.evidence_source,
        }

    def __str__(self):
        return (
            f"{self.symbol} — {self.name} "
            f"(association score: {self.association_score:.4f})"
        )
