class Disease:
    """Represents a disease."""

    def __init__(self, name, disease_id):
        self.name = name
        self.disease_id = disease_id

    def to_dict(self):
        return {
            "name": self.name,
            "disease_id": self.disease_id,
        }

    def __str__(self):
        return f"{self.name} ({self.disease_id})"
