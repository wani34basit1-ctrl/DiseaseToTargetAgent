from models.disease import Disease
from models.target import Target


def test_disease_creation():
    disease = Disease("Alzheimer disease", "EFO_0000249")
    assert disease.name == "Alzheimer disease"
    assert disease.disease_id == "EFO_0000249"


def test_target_creation():
    target = Target("ENSG00000142192", "APP", "Example target", 0.85)
    assert target.symbol == "APP"
    assert target.association_score == 0.85
    assert target.to_dict()["evidence_source"] == "Open Targets"


def test_target_string():
    target = Target("ENSG00000142192", "APP", "Example target", 0.85)
    assert "APP" in str(target)
