from agent.target_agent import TargetPrioritizationAgent


class FakeService:
    def search_disease(self, name):
        return [{"id": "EFO_0000249", "name": "Alzheimer disease"}]

    def get_disease_targets(self, disease_id, page_size=20):
        return {
            "total_associations": 2,
            "rows": [
                {
                    "score": 0.4,
                    "target": {
                        "id": "ENSG00000142192",
                        "approvedSymbol": "APP",
                        "approvedName": "Example A",
                    },
                },
                {
                    "score": 0.9,
                    "target": {
                        "id": "ENSG00000141510",
                        "approvedSymbol": "TP53",
                        "approvedName": "Example B",
                    },
                },
            ],
        }


def test_targets_are_ranked():
    agent = TargetPrioritizationAgent(service=FakeService())
    result = agent.run("Alzheimer disease", 2)
    assert result["targets"][0].symbol == "TP53"
    assert result["targets"][1].symbol == "APP"


def test_top_n():
    agent = TargetPrioritizationAgent(service=FakeService())
    result = agent.run("Alzheimer disease", 1)
    assert len(result["targets"]) == 1


def test_empty_name():
    agent = TargetPrioritizationAgent(service=FakeService())
    try:
        agent.run("   ")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError")
