import requests


class OpenTargetsService:
    """Retrieves disease and target information from Open Targets."""

    API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

    SEARCH_QUERY = """
    query SearchDisease($queryString: String!) {
        search(queryString: $queryString, entityNames: ["disease"]) {
            hits {
                id
                name
            }
        }
    }
    """

    TARGETS_QUERY = """
    query DiseaseTargets($efoId: String!, $index: Int!, $size: Int!) {
        disease(efoId: $efoId) {
            id
            name
            associatedTargets(page: {index: $index, size: $size}) {
                count
                rows {
                    score
                    target {
                        id
                        approvedSymbol
                        approvedName
                    }
                }
            }
        }
    }
    """

    def __init__(self, timeout=30):
        self.timeout = timeout
        self.session = requests.Session()

    def _post_query(self, query, variables):
        response = self.session.post(
            self.API_URL,
            json={"query": query, "variables": variables},
            timeout=self.timeout,
        )
        response.raise_for_status()
        payload = response.json()

        if payload.get("errors"):
            raise RuntimeError(
                f"Open Targets API error: {payload['errors']}"
            )

        if "data" not in payload:
            raise RuntimeError("The API returned no data.")

        return payload["data"]

    def search_disease(self, disease_name):
        data = self._post_query(
            self.SEARCH_QUERY,
            {"queryString": disease_name},
        )
        search_data = data.get("search") or {}
        return search_data.get("hits") or []

    def get_disease_targets(self, disease_id, page_size=20):
        data = self._post_query(
            self.TARGETS_QUERY,
            {
                "efoId": disease_id,
                "index": 0,
                "size": page_size,
            },
        )

        disease = data.get("disease")

        if disease is None:
            raise ValueError(f"No disease found for ID: {disease_id}")

        associations = disease.get("associatedTargets") or {}

        return {
            "disease": disease,
            "total_associations": associations.get("count", 0),
            "rows": associations.get("rows") or [],
        }
