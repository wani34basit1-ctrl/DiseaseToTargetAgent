import re
import requests

from datetime import date, datetime


class TargetEnrichmentService:
    """Adds protein, publication, and pathway evidence to a target."""

    UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
    PUBMED_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    REACTOME_URL = "https://reactome.org/ContentService/data"

    def __init__(self, timeout=25):
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Disease2Target/1.0 (educational research project)"
        })

    def get_uniprot(self, gene_symbol):
        """Find reviewed human UniProt entries for a gene symbol."""

        response = self.session.get(
            self.UNIPROT_URL,
            params={
                "query": (
                    f"(gene_exact:{gene_symbol}) AND "
                    "(organism_id:9606) AND (reviewed:true)"
                ),
                "fields": "accession,protein_name,gene_names",
                "format": "json",
                "size": 5,
            },
            timeout=self.timeout,
        )
        response.raise_for_status()

        results = response.json().get("results", [])

        if not results:
            return []

        proteins = []

        for entry in results:
            description = entry.get("proteinDescription", {})
            recommended = description.get("recommendedName", {})
            protein_name = recommended.get("fullName", {}).get("value")

            proteins.append({
                "accession": entry.get("primaryAccession"),
                "protein_name": protein_name or "Name unavailable",
                "url": (
                    "https://www.uniprot.org/uniprotkb/"
                    f"{entry.get('primaryAccession')}"
                ),
            })

        return proteins

    def get_pubmed_ids(self, gene_symbol, disease_name, limit=5):
        """Return matching PubMed papers with validated publication dates."""

        if not gene_symbol or not disease_name:
            return []

        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

        query = (
            f'({gene_symbol}[Title/Abstract]) AND '
            f'("{disease_name}"[Title/Abstract])'
        )

        response = self.session.get(
            search_url,
            params={
                "db": "pubmed",
                "term": query,
                "retmax": max(limit * 3, limit),
                "retmode": "json",
                "sort": "pub date",
                "datetype": "pdat",
                "mindate": "1990/01/01",
                "maxdate": date.today().strftime("%Y/%m/%d"),
            },
            timeout=self.timeout,
        )
        response.raise_for_status()

        ids = response.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            return []

        response = self.session.get(
            summary_url,
            params={
                "db": "pubmed",
                "id": ",".join(ids),
                "retmode": "json",
            },
            timeout=self.timeout,
        )
        response.raise_for_status()

        summaries = response.json().get("result", {})
        papers = []
        today = date.today()

        for pmid in ids:
            item = summaries.get(pmid, {})
            if not item:
                continue

            pubdate = item.get("pubdate", "")
            match = re.match(
                r"^\s*(\d{4})(?:\s+([A-Za-z]{3,9}))?(?:\s+(\d{1,2}))?",
                pubdate,
            )

            if not match:
                continue

            year = int(match.group(1))
            month_text = match.group(2)
            day_text = match.group(3)

            if year > today.year:
                continue

            if year == today.year:
                if not month_text:
                    continue

                try:
                    month = datetime.strptime(
                        month_text[:3], "%b"
                    ).month
                except ValueError:
                    continue

                if month > today.month:
                    continue

                if month == today.month:
                    if not day_text or int(day_text) > today.day:
                        continue

            papers.append({
                "pmid": pmid,
                "title": item.get("title", "Title unavailable"),
                "pubdate": pubdate,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })

            if len(papers) >= limit:
                break

        return papers


    def get_reactome_pathways(self, uniprot_accession, limit=10):
        """Retrieve pathways associated with a UniProt accession."""

        if not uniprot_accession:
            return []

        url = (
            f"{self.REACTOME_URL}/mapping/UniProt/"
            f"{uniprot_accession}/pathways"
        )

        response = self.session.get(
            url,
            headers={"Accept": "application/json"},
            timeout=self.timeout,
        )
        if response.status_code == 404:
            return []

        response.raise_for_status()
        data = response.json()

        pathways = []

        for item in data[:limit]:
            pathway_id = item.get("stId")

            pathways.append({
                "pathway_id": pathway_id,
                "name": item.get("displayName", "Name unavailable"),
                "url": (
                    f"https://reactome.org/content/detail/{pathway_id}"
                    if pathway_id else None
                ),
            })

        return pathways


    def get_known_drugs(self, ensembl_id):
        """Retrieve clinical candidates and mechanisms for a target."""

        query = """
        query TargetClinicalCandidates($ensemblId: String!) {
            target(ensemblId: $ensemblId) {
                drugAndClinicalCandidates {
                    count
                    rows {
                        maxClinicalStage
                        drug {
                            id
                            name
                            maximumClinicalStage
                            mechanismsOfAction {
                                rows {
                                    mechanismOfAction
                                    actionType
                                    targetName
                                }
                            }
                        }
                        diseases {
                            diseaseFromSource
                            disease {
                                id
                                name
                            }
                        }
                        clinicalReports {
                            id
                            title
                            url
                            clinicalStage
                            trialOverallStatus
                        }
                    }
                }
            }
        }
        """

        url = "https://api.platform.opentargets.org/api/v4/graphql"

        response = self.session.post(
            url,
            json={
                "query": query,
                "variables": {"ensemblId": ensembl_id},
            },
            timeout=self.timeout,
        )
        response.raise_for_status()

        payload = response.json()

        if payload.get("errors"):
            raise RuntimeError(
                f"Open Targets GraphQL error: {payload['errors']}"
            )

        target_data = (payload.get("data") or {}).get("target")
        if not target_data:
            return []

        candidates = (
            target_data.get("drugAndClinicalCandidates") or {}
        ).get("rows") or []

        results = []

        for row in candidates:
            drug = row.get("drug") or {}
            if not drug:
                continue

            mechanism_data = drug.get("mechanismsOfAction") or {}
            mechanisms = mechanism_data.get("rows") or []
            diseases = row.get("diseases") or []
            reports = row.get("clinicalReports") or []

            results.append({
                "evidence_level": "target-associated",
                "evidence_source": "Open Targets Platform",
                "drug": {
                    "id": drug.get("id"),
                    "name": drug.get("name"),
                    "maximumClinicalStage": drug.get(
                        "maximumClinicalStage"
                    ),
                    "isApproved": None,
                    "yearOfFirstApproval": None,
                    "hasBeenWithdrawn": None,
                },
                "phase": row.get("maxClinicalStage"),
                "status": next(
                    (
                        report.get("trialOverallStatus")
                        for report in reports
                        if report.get("trialOverallStatus")
                    ),
                    "Not reported",
                ),
                "disease": next(
                    (
                        item.get("disease")
                        for item in diseases
                        if item.get("disease")
                    ),
                    None,
                ),
                "mechanismOfAction": [
                    {
                        "mechanism": item.get("mechanismOfAction"),
                        "actionType": item.get("actionType"),
                        "targetName": item.get("targetName"),
                    }
                    for item in mechanisms
                ],
                "clinicalReports": reports,
            })

        return results

    STRING_API_URL = "https://string-db.org/api/json"

    def get_string_network(
        self, gene_symbol, species=9606, required_score=400, limit=10
    ):
        """Return a STRING protein-association network around a gene symbol.

        STRING scores are confidence scores for protein associations, not proof
        of direct physical binding. Scores are retained on a 0-1000 scale.
        """
        if not gene_symbol or not str(gene_symbol).strip():
            return {"nodes": [], "edges": []}

        response = self.session.get(
            f"{self.STRING_API_URL}/network",
            params={
                "identifiers": str(gene_symbol).strip(),
                "species": species,
                "required_score": required_score,
                "limit": limit,
                "caller_identity": "Disease2Target",
            },
            timeout=self.timeout,
        )
        response.raise_for_status()
        records = response.json()

        nodes_by_id = {}
        edges = []
        for item in records:
            protein_a = item.get("stringId_A")
            protein_b = item.get("stringId_B")
            name_a = item.get("preferredName_A")
            name_b = item.get("preferredName_B")
            if protein_a and name_a:
                nodes_by_id[protein_a] = {
                    "string_id": protein_a,
                    "preferred_name": name_a,
                }
            if protein_b and name_b:
                nodes_by_id[protein_b] = {
                    "string_id": protein_b,
                    "preferred_name": name_b,
                }
            if protein_a and protein_b:
                edges.append({
                    "protein_a_id": protein_a,
                    "protein_a": name_a,
                    "protein_b_id": protein_b,
                    "protein_b": name_b,
                    "score": item.get("score"),
                    "score_0_1000": (
                        round(float(item["score"]) * 1000)
                        if item.get("score") is not None else None
                    ),
                    "evidence_scores": {
                        key: item.get(key)
                        for key in ("nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore")
                        if item.get(key) is not None
                    },
                })

        return {
            "nodes": sorted(nodes_by_id.values(), key=lambda node: node["preferred_name"]),
            "edges": edges,
        }

    def get_string_functional_enrichment(self, identifiers, species=9606):
        """Fetch STRING functional enrichment for a list of protein identifiers."""
        unique_ids = list(dict.fromkeys(
            str(identifier).strip() for identifier in identifiers
            if identifier and str(identifier).strip()
        ))
        if not unique_ids:
            return []

        response = self.session.get(
            f"{self.STRING_API_URL}/enrichment",
            params={
                "identifiers": "\r".join(unique_ids),
                "species": species,
                "caller_identity": "Disease2Target",
            },
            timeout=self.timeout,
        )
        response.raise_for_status()
        results = response.json()
        normalized = []
        for item in results:
            category = item.get("category", "")
            normalized.append({
                "category": category,
                "term_id": item.get("term", item.get("termId")),
                "term_name": item.get("description", item.get("termName", "Name unavailable")),
                "number_of_genes": item.get("number_of_genes"),
                "number_of_genes_in_background": item.get("number_of_genes_in_background"),
                "input_gene_count": item.get("inputGenes", item.get("input_gene_count")),
                "preferred_names": item.get("preferredNames", item.get("preferred_names", "")),
                "p_value": item.get("p_value", item.get("pValue")),
                "fdr": item.get("fdr"),
                "source": "STRING",
            })
        return normalized

    def enrich_target(self, target, disease_name):
        """Collect protein, literature, pathway, interaction, and enrichment evidence.

        Failures in the newly added STRING endpoints are isolated so that they
        do not discard existing UniProt, PubMed, or Reactome results.
        """
        proteins = self.get_uniprot(target.symbol)
        publications = self.get_pubmed_ids(target.symbol, disease_name)

        # Use the first reviewed protein entry as the initial mapping.
        # Verify this mapping when a gene has multiple protein isoforms.
        selected_protein = proteins[0] if proteins else None
        pathways = self.get_reactome_pathways(
            selected_protein["accession"] if selected_protein else None
        )

        string_network = {"nodes": [], "edges": []}
        go_enrichment = []
        kegg_enrichment = []
        enrichment_errors = []
        try:
            string_network = self.get_string_network(target.symbol)
        except (requests.RequestException, RuntimeError, ValueError, TypeError) as exc:
            enrichment_errors.append({"service": "STRING network", "error": str(exc)})

        # Enrich the network's proteins; include the seed gene if STRING returns
        # no node for it. STRING's network scores are association confidence,
        # not a claim that every edge is a direct physical interaction.
        identifiers = [node["preferred_name"] for node in string_network["nodes"]]
        if target.symbol and target.symbol not in identifiers:
            identifiers.insert(0, target.symbol)
        try:
            functional_terms = self.get_string_functional_enrichment(identifiers)
            for term in functional_terms:
                category = str(term.get("category") or "").upper()
                if category in {"PROCESS", "GO:BP", "BIOLOGICAL_PROCESS"}:
                    term["category"] = "GO:BP"
                    go_enrichment.append(term)
                elif category in {"FUNCTION", "GO:MF", "MOLECULAR_FUNCTION"}:
                    term["category"] = "GO:MF"
                    go_enrichment.append(term)
                elif category in {"COMPONENT", "GO:CC", "CELLULAR_COMPONENT"}:
                    term["category"] = "GO:CC"
                    go_enrichment.append(term)
                elif category == "KEGG" or category.startswith("KEGG"):
                    term["category"] = "KEGG"
                    kegg_enrichment.append(term)
        except (requests.RequestException, RuntimeError, ValueError, TypeError) as exc:
            enrichment_errors.append({"service": "STRING functional enrichment", "error": str(exc)})

        return {
            "uniprot_entries": proteins,
            "publications": publications,
            "pathways": pathways,
            "string_network": string_network,
            "go_enrichment": go_enrichment,
            "kegg_enrichment": kegg_enrichment,
            "enrichment_errors": enrichment_errors,
        }
