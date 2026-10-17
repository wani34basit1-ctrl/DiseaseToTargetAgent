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

    def enrich_target(self, target, disease_name):
        """Collect supporting protein, literature, and pathway evidence."""

        proteins = self.get_uniprot(target.symbol)
        publications = self.get_pubmed_ids(
            target.symbol, disease_name
        )

        # Use the first reviewed protein entry as the initial mapping.
        # This mapping should be verified if multiple proteins are returned.
        selected_protein = proteins[0] if proteins else None

        pathways = self.get_reactome_pathways(
            selected_protein["accession"] if selected_protein else None
        )

        return {
            "uniprot_entries": proteins,
            "publications": publications,
            "pathways": pathways,
        }
