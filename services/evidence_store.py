"""Normalize biological evidence for retrieval-augmented generation."""

from datetime import datetime, timezone
import hashlib
import json


SOURCE_URLS = {
    "Open Targets": "https://platform.opentargets.org/",
    "UniProt": "https://www.uniprot.org/",
    "PubMed": "https://pubmed.ncbi.nlm.nih.gov/",
    "STRING": "https://string-db.org/",
    "Reactome": "https://reactome.org/",
    "GO": "https://geneontology.org/",
    "KEGG": "https://www.kegg.jp/",
}


def _text(value):
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value).strip()


def _record(
    source,
    evidence_type,
    title,
    content,
    disease_name=None,
    gene_symbol=None,
    source_id=None,
    url=None,
    metadata=None,
):
    title = _text(title)
    content = _text(content)

    if not title and not content:
        return None

    identity = "|".join((
        source,
        evidence_type,
        _text(source_id),
        _text(gene_symbol),
        title,
        content,
    ))
    evidence_id = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:20]

    return {
        "evidence_id": evidence_id,
        "source": source,
        "evidence_type": evidence_type,
        "disease_name": _text(disease_name) or None,
        "gene_symbol": _text(gene_symbol) or None,
        "title": title or "Untitled evidence",
        "content": content,
        "source_id": _text(source_id) or None,
        "url": _text(url) or None,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "metadata": metadata or {},
    }


def normalize_target_evidence(target, enrichment, disease_name=None):
    """Transform existing backend records without making network requests."""
    if not isinstance(enrichment, dict):
        enrichment = {}

    records = []
    symbol = _text(getattr(target, "symbol", None)) or None
    target_name = _text(getattr(target, "name", None))
    target_id = _text(getattr(target, "id", None))
    score = getattr(target, "association_score", None)

    if score is not None or target_name or target_id:
        details = []
        if target_name:
            details.append(f"Target name: {target_name}")
        if target_id:
            details.append(f"Target identifier: {target_id}")
        if score is not None:
            details.append(f"Association score: {score}")

        record = _record(
            source="Open Targets",
            evidence_type="target_association",
            title=f"Target association: {symbol or target_name or target_id}",
            content="; ".join(details),
            disease_name=disease_name,
            gene_symbol=symbol,
            source_id=target_id,
            url=(
                f"https://platform.opentargets.org/target/{target_id}"
                if target_id else SOURCE_URLS["Open Targets"]
            ),
            metadata={"association_score": score},
        )
        if record:
            records.append(record)

    for protein in enrichment.get("uniprot_entries", []) or []:
        accession = protein.get("accession")
        name = _text(protein.get("protein_name"))
        record = _record(
            source="UniProt",
            evidence_type="protein",
            title=name or f"UniProt entry {accession or ''}",
            content=f"Protein: {name}; Accession: {accession or 'Not provided'}",
            disease_name=disease_name,
            gene_symbol=symbol,
            source_id=accession,
            url=protein.get("url") or (
                f"https://www.uniprot.org/uniprotkb/{accession}"
                if accession else SOURCE_URLS["UniProt"]
            ),
            metadata={"accession": accession},
        )
        if record:
            records.append(record)

    for paper in enrichment.get("publications", []) or []:
        pmid = paper.get("pmid")
        title = _text(paper.get("title")) or f"PubMed record {pmid or ''}"
        content = title
        if paper.get("abstract"):
            content += "\nAbstract: " + _text(paper["abstract"])

        record = _record(
            source="PubMed",
            evidence_type="publication",
            title=title,
            content=content,
            disease_name=disease_name,
            gene_symbol=symbol,
            source_id=pmid,
            url=paper.get("url") or (
                f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                if pmid else SOURCE_URLS["PubMed"]
            ),
            metadata={
                "pmid": pmid,
                "publication_date": paper.get("pubdate"),
                "abstract_available": bool(paper.get("abstract")),
            },
        )
        if record:
            records.append(record)

    for pathway in enrichment.get("pathways", []) or []:
        pathway_id = pathway.get("pathway_id")
        name = _text(pathway.get("name")) or "Unnamed pathway"
        record = _record(
            source="Reactome",
            evidence_type="pathway",
            title=name,
            content=f"Reactome pathway: {name}; ID: {pathway_id or 'Not provided'}",
            disease_name=disease_name,
            gene_symbol=symbol,
            source_id=pathway_id,
            url=pathway.get("url") or (
                f"https://reactome.org/content/detail/{pathway_id}"
                if pathway_id else SOURCE_URLS["Reactome"]
            ),
            metadata={"pathway_id": pathway_id},
        )
        if record:
            records.append(record)

    network = enrichment.get("string_network", {}) or {}
    for edge in network.get("edges", []) or []:
        protein_a = edge.get("protein_a") or "unknown protein"
        protein_b = edge.get("protein_b") or "unknown protein"
        score_value = edge.get("score")

        record = _record(
            source="STRING",
            evidence_type="protein_association",
            title=f"Protein association: {protein_a} — {protein_b}",
            content=(
                f"STRING reports an association between {protein_a} and {protein_b}. "
                f"Association score: {score_value if score_value is not None else 'not reported'}. "
                "This association does not by itself prove direct physical binding."
            ),
            disease_name=disease_name,
            gene_symbol=symbol,
            source_id=(
                f"{edge.get('protein_a_id', '')}:"
                f"{edge.get('protein_b_id', '')}"
            ),
            url=SOURCE_URLS["STRING"],
            metadata={
                "score": score_value,
                "score_0_1000": edge.get("score_0_1000"),
                "evidence_scores": edge.get("evidence_scores", {}),
                "protein_a_id": edge.get("protein_a_id"),
                "protein_b_id": edge.get("protein_b_id"),
            },
        )
        if record:
            records.append(record)

    categories = (
        ("go_enrichment", "GO", "functional_enrichment"),
        ("kegg_enrichment", "KEGG", "pathway_enrichment"),
    )

    for key, source, evidence_type in categories:
        for term in enrichment.get(key, []) or []:
            term_id = term.get("term_id")
            name = _text(term.get("term_name")) or "Unnamed enrichment term"

            details = [
                f"Term: {name}",
                f"Category: {term.get('category') or 'Not reported'}",
            ]
            if term_id:
                details.append(f"Term ID: {term_id}")
            if term.get("p_value") is not None:
                details.append(f"P-value: {term['p_value']}")
            if term.get("fdr") is not None:
                details.append(f"FDR: {term['fdr']}")

            record = _record(
                source=source,
                evidence_type=evidence_type,
                title=name,
                content="; ".join(details),
                disease_name=disease_name,
                gene_symbol=symbol,
                source_id=term_id,
                url=SOURCE_URLS[source],
                metadata={
                    "term_id": term_id,
                    "category": term.get("category"),
                    "p_value": term.get("p_value"),
                    "fdr": term.get("fdr"),
                    "number_of_genes": term.get("number_of_genes"),
                    "preferred_names": term.get("preferred_names"),
                },
            )
            if record:
                records.append(record)

    return records


def normalize_drug_evidence(
    candidates,
    disease_name=None,
    gene_symbol=None,
):
    """Normalize existing drug records without making network requests."""
    records = []

    for candidate in candidates or []:
        if not isinstance(candidate, dict):
            continue

        drug = candidate.get("drug") or {}
        if not isinstance(drug, dict):
            drug = {}

        name = _text(drug.get("name") or candidate.get("name"))
        drug_id = _text(drug.get("id") or candidate.get("drug_id"))
        if not name and not drug_id:
            continue

        mechanisms = candidate.get("mechanismOfAction") or []
        mechanism_lines = []

        for mechanism in mechanisms:
            if not isinstance(mechanism, dict):
                continue
            description = _text(mechanism.get("mechanism"))
            action = _text(mechanism.get("actionType"))
            target = _text(mechanism.get("targetName"))

            parts = []
            if description:
                parts.append(description)
            if action:
                parts.append(f"Action type: {action}")
            if target:
                parts.append(f"Reported target: {target}")
            if parts:
                mechanism_lines.append("; ".join(parts))

        phase = candidate.get("phase")
        if phase is None:
            phase = drug.get("maximumClinicalStage")

        disease = candidate.get("disease")
        reports = candidate.get("clinicalReports") or []
        report_lines = []

        for report in reports:
            if not isinstance(report, dict):
                continue
            title = _text(report.get("title"))
            url = _text(report.get("url"))
            status = _text(report.get("trialOverallStatus"))
            stage = report.get("clinicalStage")

            parts = []
            if title:
                parts.append(title)
            if status:
                parts.append(f"Trial status: {status}")
            if stage is not None:
                parts.append(f"Clinical stage: {stage}")
            if url:
                parts.append(f"URL: {url}")
            if parts:
                report_lines.append("; ".join(parts))

        content_parts = [
            f"Drug: {name or 'Name unavailable'}",
            f"Drug ID: {drug_id or 'Not reported'}",
            f"Evidence level: {_text(candidate.get('evidence_level')) or 'Not reported'}",
            f"Evidence source: {_text(candidate.get('evidence_source')) or 'Not reported'}",
            f"Clinical stage: {phase if phase is not None else 'Not reported'}",
            f"Status: {_text(candidate.get('status')) or 'Not reported'}",
        ]

        if disease:
            content_parts.append(f"Associated disease: {_text(disease)}")
        if mechanism_lines:
            content_parts.append(
                "Reported mechanism(s): " + " | ".join(mechanism_lines)
            )
        else:
            content_parts.append(
                "Mechanism of action: not supplied in this record"
            )
        if report_lines:
            content_parts.append(
                "Clinical reports: " + " | ".join(report_lines)
            )

        url = (
            f"https://platform.opentargets.org/drug/{drug_id}"
            if drug_id else "https://platform.opentargets.org/"
        )

        record = _record(
            source="Open Targets Platform",
            evidence_type="drug_target_evidence",
            title=name or f"Drug candidate {drug_id}",
            content="; ".join(content_parts),
            disease_name=disease_name,
            gene_symbol=gene_symbol,
            source_id=drug_id or None,
            url=url,
            metadata={
                "drug_id": drug_id or None,
                "evidence_level": candidate.get("evidence_level"),
                "clinical_stage": phase,
                "status": candidate.get("status"),
                "mechanisms": mechanisms,
                "associated_disease": disease,
                "clinical_reports": reports,
            },
        )
        if record:
            records.append(record)

    return records
