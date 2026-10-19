import argparse
import json
import os

import pandas as pd
import requests

from agent.target_agent import TargetPrioritizationAgent
from services.target_enrichment import TargetEnrichmentService
from services.disease_treatments import get_disease_treatments


def main():
    parser = argparse.ArgumentParser(
        description="Disease2Target: Evidence-Enriched Agent"
    )
    parser.add_argument("--disease", required=True)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument(
        "--output",
        default="results/target_report.csv",
    )

    args = parser.parse_args()

    if args.top < 1:
        parser.error("--top must be at least 1")

    agent = TargetPrioritizationAgent()

    try:
        result = agent.run(args.disease, args.top)
    except (requests.RequestException, RuntimeError, ValueError) as exc:
        parser.exit(1, f"Analysis failed: {exc}\n")

    disease = result["disease"]
    top_candidate = result["top_candidate"]

    print("\n" + "=" * 65)
    print("DISEASE2TARGET AGENT")
    print("=" * 65)
    print(f"Disease: {disease.name}")
    print(f"Disease ID: {disease.disease_id}")

    if top_candidate is None:
        print("No target candidates were returned.")
        return

    print("\nTOP PRIORITIZED TARGET")
    print("-" * 65)
    print(f"Gene: {top_candidate.symbol}")
    print(f"Ensembl ID: {top_candidate.gene_id}")
    print(f"Name: {top_candidate.name}")
    print(f"Association score: {top_candidate.association_score:.4f}")

    enrichment_service = TargetEnrichmentService()

    # Keep the core target result even if an enrichment API fails.
    try:
        enrichment = enrichment_service.enrich_target(
            top_candidate,
            disease.name,
        )
    except (requests.RequestException, RuntimeError, ValueError) as exc:
        print(f"\nWarning: protein/literature enrichment failed: {exc}")
        enrichment = {
            "uniprot_entries": [],
            "publications": [],
            "pathways": [],
            "string_network": {"nodes": [], "edges": []},
            "go_enrichment": [],
            "kegg_enrichment": [],
            "enrichment_errors": [{"service": "target enrichment", "error": str(exc)}],
        }

    print("\nUNIPROT PROTEIN INFORMATION")
    print("-" * 65)

    for protein in enrichment["uniprot_entries"]:
        print(f"Accession: {protein['accession']}")
        print(f"Protein: {protein['protein_name']}")
        print(f"Link: {protein['url']}")

    if not enrichment["uniprot_entries"]:
        print("No reviewed UniProt entry found.")

    print("\nPUBMED — LATEST MATCHING PUBLICATIONS")
    print("-" * 65)

    for i, publication in enumerate(enrichment["publications"], start=1):
        print(f"\n{i}. {publication['title']}")
        print(f"   PMID: {publication['pmid']}")
        print(f"   Published: {publication['pubdate']}")
        print(f"   Link: {publication['url']}")

    if not enrichment["publications"]:
        print("No matching PubMed records found.")

    print("\nREACTOME BIOLOGICAL PATHWAYS")
    print("-" * 65)

    for pathway in enrichment["pathways"]:
        print(f"{pathway['pathway_id']}: {pathway['name']}")
        print(f"Link: {pathway['url']}")

    if not enrichment["pathways"]:
        print("No pathway records found for the selected protein.")


    print("\nSTRING PROTEIN-ASSOCIATION NETWORK")
    print("-" * 65)
    network = enrichment.get("string_network") or {"nodes": [], "edges": []}
    print(f"Network nodes: {len(network.get('nodes', []))}")
    print(f"Network edges: {len(network.get('edges', []))}")
    for edge in network.get("edges", [])[:10]:
        score = edge.get("score_0_1000")
        print(
            f"{edge.get('protein_a', 'Unknown')} -- "
            f"{edge.get('protein_b', 'Unknown')} | "
            f"confidence score: {score if score is not None else 'Unavailable'}/1000"
        )
    if not network.get("edges"):
        print("No STRING protein-association edges returned.")

    print("\nGENE ONTOLOGY (GO) ENRICHMENT")
    print("-" * 65)
    go_terms = enrichment.get("go_enrichment", [])
    for term in go_terms[:15]:
        print(
            f"{term.get('category', 'GO')}: {term.get('term_id', 'ID unavailable')} — "
            f"{term.get('term_name', 'Name unavailable')} | FDR: {term.get('fdr', 'Unavailable')}"
        )
    if not go_terms:
        print("No GO enrichment terms returned.")

    print("\nKEGG PATHWAY ENRICHMENT")
    print("-" * 65)
    kegg_terms = enrichment.get("kegg_enrichment", [])
    for term in kegg_terms[:15]:
        print(
            f"{term.get('term_id', 'ID unavailable')} — "
            f"{term.get('term_name', 'Name unavailable')} | FDR: {term.get('fdr', 'Unavailable')}"
        )
    if not kegg_terms:
        print("No KEGG enrichment terms returned.")

    for issue in enrichment.get("enrichment_errors", []):
        print(f"Enrichment warning ({issue.get('service', 'service')}): {issue.get('error', 'Unknown error')}")

    print("\nKNOWN DRUGS / CLINICAL CANDIDATES")
    print("-" * 65)

    try:
        drugs = enrichment_service.get_known_drugs(
            top_candidate.gene_id
        )

        if not drugs:
            print("No known drug records returned for this target.")

        for record in drugs[:10]:
            drug = record.get("drug") or {}
            indication = record.get("disease") or {}

            print(f"Drug: {drug.get('name', 'Name unavailable')}")
            print(f"Drug ID: {drug.get('id', 'Unavailable')}")
            print(f"Approved: {drug.get('isApproved', 'Unknown')}")
            print(f"Clinical phase: {record.get('phase', 'Unknown')}")
            print(f"Status: {record.get('status', 'Unknown')}")
            print(
                "Associated indication: "
                f"{indication.get('name', 'Unavailable')}"
            )
            print("-" * 30)

    except (requests.RequestException, RuntimeError, ValueError) as exc:
        drugs = []
        print(f"Drug information unavailable: {exc}")

    # Display disease-level treatments separately from target-associated drugs.
    treatments = get_disease_treatments(disease.name)

    print("\\nDISEASE-LEVEL TREATMENT OPTIONS")
    print("-" * 65)

    if treatments:
        for treatment in treatments:
            print(f"Drug: {treatment['drug']}")
            print(f"Class: {treatment['class']}")
            print(f"Mechanism: {treatment['mechanism']}")
            print(f"LDLR dependence: {treatment['ldlr_dependence']}")
            print(f"Clinical context: {treatment['scope']}")
            print("-" * 30)
    else:
        print("No curated treatment options available for this disease.")

    # Save the ranked target list.
    rows = []

    for rank, target in enumerate(result["targets"], start=1):
        rows.append({
            "rank": rank,
            "disease": disease.name,
            "disease_id": disease.disease_id,
            **target.to_dict(),
            "selected_top_candidate": (
                target.gene_id == top_candidate.gene_id
            ),
        })

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    columns = [
        "rank", "disease", "disease_id", "gene_id",
        "symbol", "name", "association_score",
        "evidence_source", "selected_top_candidate",
    ]

    pd.DataFrame(rows, columns=columns).to_csv(
        args.output, index=False
    )

    # Save detailed enrichment separately as JSON.
    detail_path = os.path.splitext(args.output)[0] + "_details.json"

    details = {
        "disease_treatments": treatments,
        "disease": disease.to_dict(),
        "selected_target": top_candidate.to_dict(),
        "uniprot": enrichment["uniprot_entries"],
        "publications": enrichment["publications"],
        "pathways": enrichment["pathways"],
        "string_network": enrichment.get("string_network", {"nodes": [], "edges": []}),
        "go_enrichment": enrichment.get("go_enrichment", []),
        "kegg_enrichment": enrichment.get("kegg_enrichment", []),
        "enrichment_errors": enrichment.get("enrichment_errors", []),
        "known_drugs": drugs,
        "interpretation": (
            "Association scores are evidence-ranking measures, "
            "not probabilities of clinical efficacy. Reactome "
            "pathways describe biological context; they do not "
            "necessarily represent a drug's mechanism of action."
        ),
    }

    with open(detail_path, "w", encoding="utf-8") as file:
        json.dump(details, file, indent=2, ensure_ascii=False)

    output_stem = os.path.splitext(args.output)[0]
    network = enrichment.get("string_network") or {"nodes": [], "edges": []}
    export_tables = {
        "_string_nodes.csv": network.get("nodes", []),
        "_string_edges.csv": network.get("edges", []),
        "_go_enrichment.csv": enrichment.get("go_enrichment", []),
        "_kegg_enrichment.csv": enrichment.get("kegg_enrichment", []),
        "_reactome_pathways.csv": enrichment.get("pathways", []),
    }
    for suffix, records in export_tables.items():
        table_path = output_stem + suffix
        pd.DataFrame(records).to_csv(table_path, index=False)
        print(f"Additional results saved to: {table_path}")

    print(f"\nRanked targets saved to: {args.output}")
    print(f"Enriched details saved to: {detail_path}")


if __name__ == "__main__":
    main()
