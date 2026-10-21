import json
import os
import math
from html import escape
from datetime import datetime

import gradio as gr
import pandas as pd
import requests

from agent.target_agent import TargetPrioritizationAgent
from services.target_enrichment import TargetEnrichmentService
from services.disease_treatments import get_disease_treatments
from services.agentic_research import AgenticResearchService
from services.evidence_store import normalize_target_evidence, normalize_drug_evidence
from services.research_assistant import answer_research_question
from main import build_target_summary, build_evidence_summary




def empty_outputs():
   
    return (
        pd.DataFrame(),  # Ranked targets
        pd.DataFrame(),  # UniProt
        pd.DataFrame(),  # PubMed
        pd.DataFrame(),  # Reactome
        None,             # STRING network HTML
        pd.DataFrame(),  # GO enrichment
        pd.DataFrame(),  # KEGG pathways
        pd.DataFrame(),  # Clinical candidates
        pd.DataFrame(),  # Disease treatments
        {},               # Detailed JSON
        None,             # CSV download
        None,             # JSON download
    )


def render_string_network(network):
    """Render a confidence-filterable STRING association network as SVG."""
    network = network or {}
    nodes = network.get("nodes") or []
    edges = network.get("edges") or []

    if not nodes or not edges:
        return (
            "<div style='padding:24px;border:1px solid #17445a;"
            "border-radius:14px;color:#9bb3c5;background:#081b2b'>"
            "<strong>STRING network unavailable</strong>"
            "<p>No protein associations were returned for this target.</p>"
            "</div>"
        )

    channels = [
        ("nscore", "Neighborhood", "#2878b5"),
        ("fscore", "Gene fusion", "#8b5cf6"),
        ("pscore", "Co-occurrence", "#0f9d91"),
        ("ascore", "Co-expression", "#39945b"),
        ("escore", "Experiments", "#d97732"),
        ("dscore", "Databases", "#d34d5d"),
        ("tscore", "Text mining", "#77838f"),
    ]

    def number(value, default=0.0):
        try:
            return float(value) if value is not None else default
        except (TypeError, ValueError):
            return default

    valid_nodes = [
        node for node in nodes
        if node.get("string_id")
    ]
    node_ids = {node["string_id"] for node in valid_nodes}
    valid_edges = [
        edge for edge in edges
        if edge.get("protein_a_id") in node_ids
        and edge.get("protein_b_id") in node_ids
    ]

    if not valid_nodes or not valid_edges:
        return (
            "<div style='padding:20px;border:1px solid #17445a;"
            "border-radius:14px;color:#9bb3c5;background:#081b2b'>"
            "The network has no valid node-edge connections to display."
            "</div>"
        )

    degree = {node["string_id"]: 0 for node in valid_nodes}
    for edge in valid_edges:
        degree[edge["protein_a_id"]] += 1
        degree[edge["protein_b_id"]] += 1

    node_by_id = {node["string_id"]: node for node in valid_nodes}
    width, height = 900, 650
    cx, cy = width / 2, height / 2
    radius = min(width, height) * 0.35

    
    ordered_nodes = sorted(
        valid_nodes,
        key=lambda node: (-degree[node["string_id"]],
                          str(node.get("preferred_name") or ""))
    )
    positions = {}
    positions[ordered_nodes[0]["string_id"]] = (cx, cy)

    remaining = ordered_nodes[1:]
    for index, node in enumerate(remaining):
        angle = 2 * math.pi * index / max(1, len(remaining)) - math.pi / 2
        ring_radius = radius if index % 2 == 0 else radius * 0.78
        positions[node["string_id"]] = (
            cx + ring_radius * math.cos(angle),
            cy + ring_radius * math.sin(angle),
        )

    def channel_for(edge):
        evidence = edge.get("evidence_scores") or {}
        available = [
            (number(evidence.get(key)), key, label, color)
            for key, label, color in channels
            if evidence.get(key) is not None
        ]
        if not available:
            return None, "Unspecified evidence", "#a3adb5"
        _, key, label, color = max(available, key=lambda item: item[0])
        return key, label, color

   
    svg_id = "string-network-" + str(abs(hash(
        tuple(sorted(str(node["string_id"]) for node in valid_nodes))
    )))
    edge_markup = []

    for index, edge in enumerate(valid_edges):
        a = positions.get(edge.get("protein_a_id"))
        b = positions.get(edge.get("protein_b_id"))
        if not a or not b:
            continue

        score = number(edge.get("score"), 0)
        score_1000 = number(
            edge.get("score_0_1000"),
            score * 1000 if score <= 1 else score
        )
        confidence = max(0, min(1000, score_1000))
        _, channel_label, color = channel_for(edge)
        protein_a = str(edge.get("protein_a") or "?")
        protein_b = str(edge.get("protein_b") or "?")

        tooltip = escape(
            f"{protein_a} — {protein_b}; "
            f"STRING confidence: {confidence:.0f}/1000; "
            f"strongest supplied channel: {channel_label}"
        )
        edge_markup.append(
            f"<line class='network-edge' data-confidence='{confidence:.1f}' "
            f"x1='{a[0]:.1f}' y1='{a[1]:.1f}' "
            f"x2='{b[0]:.1f}' y2='{b[1]:.1f}' "
            f"stroke='{color}' stroke-opacity='.68' "
            f"stroke-width='{1 + confidence / 1000 * 2.6:.2f}'>"
            f"<title>{tooltip}</title></line>"
        )

    node_markup = []
    for node in ordered_nodes:
        node_id = node["string_id"]
        x, y = positions[node_id]
        name = str(node.get("preferred_name") or "Unknown")
        safe_name = escape(name)
        hub_degree = degree[node_id]
        node_radius = min(25, 10 + hub_degree * 1.3)
        is_hub = node_id == ordered_nodes[0]["string_id"]
        fill = "#167d78" if is_hub else "#244c61"
        tooltip = escape(
            f"{name}; connections in displayed network: {hub_degree}"
            + ("; highest-degree node" if is_hub else "")
        )
        node_markup.append(
            f"<g class='network-node' tabindex='0' role='img' "
            f"aria-label='{tooltip}'>"
            f"<title>{tooltip}</title>"
            f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{node_radius:.1f}' "
            f"fill='{fill}' stroke='#ffffff' stroke-width='2.5'/>"
            f"<text x='{x:.1f}' y='{y + node_radius + 16:.1f}' "
            f"text-anchor='middle' font-size='12' font-weight='600' "
            f"fill='#172b3a'>{safe_name}</text></g>"
        )

    legend_items = "".join(
        "<span style='display:inline-flex;align-items:center;gap:6px;"
        "font-size:12px;color:#465966;margin:4px 12px 4px 0'>"
        f"<span style='width:9px;height:9px;border-radius:50%;"
        f"background:{color};display:inline-block'></span>{escape(label)}</span>"
        for _, label, color in channels
    )
    max_score = max(
        number(edge.get("score_0_1000"),
               number(edge.get("score")) * 1000)
        for edge in valid_edges
    )
    edge_count = len(edge_markup)
    node_count = len(ordered_nodes)
    hub_name = str(ordered_nodes[0].get("preferred_name") or "Unknown")

    return (
        "<div style='background:#081b2b;border:1px solid #17445a;"
        "border-radius:16px;padding:16px;color:#edf7ff'>"
        "<div style='display:flex;flex-wrap:wrap;gap:18px;"
        "align-items:center;justify-content:space-between'>"
        "<div><div style='font-size:11px;letter-spacing:1.2px;"
        "color:#9bb3c5;font-weight:700'>PROTEIN ASSOCIATION MAP</div>"
        f"<div style='font-size:22px;font-weight:700;margin-top:4px'>"
        f"{node_count} proteins · {edge_count} associations</div></div>"
        "<div style='font-size:12px;color:#9bb3c5'>"
        f"Highest degree: <strong>{escape(hub_name)}</strong>"
        "</div></div>"
        "<div style='margin:16px 0 8px;display:flex;flex-wrap:wrap;"
        "align-items:center;gap:12px'>"
        "<label style='font-size:13px;font-weight:600' "
        f"for='{svg_id}-filter'>Minimum confidence: "
        "<span id='" + svg_id + "-value'>0</span>/1000</label>"
        f"<input id='{svg_id}-filter' type='range' min='0' max='1000' "
        "step='50' value='0' style='flex:1;min-width:160px' "
        f"oninput=\"(function(s){{var root=document.getElementById('{svg_id}');"
        f"document.getElementById('{svg_id}-value').textContent=s.value;"
        "root.querySelectorAll('.network-edge').forEach(function(e){"
        "e.style.display=Number(e.dataset.confidence)>=Number(s.value)"
        "?'':'none';});})(this)\">"
        "</div>"
        f"<svg id='{svg_id}' viewBox='0 0 {width} {height}' "
        "role='img' aria-label='STRING protein association network' "
        "style='width:100%;max-height:650px;background:#fafbf8;"
        "border:1px solid #edf0ed;border-radius:12px'>"
        + "".join(edge_markup)
        + "".join(node_markup)
        + "</svg>"
        "<div style='margin-top:12px'><div style='font-size:11px;"
        "font-weight:700;letter-spacing:.8px;color:#9bb3c5;"
        "margin-bottom:6px'>EVIDENCE CHANNELS</div>"
        + legend_items
        + "</div>"
        "<p style='font-size:12px;line-height:1.6;color:#9bb3c5;"
        "margin:12px 0 0'>Edge color indicates the strongest supplied "
        "STRING evidence channel; edges may have support from multiple "
        "channels. Edge width reflects the combined STRING confidence score. "
        "Node size reflects degree within this returned network. STRING "
        "associations do not necessarily represent direct physical binding "
        "or clinical efficacy.</p></div>"
    )

def analyze_disease(disease_name, top_n):
    if not disease_name or not disease_name.strip():
        return ("Enter a disease name to begin.", *empty_outputs())

    try:
        agent = TargetPrioritizationAgent()
        result = agent.run(disease_name.strip(), int(top_n))
    except (requests.RequestException, RuntimeError, ValueError) as exc:
        return (f"**Analysis failed:** {exc}", *empty_outputs())

    disease = result["disease"]
    targets = result["targets"]
    candidate = result["top_candidate"]

    targets_df = pd.DataFrame([
        {
            "Rank": rank,
            "Gene": target.symbol,
            "Name": target.name,
            "Ensembl ID": target.gene_id,
            "Association score": (
                round(float(target.association_score), 4)
                if target.association_score is not None else None
            ),
        }
        for rank, target in enumerate(targets, start=1)
    ])

    if candidate is None:
        summary = (
            f"## {disease.name}\n\n"
            f"Disease ID: `{disease.disease_id}`\n\n"
            "No target candidates were returned."
        )
        return (summary, *empty_outputs())

    service = TargetEnrichmentService()
    warnings = []

    try:
        enrichment = service.enrich_target(candidate, disease.name)
    except (requests.RequestException, RuntimeError, ValueError) as exc:
        warnings.append(f"Protein/publication/pathway enrichment failed: {exc}")
        enrichment = {
            "uniprot_entries": [],
            "publications": [],
            "pathways": [],
        }

    drug_status = {
        "status": "failed",
        "source": "Open Targets Platform",
        "record_count": None,
    }

    try:
        drugs = service.get_known_drugs(candidate.gene_id)
        drug_status = {
            "status": "retrieved" if drugs else "empty",
            "source": "Open Targets Platform",
            "record_count": len(drugs),
        }
    except (requests.RequestException, RuntimeError, ValueError) as exc:
        drugs = []
        drug_status["error"] = str(exc)
        warnings.append(f"Clinical candidate lookup failed: {exc}")

    treatments = get_disease_treatments(disease.name)

    target_summary = build_target_summary(candidate, targets)

    evidence_summary = build_evidence_summary(
        enrichment=enrichment,
        drugs=drugs,
        treatments=treatments,
        drug_status=drug_status,
    )

    proteins_df = pd.DataFrame(enrichment["uniprot_entries"])
    publications_df = pd.DataFrame(enrichment["publications"])
    pathways_df = pd.DataFrame(enrichment["pathways"])
    go_df = pd.DataFrame(enrichment.get("go_enrichment", []))
    kegg_df = pd.DataFrame(enrichment.get("kegg_enrichment", []))
    network_html = render_string_network(enrichment.get("string_network"))
    drugs_df = pd.json_normalize(drugs) if drugs else pd.DataFrame()
    treatments_df = pd.DataFrame(treatments)

    # Compact overview cards, using the existing analysis results.
    evidence_categories = evidence_summary.get("categories", {})
    retrieved_count = sum(
        1 for item in evidence_categories.values()
        if item.get("status") == "retrieved"
    )
    failed_count = sum(
        1 for item in evidence_categories.values()
        if item.get("status") == "failed"
    )
    evidence_total = len(evidence_categories)

    score_text = f"{candidate.association_score:.4f}"
    safe_disease = escape(str(disease.name))
    safe_symbol = escape(str(candidate.symbol))
    safe_target_name = escape(str(candidate.name or "Name unavailable"))
    safe_disease_id = escape(str(disease.disease_id or "Not available"))
    safe_gene_id = escape(str(candidate.gene_id or "Not available"))

    summary = (
        '<div style="margin:0 0 18px">'
        '<div style="font-size:11px;letter-spacing:1.4px;'
        'font-weight:700;color:#9bb3c5">D Y S C O V E R Y / ANALYSIS</div>'
        f'<h2 style="margin:6px 0 4px;font-size:28px;line-height:1.2;'
        f'color:#edf7ff">{safe_disease}</h2>'
        f'<div style="font-size:12px;color:#9bb3c5">Disease ID: '
        f'<code>{safe_disease_id}</code></div>'
        '</div>'
        '<div class="overview-card-grid" style="display:grid;'
        'gap:12px;margin:18px 0 24px">'
        '<div style="background:#081b2b;border:1px solid #17445a;'
        'border-radius:14px;padding:17px">'
        '<div style="font-size:11px;letter-spacing:.7px;color:#9bb3c5;'
        'font-weight:700">TOP PRIORITIZED TARGET</div>'
        f'<div style="font-size:25px;font-weight:750;color:#17e0c0;'
        f'margin-top:8px">{safe_symbol}</div>'
        f'<div style="font-size:13px;color:#edf7ff;margin-top:3px">'
        f'{safe_target_name}</div>'
        f'<div style="font-size:11px;color:#9bb3c5;margin-top:8px">'
        f'Ensembl: {safe_gene_id}</div>'
        '</div>'
        '<div style="background:#081b2b;border:1px solid #17445a;'
        'border-radius:14px;padding:17px">'
        '<div style="font-size:11px;letter-spacing:.7px;color:#9bb3c5;'
        'font-weight:700">ASSOCIATION SCORE</div>'
        f'<div style="font-size:29px;font-weight:750;color:#edf7ff;'
        f'margin-top:8px">{escape(score_text)}</div>'
        '<div style="font-size:11px;color:#9bb3c5;margin-top:7px">'
        'Ranking signal, not a probability of clinical efficacy.</div>'
        '</div>'
        '<div style="background:#081b2b;border:1px solid #17445a;'
        'border-radius:14px;padding:17px">'
        '<div style="font-size:11px;letter-spacing:.7px;color:#9bb3c5;'
        'font-weight:700">TARGETS DISPLAYED</div>'
        f'<div style="font-size:29px;font-weight:750;color:#edf7ff;'
        f'margin-top:8px">{len(targets)}</div>'
        '<div style="font-size:11px;color:#9bb3c5;margin-top:7px">'
        'Ranked candidates returned for this disease.</div>'
        '</div>'
        '<div style="background:#081b2b;border:1px solid #17445a;'
        'border-radius:14px;padding:17px">'
        '<div style="font-size:11px;letter-spacing:.7px;color:#9bb3c5;'
        'font-weight:700">EVIDENCE SOURCES RETRIEVED</div>'
        f'<div style="font-size:29px;font-weight:750;color:#17e0c0;'
        f'margin-top:8px">{retrieved_count}<span style="font-size:15px;'
        f'font-weight:500;color:#9bb3c5"> / {evidence_total}</span></div>'
        f'<div style="font-size:11px;color:#9bb3c5;margin-top:7px">'
        f'{failed_count} failed source(s); empty sources are counted '
        'separately.</div>'
        '</div>'
        '</div>'
        '<h3 style="margin:22px 0 8px;color:#edf7ff">'
        'Evidence retrieval status</h3>'
        '<div style="overflow-x:auto;margin:12px 0">'
        '<table style="width:100%;border-collapse:collapse;background:#081b2b">'
        '<thead><tr>'
        '<th style="text-align:left;padding:10px;border-bottom:1px solid #17445a">Evidence category</th>'
        '<th style="text-align:left;padding:10px;border-bottom:1px solid #17445a">Status</th>'
        '<th style="text-align:right;padding:10px;border-bottom:1px solid #17445a">Records</th>'
        '</tr></thead><tbody>'
    )

    category_labels = {
        "uniprot": "UniProt protein information",
        "publications": "PubMed publications",
        "reactome_pathways": "Reactome pathways",
        "string_network": "STRING protein network",
        "go_enrichment": "GO enrichment",
        "kegg_enrichment": "KEGG enrichment",
        "known_drugs": "Target-associated clinical candidates",
        "disease_treatments": "Curated disease treatments",
    }

    for key, label in category_labels.items():
        item = evidence_summary.get("categories", {}).get(key, {})
        status = item.get("status", "unavailable")
        count = item.get("record_count")
        count_text = str(count) if count is not None else "Not available"

        if key == "string_network" and item.get("node_count") is not None:
            count_text = (
                f"{item.get('node_count')} nodes, "
                f"{item.get('record_count', 0)} edges"
            )

        safe_label = escape(str(label))
        safe_status = escape(str(status))
        safe_count = escape(str(count_text))
        status_color = {
            "retrieved": "#167d78",
            "empty": "#647581",
            "failed": "#b54747",
        }.get(status, "#647581")
        summary += (
            "<tr>"
            f'<td style="padding:10px;border-bottom:1px solid #17445a">{safe_label}</td>'
            f'<td style="padding:10px;border-bottom:1px solid #17445a;color:{status_color};font-weight:600">{safe_status}</td>'
            f'<td style="padding:10px;text-align:right;border-bottom:1px solid #17445a">{safe_count}</td>'
            "</tr>"
        )

    summary += (
        "</tbody></table></div>"
        "<p style='margin-top:12px;color:#9bb3c5;font-size:13px'>"
        "<strong>Status meanings:</strong> <code>retrieved</code> means data was returned; "
        "`empty` means the source returned no records; "
        "`failed` means retrieval failed.\n\n"
        "Disease-level treatments are not necessarily validated against "
        "the selected target."
    )

    if warnings:
        summary += "\n\n### Warnings\n" + "\n".join(
            f"- {warning}" for warning in warnings
        )

    details = {
        "disease": disease.to_dict(),
        "selected_target": candidate.to_dict(),
        "ranked_targets": [target.to_dict() for target in targets],
        "uniprot": enrichment["uniprot_entries"],
        "publications": enrichment["publications"],
        "pathways": enrichment["pathways"],
        "string_network": enrichment.get("string_network", {"nodes": [], "edges": []}),
        "go_enrichment": enrichment.get("go_enrichment", []),
        "kegg_enrichment": enrichment.get("kegg_enrichment", []),
        "known_drugs": drugs,
        "disease_treatments": treatments,
        "target_summary": target_summary,
        "evidence_summary": evidence_summary,
        "warnings": warnings,
        "interpretation": (
            "Association scores rank evidence and are not probabilities "
            "of clinical efficacy. Curated disease-level treatments are "
            "separate from target-associated clinical candidates."
        ),
    }

    columns = [
        "rank", "disease", "disease_id", "gene_id", "symbol", "name",
        "association_score", "evidence_source", "selected_top_candidate",
    ]

    csv_rows = []
    for rank, target in enumerate(targets, start=1):
        csv_rows.append({
            "rank": rank,
            "disease": disease.name,
            "disease_id": disease.disease_id,
            **target.to_dict(),
            "selected_top_candidate": target.gene_id == candidate.gene_id,
        })

    csv_df = pd.DataFrame(csv_rows, columns=columns)

    os.makedirs("results", exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    csv_path = os.path.abspath(f"results/web_report_{stamp}.csv")
    json_path = os.path.abspath(f"results/web_report_{stamp}_details.json")

    csv_df.to_csv(csv_path, index=False)
    with open(json_path, "w", encoding="utf-8") as file:
        json.dump(details, file, indent=2, ensure_ascii=False)

    return (
        summary,
        targets_df,
        proteins_df,
        publications_df,
        pathways_df,
        network_html,
        go_df,
        kegg_df,
        drugs_df,
        treatments_df,
        details,
        csv_path,
        json_path,
    )




CSS = """
:root {
  color-scheme: dark;
  --bg: #03111f;
  --surface: #081b2b;
  --surface-2: #0c2436;
  --ink: #edf7ff;
  --muted: #9bb3c5;
  --cyan: #19c9e5;
  --teal: #17e0c0;
  --line: #17445a;
  --paper: #071a29;
}

html, body, .main, .wrap {
  background: var(--bg) !important;
  color: var(--ink) !important;
}

.gradio-container {
  max-width: 1500px !important;
  margin: auto !important;
  padding: 0 clamp(14px, 3vw, 42px) 38px !important;
  background:
    radial-gradient(ellipse at 12% 8%, rgba(0, 143, 190, .12), transparent 32%),
    radial-gradient(ellipse at 92% 18%, rgba(0, 90, 160, .12), transparent 30%),
    var(--bg) !important;
  color: var(--ink) !important;
  font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif !important;
}

.gradio-container p,
.gradio-container label,
.gradio-container .prose,
.gradio-container .markdown,
.gradio-container .block {
  color: var(--ink) !important;
}

#brand-header {
  position: relative;
  overflow: hidden;
  padding: clamp(25px, 4vw, 48px) clamp(2px, 1vw, 12px) 30px;
  margin: 0 0 30px;
  border-bottom: 1px solid var(--line);
  background:
    radial-gradient(ellipse at 86% 0%, rgba(0, 174, 220, .12), transparent 35%),
    linear-gradient(115deg, rgba(3, 17, 31, .15), rgba(8, 38, 57, .25));
}

#brand-header:after {
  content: "";
  position: absolute;
  pointer-events: none;
  right: -5%;
  bottom: -55px;
  width: 58%;
  height: 100px;
  opacity: .48;
  background: repeating-radial-gradient(
    ellipse at 50% 100%,
    transparent 0 13px,
    rgba(0, 177, 226, .22) 14px 15px,
    transparent 16px 24px
  );
}

.brand-top {
  position: relative;
  z-index: 1;
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 28px;
}

.brand-identity {
  display: flex;
  align-items: center;
  gap: 16px;
  min-width: 0;
}

.brand-symbol {
  width: clamp(54px, 6vw, 76px);
  height: clamp(54px, 6vw, 76px);
  flex: 0 0 auto;
  filter: drop-shadow(0 0 14px rgba(25, 201, 229, .22));
}

.brand-title {
  color: #f4f9ff !important;
  font-size: clamp(35px, 5.4vw, 66px);
  font-weight: 780;
  letter-spacing: -2.8px;
  line-height: 1;
  margin: 0;
}

.brand-title .brand-o {
  color: var(--cyan);
  text-shadow: 0 0 18px rgba(25, 201, 229, .35);
}

.brand-kicker {
  color: #d5e7f2 !important;
  font-size: clamp(10px, 1.25vw, 15px);
  letter-spacing: clamp(1.5px, .28vw, 3.4px);
  font-weight: 500;
  margin-top: 11px;
}

.brand-tagline {
  flex: 0 0 auto;
  border-left: 2px solid var(--cyan);
  padding: 8px 0 8px 24px;
  color: #e5f5ff !important;
  font-size: clamp(12px, 1.5vw, 19px);
  letter-spacing: clamp(1px, .22vw, 2.5px);
  line-height: 1.7;
  white-space: nowrap;
}

.brand-tagline span {
  color: var(--teal) !important;
  padding: 0 7px;
}

.section-heading {
  color: var(--cyan) !important;
  font-size: 11px;
  letter-spacing: 1.8px;
  font-weight: 750;
  margin: 18px 0 10px;
}

.gradio-container .panel,
.gradio-container .form,
.gradio-container fieldset,
.gradio-container .tabitem,
.gradio-container .block,
.gradio-container .wrap {
  border-color: var(--line) !important;
}

.gradio-container input,
.gradio-container textarea,
.gradio-container select {
  background: var(--surface) !important;
  color: var(--ink) !important;
  border: 1px solid #20516a !important;
  border-radius: 9px !important;
}

.gradio-container input:focus,
.gradio-container textarea:focus {
  border-color: var(--cyan) !important;
  box-shadow: 0 0 0 2px rgba(25, 201, 229, .15) !important;
}

.gradio-container input::placeholder,
.gradio-container textarea::placeholder {
  color: #819db1 !important;
  opacity: 1 !important;
}

.gradio-container label,
.gradio-container .label-wrap {
  color: #d8eaf6 !important;
  font-weight: 600 !important;
}

.gradio-container button.primary {
  background: linear-gradient(100deg, #08aeca, #13cbbd) !important;
  color: #031522 !important;
  border: 1px solid #19c9e5 !important;
  border-radius: 9px !important;
  font-weight: 750 !important;
  box-shadow: 0 0 20px rgba(0, 190, 220, .12) !important;
}

.gradio-container button.primary:hover {
  filter: brightness(1.12);
}

.gradio-container button:not(.primary) {
  border-color: var(--line) !important;
}

.gradio-container .tab-nav {
  border-bottom: 1px solid var(--line) !important;
  gap: 5px;
}

.gradio-container .tab-nav button {
  color: var(--muted) !important;
  font-weight: 600 !important;
  background: transparent !important;
}

.gradio-container .tab-nav button.selected {
  color: var(--cyan) !important;
  border-color: var(--cyan) !important;
}

.gradio-container table,
.gradio-container .table-wrap,
.gradio-container .dataframe {
  background: var(--surface) !important;
  color: var(--ink) !important;
  border-color: var(--line) !important;
}

.gradio-container th {
  background: #102c40 !important;
  color: #e5f6ff !important;
  font-weight: 700 !important;
}

.gradio-container td {
  color: #d4e6f2 !important;
}

.gradio-container .prose h1,
.gradio-container .prose h2,
.gradio-container .prose h3,
.gradio-container .prose p,
.gradio-container .prose li {
  color: var(--ink) !important;
}

#developer-card {
  background: var(--surface);
  border: 1px solid var(--line);
  border-radius: 14px;
  padding: clamp(20px, 3vw, 30px);
  margin-top: 18px;
  box-shadow: 0 12px 36px rgba(0, 0, 0, .12);
}

.developer-name {
  color: #f0f8ff !important;
  font-size: 23px;
  font-weight: 750;
  letter-spacing: -.5px;
}

.developer-title {
  color: var(--cyan) !important;
  font-weight: 650;
  margin-top: 5px;
}

.developer-link {
  display: inline-block;
  color: #031522 !important;
  background: var(--cyan);
  border: 1px solid var(--cyan);
  border-radius: 8px;
  padding: 10px 13px;
  text-decoration: none !important;
  overflow-wrap: anywhere;
  margin-top: 10px;
  font-weight: 700;
}

.disclaimer {
  color: var(--muted) !important;
  font-size: 12px;
  line-height: 1.8;
  border-top: 1px solid var(--line);
  padding-top: 18px;
  margin-top: 28px;
}

footer { display: none !important; }

.overview-card-grid {
  grid-template-columns: repeat(4, minmax(0, 1fr));
}

@media (max-width: 850px) {
  .brand-top {
    align-items: flex-start;
    flex-direction: column;
    gap: 18px;
  }
  .brand-tagline {
    border-left: 0;
    border-top: 1px solid var(--line);
    padding: 12px 0 0;
    white-space: normal;
  }
  .overview-card-grid {
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }
}

@media (max-width: 600px) {
  .brand-identity { gap: 11px; }
  .brand-title {
    font-size: clamp(31px, 9vw, 42px);
    letter-spacing: -1.8px;
  }
  .brand-symbol {
    width: 48px;
    height: 48px;
  }
  .brand-kicker {
    font-size: 9px;
    letter-spacing: 1.2px;
  }
  #brand-header { margin-bottom: 24px; }
  .overview-card-grid {
    grid-template-columns: minmax(0, 1fr);
  }
}

/* DYScovery UI refinement overrides */
.gradio-container,
.gradio-container .contain,
.gradio-container .form,
.gradio-container .tabitem,
.gradio-container .panel,
.gradio-container [data-testid="block"],
.gradio-container [data-testid="tabitem"] {
  background-color: var(--bg) !important;
  color: var(--ink) !important;
}

.gradio-container .panel,
.gradio-container [data-testid="form"],
.gradio-container fieldset {
  background: linear-gradient(145deg, #071c2c, #061522) !important;
  border: 1px solid var(--line) !important;
  border-radius: 14px !important;
}

.gradio-container .tab-nav {
  background: transparent !important;
  border-bottom: 1px solid var(--line) !important;
  overflow-x: auto !important;
}

.gradio-container .tab-nav button {
  color: #a8bfd0 !important;
  border-bottom: 2px solid transparent !important;
  white-space: nowrap !important;
}

.gradio-container .tab-nav button:hover {
  color: #eafaff !important;
  background: rgba(25, 201, 229, .06) !important;
}

.gradio-container .tab-nav button.selected {
  color: var(--cyan) !important;
  border-bottom-color: var(--cyan) !important;
  background: rgba(25, 201, 229, .07) !important;
}

.gradio-container table,
.gradio-container .table-wrap,
.gradio-container .dataframe,
.gradio-container [data-testid="dataframe"] {
  background: #071a29 !important;
  color: var(--ink) !important;
  border: 1px solid var(--line) !important;
  border-radius: 10px !important;
  overflow: hidden !important;
}

.gradio-container th,
.gradio-container thead,
.gradio-container thead tr {
  background: #103047 !important;
  color: #e9f8ff !important;
}

.gradio-container td {
  background: #081c2b !important;
  color: #d8e9f5 !important;
}

.gradio-container tbody tr:nth-child(even) td {
  background: #0b2233 !important;
}

.gradio-container tbody tr:hover td {
  background: #103047 !important;
}

.gradio-container input,
.gradio-container textarea,
.gradio-container select,
.gradio-container [role="textbox"] {
  background: #081c2b !important;
  color: #f0f8ff !important;
  border-color: #20516a !important;
}

.gradio-container input::placeholder,
.gradio-container textarea::placeholder {
  color: #8ba8bb !important;
}

.gradio-container button.primary {
  background: linear-gradient(100deg, #10d5d8, #08a9f4) !important;
  color: #031522 !important;
  border: 0 !important;
  border-radius: 10px !important;
  box-shadow: 0 5px 20px rgba(0, 190, 220, .16) !important;
}

.gradio-container button.primary:hover {
  filter: brightness(1.1);
}

.gradio-container pre,
.gradio-container code,
.gradio-container [data-testid="json"],
.gradio-container .json-holder {
  background: #061522 !important;
  color: #c8f5ff !important;
  border-color: var(--line) !important;
  border-radius: 8px !important;
}

.gradio-container .prose,
.gradio-container .markdown,
.gradio-container .prose p,
.gradio-container .prose li,
.gradio-container .prose h1,
.gradio-container .prose h2,
.gradio-container .prose h3,
.gradio-container .prose strong {
  color: #e2eef8 !important;
}

.gradio-container .prose a {
  color: #39d7ef !important;
}

.gradio-container hr {
  border-color: var(--line) !important;
}

.gradio-container input[type="range"] {
  accent-color: var(--cyan) !important;
}

@media (max-width: 700px) {
  .gradio-container {
    padding-left: 12px !important;
    padding-right: 12px !important;
  }
}


/* DYScovery inline-summary dark theme */
.gradio-container [style*="background:#081b2b"],
.gradio-container [style*="background: #081b2b"],
.gradio-container [style*="background:#081b2bfff"],
.gradio-container [style*="background: #081b2bfff"] {
  background: #081b2b !important;
  color: #edf7ff !important;
  border-color: #17445a !important;
}

.gradio-container [style*="color:#edf7ff"],
.gradio-container [style*="color: #edf7ff"],
.gradio-container [style*="color:#294541"],
.gradio-container [style*="color: #294541"] {
  color: #edf7ff !important;
}

.gradio-container [style*="color:#9bb3c5"],
.gradio-container [style*="color: #9bb3c5"] {
  color: #9bb3c5 !important;
}

.gradio-container [style*="border-bottom:1px solid #17445a"],
.gradio-container [style*="border-bottom:1px solid #17445a"] {
  border-bottom-color: #17445a !important;
}

.gradio-container [style*="border:1px solid #17445a"],
.gradio-container [style*="border: 1px solid #17445a"] {
  border-color: #17445a !important;
}

.gradio-container [data-testid="textbox"],
.gradio-container [data-testid="slider"],
.gradio-container [data-testid="number"],
.gradio-container [data-testid="form"],
.gradio-container .panel {
  background-color: #081b2b !important;
  color: #edf7ff !important;
  border-color: #17445a !important;
}

.gradio-container [data-testid="textbox"] input,
.gradio-container [data-testid="number"] input {
  background-color: #071a29 !important;
  color: #edf7ff !important;
}

.gradio-container [data-testid="slider"] input {
  accent-color: #19c9e5 !important;
}

.gradio-container [data-testid="form"] > *,
.gradio-container .panel > * {
  color: #edf7ff !important;
}


/* Integrated Dyscovery hero and search layout */
#brand-header {
  padding: 34px 18px 24px !important;
  margin-bottom: 24px !important;
  border-bottom: 1px solid #17445a !important;
  background:
    radial-gradient(ellipse at 0% 15%, rgba(0, 181, 220, .16), transparent 35%),
    radial-gradient(ellipse at 100% 5%, rgba(0, 106, 200, .15), transparent 38%),
    linear-gradient(145deg, #03111f, #061a2a) !important;
}

.brand-hero {
  text-align: center;
  position: relative;
  z-index: 1;
}

.brand-wordmark {
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 0;
  color: #f4f9ff;
  font-size: clamp(46px, 7vw, 88px);
  line-height: 1;
  font-weight: 780;
  letter-spacing: -4px;
  text-shadow: 0 0 24px rgba(25, 201, 229, .08);
}

.wordmark-d {
  width: clamp(45px, 6vw, 72px);
  height: clamp(50px, 6.4vw, 78px);
  margin-right: 2px;
  flex: 0 0 auto;
}

.wordmark-o {
  width: .77em;
  height: .77em;
  margin: 0 .015em;
  position: relative;
  display: inline-block;
  flex: 0 0 auto;
  border: 3px solid #19c9e5;
  border-radius: 50%;
  box-shadow: 0 0 14px rgba(25, 201, 229, .2);
}

.wordmark-o:before,
.wordmark-o:after {
  content: "";
  position: absolute;
  background: #19c9e5;
  left: 50%;
  top: -14%;
  transform: translateX(-50%);
  width: 2px;
  height: 128%;
}

.wordmark-o:after {
  left: -14%;
  top: 50%;
  transform: translateY(-50%);
  width: 128%;
  height: 2px;
}

.wordmark-crosshair {
  position: absolute;
  inset: 29%;
  border: 2px solid #19c9e5;
  border-radius: 50%;
  background: #03111f;
}

.brand-kicker {
  margin-top: 12px !important;
  color: #e3f0fa !important;
  font-size: clamp(13px, 1.7vw, 23px) !important;
  font-weight: 450 !important;
  letter-spacing: clamp(1.5px, .35vw, 4px) !important;
  text-transform: none !important;
}

.brand-divider {
  height: 2px;
  width: 100px;
  margin: 14px auto 12px;
  background: linear-gradient(90deg, #19e0d0, #347eff);
}

.brand-tagline {
  border: 0 !important;
  padding: 0 !important;
  white-space: normal !important;
  color: #e8f5ff !important;
  font-size: clamp(13px, 1.65vw, 21px) !important;
  letter-spacing: clamp(1px, .22vw, 2.5px) !important;
}

.brand-tagline span {
  color: #19c9e5 !important;
  padding: 0 8px !important;
}

.brand-features {
  display: flex;
  align-items: flex-start;
  justify-content: center;
  gap: clamp(15px, 4vw, 58px);
  flex-wrap: wrap;
  margin: 26px auto 0;
}

.brand-features > div {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 7px;
  color: #d8eaf6;
  font-size: 10px;
  letter-spacing: 1.2px;
  line-height: 1.45;
  text-transform: uppercase;
}

.feature-icon {
  color: #19c9e5;
  font-size: 30px;
  line-height: 1;
  font-weight: 400;
}

/* DYScovery high-contrast dark UI fix */
:root {
  color-scheme: dark !important;
  --body-background-fill: #03111f !important;
  --background-fill-primary: #081b2b !important;
  --background-fill-secondary: #0c2436 !important;
  --block-background-fill: #081b2b !important;
  --input-background-fill: #071a29 !important;
  --border-color-primary: #17445a !important;
  --body-text-color: #edf7ff !important;
  --body-text-color-subdued: #a9c1d2 !important;
  --color-accent: #19c9e5 !important;
}

html, body, body > div, .main, .wrap,
.gradio-container, .gradio-container .main {
  background-color: #03111f !important;
  color: #edf7ff !important;
}

.gradio-container,
.gradio-container * {
  scrollbar-color: #176078 #071a29;
}

.gradio-container .panel,
.gradio-container .search-panel,
.gradio-container [data-testid="form"],
.gradio-container [data-testid="block"],
.gradio-container .block,
.gradio-container fieldset,
.gradio-container .tabitem,
.gradio-container [data-testid="tabitem"] {
  background-color: #081b2b !important;
  color: #edf7ff !important;
  border-color: #17445a !important;
}

.gradio-container [data-testid="textbox"],
.gradio-container [data-testid="number"],
.gradio-container [data-testid="slider"],
.gradio-container [data-testid="dropdown"],
.gradio-container [data-testid="input"],
.gradio-container .wrap,
.gradio-container .input-container {
  background-color: #081b2b !important;
  color: #edf7ff !important;
  border-color: #20516a !important;
}

.gradio-container input,
.gradio-container textarea,
.gradio-container select,
.gradio-container input[type="number"] {
  background-color: #071a29 !important;
  color: #f2f9ff !important;
  border: 1px solid #20516a !important;
  border-radius: 8px !important;
  -webkit-text-fill-color: #f2f9ff !important;
}

.gradio-container input::placeholder,
.gradio-container textarea::placeholder {
  color: #8eacbf !important;
  -webkit-text-fill-color: #8eacbf !important;
  opacity: 1 !important;
}

.gradio-container label,
.gradio-container .label-wrap,
.gradio-container [data-testid="block"] label {
  color: #c8dce9 !important;
}

.gradio-container button.primary {
  background: linear-gradient(105deg, #11d5d3, #079ff0) !important;
  color: #031522 !important;
  -webkit-text-fill-color: #031522 !important;
  border: 1px solid #19c9e5 !important;
  border-radius: 10px !important;
  font-weight: 750 !important;
}

.gradio-container .tab-nav {
  background: #03111f !important;
  border-bottom: 1px solid #17445a !important;
}

.gradio-container .tab-nav button {
  background: transparent !important;
  color: #b3c8d7 !important;
  border-color: transparent !important;
}

.gradio-container .tab-nav button.selected {
  color: #27d9ef !important;
  border-bottom: 2px solid #19c9e5 !important;
  background: #092236 !important;
}

.gradio-container table,
.gradio-container .table-wrap,
.gradio-container .dataframe,
.gradio-container [data-testid="dataframe"] {
  background: #071a29 !important;
  color: #edf7ff !important;
  border-color: #17445a !important;
}

.gradio-container th,
.gradio-container thead,
.gradio-container thead tr {
  background: #103047 !important;
  color: #eaf8ff !important;
}

.gradio-container td {
  background: #081b2b !important;
  color: #d8e9f5 !important;
}

.gradio-container tbody tr:nth-child(even) td {
  background: #0b2233 !important;
}

.gradio-container pre,
.gradio-container code,
.gradio-container [data-testid="json"],
.gradio-container .json-holder {
  background: #061522 !important;
  color: #c8f5ff !important;
  border-color: #17445a !important;
}

.gradio-container .prose,
.gradio-container .markdown,
.gradio-container .prose p,
.gradio-container .prose li,
.gradio-container .prose h1,
.gradio-container .prose h2,
.gradio-container .prose h3,
.gradio-container .prose strong {
  color: #e2eef8 !important;
}

/* Fix inline-styled HTML cards and evidence summaries */
.gradio-container [style*="background:#fff"],
.gradio-container [style*="background: #fff"],
.gradio-container [style*="background:#ffffff"],
.gradio-container [style*="background: #ffffff"],
.gradio-container [style*="background: white"] {
  background: #081b2b !important;
  color: #edf7ff !important;
  border-color: #17445a !important;
}

.gradio-container [style*="color:#172b3a"],
.gradio-container [style*="color: #172b3a"],
.gradio-container [style*="color:#647581"],
.gradio-container [style*="color: #647581"] {
  color: #b4cad9 !important;
}

.brand-wordmark,
.brand-wordmark > span {
  color: #f4f9ff !important;
  text-shadow: 0 0 20px rgba(25, 201, 229, .15) !important;
}

.brand-wordmark .wordmark-o {
  color: #19c9e5 !important;
}

.brand-kicker,
.brand-tagline,
.brand-features > div {
  color: #d9eaf6 !important;
  opacity: 1 !important;
}

.brand-features .feature-icon {
  color: #19c9e5 !important;
}

#brand-header {
  min-height: 0 !important;
  padding-bottom: 24px !important;
  margin-bottom: 24px !important;
}

@media (max-width: 760px) {
  .gradio-container {
    padding-left: 12px !important;
    padding-right: 12px !important;
  }
}




/* Dyscovery: unified responsive search panel */
.gradio-container .search-panel {
  box-sizing: border-box !important;
  width: 100% !important;
  min-width: 0 !important;
  margin: 14px 0 26px !important;
  padding: clamp(14px, 2vw, 22px) !important;
  background: linear-gradient(135deg, #071b2a, #041522) !important;
  border: 1px solid #176078 !important;
  border-radius: 16px !important;
  box-shadow: 0 12px 32px rgba(0, 0, 0, .16) !important;
}

.gradio-container .search-panel > div,
.gradio-container .search-panel .form,
.gradio-container .search-panel [data-testid="block"] {
  box-sizing: border-box !important;
  min-width: 0 !important;
  max-width: 100% !important;
}

.gradio-container .search-controls {
  display: flex !important;
  flex-flow: row wrap !important;
  align-items: end !important;
  gap: 16px !important;
  width: 100% !important;
  min-width: 0 !important;
}

.gradio-container .search-controls > * {
  min-width: 0 !important;
  max-width: 100% !important;
}

.gradio-container .search-panel label,
.gradio-container .search-panel .label-wrap {
  color: #d9edf7 !important;
  font-size: 13px !important;
  font-weight: 650 !important;
  line-height: 1.4 !important;
  overflow-wrap: anywhere !important;
}

.gradio-container .search-panel input:not([type="range"]) {
  box-sizing: border-box !important;
  width: 100% !important;
  min-width: 0 !important;
  min-height: 48px !important;
  padding: 11px 14px !important;
  background: #061725 !important;
  border: 1px solid #20516a !important;
  border-radius: 9px !important;
  color: #f2faff !important;
  font-size: 15px !important;
}

.gradio-container .search-panel input:focus {
  border-color: #19c9e5 !important;
  box-shadow: 0 0 0 3px rgba(25, 201, 229, .10) !important;
}

.gradio-container .search-panel input[type="range"] {
  width: 100% !important;
  min-width: 0 !important;
  accent-color: #19c9e5 !important;
}

.gradio-container .search-panel button.primary {
  box-sizing: border-box !important;
  min-width: 0 !important;
  min-height: 48px !important;
  padding: 12px 16px !important;
  white-space: normal !important;
  background: linear-gradient(105deg, #13d6d1, #08aaf4) !important;
  border: 1px solid rgba(64, 232, 244, .75) !important;
  border-radius: 10px !important;
  color: #031522 !important;
  font-size: 15px !important;
  font-weight: 800 !important;
  box-shadow: 0 5px 18px rgba(0, 183, 222, .15) !important;
}

@media (max-width: 760px) {
  .gradio-container .search-controls {
    align-items: stretch !important;
    gap: 14px !important;
  }

  .gradio-container .search-controls > * {
    flex: 1 1 100% !important;
    width: 100% !important;
  }

  .gradio-container .search-panel {
    padding: 14px !important;
    border-radius: 12px !important;
  }

  .gradio-container .search-panel button.primary {
    width: 100% !important;
  }
}


/* Search panel: remove nested ring and uneven control backgrounds */
.gradio-container .search-panel {
  border: 0 !important;
  box-shadow: none !important;
  background: linear-gradient(135deg, #071b2a, #041522) !important;
}

.gradio-container .search-panel > div,
.gradio-container .search-controls,
.gradio-container .search-controls > *,
.gradio-container .search-controls [data-testid="block"] {
  background-color: transparent !important;
  box-shadow: none !important;
}

/* Keep slider label and numeric value from colliding */
.gradio-container .search-panel [data-testid="slider"] {
  overflow: visible !important;
  min-width: 0 !important;
}

.gradio-container .search-panel [data-testid="slider"] input[type="number"] {
  position: static !important;
  display: inline-block !important;
  box-sizing: border-box !important;
  width: 76px !important;
  min-width: 76px !important;
  height: 38px !important;
  margin: 0 0 8px auto !important;
  padding: 4px 8px !important;
  background: #061725 !important;
  color: #f2faff !important;
  border: 1px solid #20516a !important;
  border-radius: 8px !important;
  text-align: center !important;
}

.gradio-container .search-panel [data-testid="slider"] button {
  max-width: 38px !important;
}

.gradio-container .search-panel button.primary {
  align-self: end !important;
  margin-top: 0 !important;
}

/* Retain a single subtle outline around the search area */
.gradio-container .search-panel {
  outline: 1px solid #176078 !important;
  outline-offset: 0 !important;
}

@media (max-width: 760px) {
  .gradio-container .search-panel [data-testid="slider"] input[type="number"] {
    margin-left: auto !important;
  }
}

"""

HEADER = """
<header id="brand-header">
  <div class="brand-hero">
    <div class="brand-wordmark" aria-label="Dyscovery">
      <svg class="wordmark-d" viewBox="0 0 70 76" role="img" aria-label="DNA D">
        <defs>
          <linearGradient id="wordmark-d-gradient" x1="0" y1="0" x2="1" y2="1">
            <stop offset="0%" stop-color="#19e0d0"/>
            <stop offset="60%" stop-color="#19c9e5"/>
            <stop offset="100%" stop-color="#347eff"/>
          </linearGradient>
        </defs>
        <path d="M8 5 H28 C52 5 65 18 65 38 C65 58 52 71 28 71 H8 Z"
              fill="none" stroke="url(#wordmark-d-gradient)" stroke-width="5"/>
        <path d="M12 13 C40 20 12 30 40 38 C12 46 40 56 12 63"
              fill="none" stroke="url(#wordmark-d-gradient)" stroke-width="3.5" stroke-linecap="round"/>
        <path d="M15 18 H29 M20 27 H34 M20 49 H34 M15 58 H29"
              stroke="#eafcff" stroke-width="2.5" stroke-linecap="round"/>
      </svg>
      <span>ysc</span>
      <span class="wordmark-o" aria-label="o">
        <span class="wordmark-crosshair"></span>
      </span>
      <span>very</span>
    </div>
    <div class="brand-kicker">Disease-to-Target Research Lab</div>
    <div class="brand-divider"></div>
    <div class="brand-tagline">Bioinformatics <span>×</span> AI <span>×</span> Drug Discovery</div>
    <div class="brand-features">
      <div><span class="feature-icon">⌕</span><span>Disease<br/>Exploration</span></div>
      <div><span class="feature-icon">◎</span><span>Target<br/>Prioritisation</span></div>
      <div><span class="feature-icon">◇</span><span>Protein<br/>Analysis</span></div>
      <div><span class="feature-icon">▤</span><span>Literature<br/>Insights</span></div>
      <div><span class="feature-icon">⌘</span><span>Pathway<br/>Analysis</span></div>
      <div><span class="feature-icon">⬡</span><span>Drug Discovery<br/>Opportunities</span></div>
    </div>
  </div>
</header>
"""

ABOUT = """
<div id="developer-card">
  <div class="section-heading">ABOUT THE DEVELOPER</div>
  <div style="display:flex;align-items:center;gap:16px;margin-top:14px">
    <div style="width:58px;height:58px;flex:0 0 58px;border-radius:15px;
                background:#103248;display:grid;place-items:center;
                border:1px solid #20536a">
      <svg viewBox="0 0 48 48" width="42" height="42" aria-hidden="true">
        <path d="M11 4 C37 15 11 33 37 44 M37 4 C11 15 37 33 11 44"
              fill="none" stroke="#39d8e6" stroke-width="2.5"/>
        <g stroke="#35c5a0" stroke-width="2">
          <path d="M16 8 L32 8 M19 17 L29 17 M16 26 L32 26
                   M19 35 L29 35 M16 42 L32 42"/>
        </g>
      </svg>
    </div>
    <div>
      <div class="developer-name">Basit Mushtaq Wani</div>
      <div class="developer-title">Bioinformatician</div>
    </div>
  </div>
  <div style="color:#9cb6cb;margin-top:20px;font-size:13px">
    RESEARCH INTERESTS
  </div>
  <div class="brand-tags">
    <span class="brand-tag">Bioinformatics</span>
    <span class="brand-tag">AI/ML in Drug Discovery</span>
    <span class="brand-tag">Target Prioritization</span>
  </div>
  <div style="color:#9cb6cb;margin-top:20px;font-size:13px">
    PROJECT REPOSITORY
  </div>
  <a class="developer-link"
     href="https://github.com/wani34basit1-ctrl/DiseaseToTargetAgent"
     target="_blank" rel="noopener noreferrer">
     GitHub ↗ &nbsp; View Disease2Target Agent source code
  </a>
</div>
"""




def run_agentic_research(question):
    """Return the research answer and source-linked evidence."""
    if not isinstance(question, str) or not question.strip():
        return "Enter a disease or research question.", pd.DataFrame()

    try:
        result = AgenticResearchService().run(question.strip())
    except Exception as exc:
        return (
            f"AI research failed ({type(exc).__name__}). "
            "Check the application terminal for details.",
            pd.DataFrame(),
        )

    answer = result.get("research_answer") or {}
    answer_text = answer.get("answer", "No answer was generated.")
    if not isinstance(answer_text, str):
        answer_text = str(answer_text)

    if result.get("citation_validation_scope") == "label_integrity_only":
        answer_text += (
            "\n\nCitation limitation: labels are checked against available "
            "records, but scientific claim support is not verified."
        )

    rows = []
    for record in result.get("evidence") or []:
        if isinstance(record, dict):
            rows.append({
                "Evidence ID": record.get("evidence_id"),
                "Source": record.get("source"),
                "Evidence Type": record.get("evidence_type"),
                "Gene": record.get("gene_symbol"),
                "Disease": record.get("disease_name"),
                "Title": record.get("title"),
                "Source ID": record.get("source_id"),
                "URL": record.get("url"),
            })

    columns = [
        "Evidence ID", "Source", "Evidence Type", "Gene",
        "Disease", "Title", "Source ID", "URL",
    ]
    return answer_text, pd.DataFrame(rows, columns=columns)


def build_synthesis_evidence(details):
    """Normalize the evidence already collected by the standard pipeline."""
    from types import SimpleNamespace

    if not isinstance(details, dict):
        return []

    disease = details.get("disease") or {}
    disease_name = (
        disease.get("name")
        if isinstance(disease, dict)
        else str(disease)
    )

    records = []
    targets = details.get("ranked_targets") or []
    selected = details.get("selected_target") or {}

    enrichment = {
        "uniprot_entries": details.get("uniprot") or [],
        "publications": details.get("publications") or [],
        "pathways": details.get("pathways") or [],
        "string_network": details.get("string_network") or {},
        "go_enrichment": details.get("go_enrichment") or [],
        "kegg_enrichment": details.get("kegg_enrichment") or [],
    }

    # Normalize each ranked target with the evidence returned for the
    # selected target. Do not claim that selected-target enrichment
    # belongs to every other ranked target.
    for index, item in enumerate(targets):
        if not isinstance(item, dict):
            continue

        target = SimpleNamespace(
            symbol=item.get("symbol"),
            name=item.get("name"),
            id=item.get("gene_id") or item.get("id"),
            association_score=item.get("association_score"),
        )
        target_enrichment = enrichment if (
            item.get("gene_id") == selected.get("gene_id")
            or item.get("id") == selected.get("gene_id")
            or item.get("symbol") == selected.get("symbol")
        ) else {}

        records.extend(
            normalize_target_evidence(
                target, target_enrichment, disease_name=disease_name
            )
        )

   
    records.extend(
        normalize_drug_evidence(
            details.get("known_drugs") or [],
            disease_name=disease_name,
            gene_symbol=selected.get("symbol"),
        )
    )

   
    from services.evidence_store import _record

    for item in details.get("disease_treatments") or []:
        if not isinstance(item, dict):
            continue
        record = _record(
            source="Disease treatments",
            evidence_type="disease_treatment",
            title=(
                item.get("name")
                or item.get("drug")
                or item.get("title")
                or "Disease treatment"
            ),
            content=item,
            disease_name=disease_name,
            gene_symbol=None,
            source_id=item.get("id") or item.get("drug_id"),
            url=item.get("url"),
            metadata={
                "evidence_level": item.get("evidence_level", "disease-level"),
                "evidence_source": item.get(
                    "evidence_source", "curated_disease_treatment_catalogue"
                ),
            },
        )
        if record:
            records.append(record)

   
    summary = details.get("evidence_summary") or {}
    record = _record(
        source="Pipeline",
        evidence_type="retrieval_status",
        title="Evidence retrieval status and limitations",
        content={
            "categories": summary.get("categories", {}),
            "warnings": details.get("warnings", []),
            "interpretation": details.get("interpretation", ""),
        },
        disease_name=disease_name,
    )
    if record:
        records.append(record)

    return records


def run_integrated_search(query, top_n):
    """Run disease analysis once, then synthesize its existing evidence."""
    standard_results = analyze_disease(query, top_n)
    details = standard_results[10] if len(standard_results) > 10 else {}
    records = build_synthesis_evidence(details)

    try:
        if records:
            question = (
                f"Prepare a comprehensive evidence-grounded disease-to-target "
                f"research synthesis for {query.strip()}. Explain disease "
                "biology only to the extent supported by the supplied evidence. "
                "Compare the ranked targets and their association scores. "
                "Describe UniProt protein records, summarize each PubMed paper "
                "using its supplied abstract when available, and discuss "
                "Reactome, STRING, GO and KEGG findings, clinical candidates, "
                "disease treatments, evidence gaps, and limitations. Clearly "
                "distinguish target associations from causal mechanisms and "
                "clinical efficacy. State when disease biology or abstracts "
                "are missing. Use source labels such as [E1] for every "
                "evidence-based claim; do not invent facts or citations."
            )
            answer_result = answer_research_question(
                question, records, top_k=min(len(records), 50)
            )
            answer_text = answer_result.get(
                "answer", "No synthesis was generated."
            )
            validation = answer_result.get("citation_validation") or {}
            if validation.get("scope") == "label_integrity_only":
                answer_text += (
                    "\n\n**Citation limitation:** citation labels were "
                    "checked, but scientific claim support was not verified."
                )
        else:
            answer_text = (
                "Synthesis unavailable because the standard analysis returned "
                "no normalized evidence records. Review the analysis summary "
                "and retrieval warnings."
            )
    except Exception as exc:
        import traceback
        traceback.print_exc()
        answer_text = (
            f"Synthesis generation failed: {type(exc).__name__}: {exc}\n\n"
            "The standard analysis results remain available. "
            "Check the application terminal for the full traceback."
        )

    return (*standard_results, answer_text)



with gr.Blocks(title="Dyscovery | Disease-to-Target Research Lab") as demo:

    gr.HTML(HEADER)

    gr.Markdown(
        '<div class="section-heading">01 / DISEASE ANALYSIS</div>'
    )
    gr.Markdown(
        "Enter a disease to explore associated targets and supporting "
        "biological evidence."
    )

    with gr.Group(elem_classes=["panel", "search-panel"]):
        with gr.Row(elem_classes=["search-controls"]):
            disease_input = gr.Textbox(
                label="Disease or condition",
                placeholder="e.g. familial hypercholesterolemia, breast cancer, Alzheimer's disease",
                scale=5,
                min_width=260,
            )
            top_input = gr.Number(
                value=10,
                label="Number of top targets (1–20)",
                precision=0,
                minimum=1,
                maximum=20,
                step=1,
                scale=2,
                min_width=180,
            )
            analyze_button = gr.Button(
                "▶  Analyze",
                variant="primary",
                size="lg",
                scale=1,
                min_width=150,
            )

    summary_output = gr.Markdown()

    gr.Markdown(
        '<div class="section-heading">02 / EVIDENCE EXPLORER</div>'
    )

    with gr.Tabs():
        with gr.Tab("Ranked Targets"):
            gr.Markdown(
                "Targets are ordered by the prioritization pipeline. "
                "Association scores are ranking signals, not probabilities "
                "of clinical efficacy."
            )
            targets_output = gr.Dataframe(
                label="Disease-associated targets",
                interactive=False,
                wrap=True,
            )
        with gr.Tab("UniProt"):
            proteins_output = gr.Dataframe(
                label="Protein information",
                interactive=False,
                wrap=True,
            )
        with gr.Tab("PubMed"):
            publications_output = gr.Dataframe(
                label="Matching publications",
                interactive=False,
                wrap=True,
            )
        with gr.Tab("Reactome"):
            pathways_output = gr.Dataframe(
                label="Biological pathways",
                interactive=False,
                wrap=True,
            )
        with gr.Tab("STRING Network"):
            network_output = gr.HTML(
                label="Protein association network"
            )
        with gr.Tab("GO Enrichment"):
            go_output = gr.Dataframe(
                label="Gene Ontology enrichment",
                interactive=False,
                wrap=True,
            )
        with gr.Tab("KEGG Pathways"):
            kegg_output = gr.Dataframe(
                label="KEGG pathway enrichment",
                interactive=False,
                wrap=True,
            )
        with gr.Tab("Clinical Candidates"):
            drugs_output = gr.Dataframe(
                label="Target-associated clinical candidates",
                interactive=False,
                wrap=True,
            )
        with gr.Tab("Disease Treatments"):
            treatments_output = gr.Dataframe(
                label="Curated disease-level treatments",
                interactive=False,
                wrap=True,
            )
        with gr.Tab("Research Synthesis"):
            research_answer_output = gr.Markdown(
                label="Evidence-grounded research synthesis"
            )
        with gr.Tab("Detailed JSON"):
            json_output = gr.JSON(label="Complete analysis details")

    gr.Markdown(
        '<div class="section-heading">03 / EXPORT RESULTS</div>'
    )
    with gr.Row():
        csv_download = gr.File(label="Ranked targets CSV")
        json_download = gr.File(label="Detailed JSON report")

    outputs = [
        summary_output,
        targets_output,
        proteins_output,
        publications_output,
        pathways_output,
        network_output,
        go_output,
        kegg_output,
        drugs_output,
        treatments_output,
        json_output,
        csv_download,
        json_download,
        research_answer_output,
    ]

    analyze_button.click(
        fn=run_integrated_search,
        inputs=[disease_input, top_input],
        outputs=outputs,
    )

    disease_input.submit(
        fn=run_integrated_search,
        inputs=[disease_input, top_input],
        outputs=outputs,
    )

    gr.Markdown(
        """
        <div class="disclaimer">
          <strong>Research-use notice:</strong> This platform provides
          computational prioritization and evidence exploration. Association
          scores are not probabilities of treatment success. Clinical
          candidates and curated disease-level treatments are distinct
          categories. Results require independent review and validation and
          are not medical advice.
        </div>
        """
    )

    gr.HTML(ABOUT)


if __name__ == "__main__":
    demo.launch(
        server_name="0.0.0.0",
        server_port=int(__import__("os").environ.get("PORT", "7860")),
        theme=gr.themes.Base(
            primary_hue="cyan",
            secondary_hue="teal",
            neutral_hue="slate",
        ),
        css=CSS,
    )
