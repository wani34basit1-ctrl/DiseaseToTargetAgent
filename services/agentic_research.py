"""Orchestrate bounded research tools and evidence-grounded answering."""

import hashlib

from agent.tool_executor import ResearchToolExecutor
from services.research_assistant import answer_research_question


class AgenticResearchService:
    """Plan tool calls, collect sourced evidence, and answer a research question."""

    def __init__(self, executor=None, answerer=None):
        self.executor = executor or ResearchToolExecutor()
        self.answerer = answerer or answer_research_question

    @staticmethod
    def _record(source, evidence_type, title, content, gene_symbol=None,
                disease_name=None, source_id=None, url=None, metadata=None):
        identity = "|".join([
            source, evidence_type, str(source_id or ""),
            str(gene_symbol or ""), str(title or ""),
        ])
        evidence_id = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:20]

        return {
            "evidence_id": evidence_id,
            "source": source,
            "evidence_type": evidence_type,
            "gene_symbol": gene_symbol,
            "disease_name": disease_name,
            "title": title or "Untitled evidence",
            "content": content or "",
            "source_id": str(source_id) if source_id is not None else None,
            "url": url,
            "metadata": metadata or {},
        }

    @classmethod
    def _collect_evidence(cls, execution):
        records = []
        seen = set()

        for item in execution["results"]:
            if item.get("status") != "completed":
                continue

            result = item.get("result") or {}
            tool = item.get("tool")

            # Enrichment results already contain normalized, source-linked records.
            for record in result.get("evidence", []) or []:
                if isinstance(record, dict) and record.get("evidence_id"):
                    if record["evidence_id"] not in seen:
                        records.append(record)
                        seen.add(record["evidence_id"])

            # Preserve PubMed metadata as individual evidence records.
            if tool == "search_literature":
                for paper in result.get("results", []) or []:
                    if not isinstance(paper, dict):
                        continue
                    title = paper.get("title") or "PubMed publication"
                    pmid = paper.get("pmid")
                    record = cls._record(
                        source="PubMed",
                        evidence_type="publication",
                        title=title,
                        content=(
                            title + (
                                "\nAbstract: " + paper["abstract"]
                                if paper.get("abstract") else ""
                            )
                        ),
                        gene_symbol=result.get("gene_symbol"),
                        disease_name=result.get("disease_name"),
                        source_id=pmid,
                        url=paper.get("url"),
                        metadata={
                            "publication_date": paper.get("pubdate"),
                            "abstract_available": bool(paper.get("abstract")),
                        },
                    )
                    if record["evidence_id"] not in seen:
                        records.append(record)
                        seen.add(record["evidence_id"])

            # Keep target-association evidence distinct from experimental evidence.
            if tool == "find_disease_targets":
                disease = result.get("disease") or {}
                for target in result.get("targets", []) or []:
                    if not isinstance(target, dict):
                        continue
                    symbol = target.get("symbol")
                    score = target.get("association_score")
                    record = cls._record(
                        source="Open Targets",
                        evidence_type="target_association",
                        title=f"Target association: {symbol or target.get('name', 'Unknown')}",
                        content=(
                            f"Open Targets associates {symbol or 'this target'} "
                            f"with {disease.get('name', 'the queried disease')}. "
                            f"Association score: {score if score is not None else 'not reported'}. "
                            "This score is not a probability of clinical efficacy."
                        ),
                        gene_symbol=symbol,
                        disease_name=disease.get("name"),
                        source_id=target.get("gene_id"),
                        url=(
                            "https://platform.opentargets.org/target/"
                            + str(target["gene_id"])
                            if target.get("gene_id")
                            else "https://platform.opentargets.org/"
                        ),
                        metadata={"association_score": score},
                    )
                    if record["evidence_id"] not in seen:
                        records.append(record)
                        seen.add(record["evidence_id"])

        return records

    def run(self, question, model=None):
        """Execute a bounded plan and answer using only collected evidence."""
        if not isinstance(question, str) or not question.strip():
            raise ValueError("Enter a research question.")

        execution = self.executor.run(question, model=model)

        if not isinstance(execution, dict):
            raise ValueError("The executor returned an invalid response.")

        results = execution.get("results")
        if not isinstance(results, list):
            raise ValueError("The executor response must contain a results list.")

        for item in results:
            if not isinstance(item, dict):
                raise ValueError("The executor returned a malformed result.")
            if item.get("status") not in {"completed", "failed"}:
                raise ValueError("The executor returned an invalid result status.")

        evidence = self._collect_evidence(execution)
        answer = self.answerer(question.strip(), evidence)

        if not isinstance(answer, dict):
            raise ValueError("The research answerer returned an invalid response.")

        citation_validation = answer.get("citation_validation")
        validation_scope = (
            citation_validation.get("validation_scope", "not_reported")
            if isinstance(citation_validation, dict)
            else "not_reported"
        )

        return {
            "question": question.strip(),
            "tool_execution": execution,
            "evidence_count": len(evidence),
            "evidence": evidence,
            "research_answer": answer,
            "citation_validation_scope": validation_scope,
            "citation_support_verified": False,
        }
