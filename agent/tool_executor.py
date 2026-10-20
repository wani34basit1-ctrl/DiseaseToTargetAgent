"""Bounded executor for approved biological research tools."""

from agent.research_tools import ResearchTools, ResearchToolError
from services.tool_selector import select_tools, MAX_TOOL_CALLS


class ResearchToolExecutor:
    """Select and execute a limited number of approved research tools."""

    def __init__(self, research_tools=None, selector=None):
        self.research_tools = research_tools or ResearchTools()
        self.selector = selector or select_tools

    def run(self, user_request, model=None):
        if not isinstance(user_request, str) or not user_request.strip():
            raise ResearchToolError("Enter a research request.")

        plan = self.selector(user_request, model=model)

        if not isinstance(plan, dict) or set(plan) != {"tool_calls"}:
            raise ResearchToolError("The tool selector returned an invalid plan.")

        calls = plan["tool_calls"]

        if not isinstance(calls, list):
            raise ResearchToolError("The tool plan must contain a list of calls.")

        if len(calls) > MAX_TOOL_CALLS:
            raise ResearchToolError("The tool plan exceeds the call limit.")

        # Validate the entire plan before executing any tool.
        for call in calls:
            if not isinstance(call, dict) or set(call) != {"tool", "arguments"}:
                raise ResearchToolError("A tool call has an invalid structure.")

            if not isinstance(call["tool"], str) or not call["tool"].strip():
                raise ResearchToolError("A tool call has an invalid tool name.")

            if not isinstance(call["arguments"], dict):
                raise ResearchToolError("Tool arguments must be an object.")

        results = []

        for index, call in enumerate(calls, start=1):
            tool_name = call["tool"]
            arguments = call["arguments"]

            try:
                result = self.research_tools.execute(tool_name, arguments)
                results.append({
                    "call_index": index,
                    "tool": tool_name,
                    "status": "completed",
                    "result": result,
                })
            except (ResearchToolError, ValueError, TypeError) as exc:
                results.append({
                    "call_index": index,
                    "tool": tool_name,
                    "status": "failed",
                    "error": str(exc),
                })
            except Exception:
                results.append({
                    "call_index": index,
                    "tool": tool_name,
                    "status": "failed",
                    "error": "The research tool failed unexpectedly.",
                })

        return {
            "request": user_request.strip(),
            "planned_call_count": len(calls),
            "completed_call_count": sum(
                item["status"] == "completed" for item in results
            ),
            "results": results,
        }
