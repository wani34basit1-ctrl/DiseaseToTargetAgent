from agent.tool_executor import ResearchToolExecutor


class FakeResearchTools:
    def __init__(self):
        self.calls = []

    def execute(self, tool_name, arguments):
        self.calls.append((tool_name, arguments))
        if tool_name == "search_disease":
            return {"results": [{"name": arguments["disease_name"]}]}
        raise ValueError("Mock tool failure")


def test_executor_runs_approved_plan():
    tools = FakeResearchTools()
    executor = ResearchToolExecutor(
        research_tools=tools,
        selector=lambda request, model=None: {
            "tool_calls": [{
                "tool": "search_disease",
                "arguments": {"disease_name": "Alzheimer disease"},
            }]
        },
    )

    result = executor.run("Find Alzheimer disease")

    assert result["planned_call_count"] == 1
    assert result["completed_call_count"] == 1
    assert result["results"][0]["status"] == "completed"
    assert len(tools.calls) == 1


def test_executor_records_tool_failure_and_continues():
    tools = FakeResearchTools()
    executor = ResearchToolExecutor(
        research_tools=tools,
        selector=lambda request, model=None: {
            "tool_calls": [
                {
                    "tool": "unknown_mock_tool",
                    "arguments": {},
                },
                {
                    "tool": "search_disease",
                    "arguments": {"disease_name": "Alzheimer disease"},
                },
            ]
        },
    )

    result = executor.run("Research Alzheimer disease")

    assert result["results"][0]["status"] == "failed"
    assert result["results"][1]["status"] == "completed"
    assert result["completed_call_count"] == 1


def test_executor_returns_empty_result_for_empty_plan():
    executor = ResearchToolExecutor(
        research_tools=FakeResearchTools(),
        selector=lambda request, model=None: {"tool_calls": []},
    )

    result = executor.run("General research question")

    assert result["planned_call_count"] == 0
    assert result["completed_call_count"] == 0
    assert result["results"] == []


def test_executor_rejects_malformed_plan_before_execution():
    tools = FakeResearchTools()
    executor = ResearchToolExecutor(
        research_tools=tools,
        selector=lambda request, model=None: {
            "tool_calls": [
                {
                    "tool": "search_disease",
                    "arguments": {"disease_name": "Alzheimer disease"},
                },
                {"tool": "search_disease"},
            ]
        },
    )

    try:
        executor.run("Research Alzheimer disease")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected malformed plan to be rejected")

    assert tools.calls == []


def test_executor_rejects_plan_over_call_limit_before_execution():
    from services.tool_selector import MAX_TOOL_CALLS

    tools = FakeResearchTools()
    executor = ResearchToolExecutor(
        research_tools=tools,
        selector=lambda request, model=None: {
            "tool_calls": [
                {"tool": "search_disease", "arguments": {}}
                for _ in range(MAX_TOOL_CALLS + 1)
            ]
        },
    )

    try:
        executor.run("Research a disease")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected excessive tool calls to be rejected")

    assert tools.calls == []


def test_executor_rejects_empty_request_before_selector():
    tools = FakeResearchTools()
    selector_calls = []

    def selector(request, model=None):
        selector_calls.append(request)
        return {"tool_calls": []}

    executor = ResearchToolExecutor(
        research_tools=tools,
        selector=selector,
    )

    try:
        executor.run("   ")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected empty request to be rejected")

    assert selector_calls == []
    assert tools.calls == []


def test_executor_rejects_non_list_tool_calls():
    tools = FakeResearchTools()
    executor = ResearchToolExecutor(
        research_tools=tools,
        selector=lambda request, model=None: {"tool_calls": "not-a-list"},
    )

    try:
        executor.run("Research a disease")
    except ValueError:
        pass
    else:
        raise AssertionError("Expected invalid tool call list to be rejected")

    assert tools.calls == []
