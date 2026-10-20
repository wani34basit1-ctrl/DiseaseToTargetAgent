import json

import pytest

from services import tool_selector
from services.tool_selector import ToolSelectionError, select_tools


def mock_response(monkeypatch, value):
    monkeypatch.setattr(
        tool_selector,
        "generate_response",
        lambda **kwargs: value,
    )


def test_selects_approved_tool(monkeypatch):
    mock_response(
        monkeypatch,
        json.dumps({
            "tool_calls": [{
                "tool": "search_disease",
                "arguments": {"disease_name": "Alzheimer disease"},
            }]
        }),
    )
    plan = select_tools("Find Alzheimer disease")
    assert plan["tool_calls"][0]["tool"] == "search_disease"


def test_rejects_unknown_tool(monkeypatch):
    mock_response(
        monkeypatch,
        '{"tool_calls":[{"tool":"run_code","arguments":{}}]}',
    )
    with pytest.raises(ToolSelectionError):
        select_tools("Find a disease")


def test_rejects_too_many_calls(monkeypatch):
    call = {
        "tool": "search_disease",
        "arguments": {"disease_name": "Alzheimer disease"},
    }
    mock_response(
        monkeypatch,
        json.dumps({"tool_calls": [call, call, call, call]}),
    )
    with pytest.raises(ToolSelectionError, match="call limit"):
        select_tools("Research Alzheimer disease")


def test_rejects_unapproved_arguments(monkeypatch):
    mock_response(
        monkeypatch,
        json.dumps({
            "tool_calls": [{
                "tool": "search_disease",
                "arguments": {
                    "disease_name": "Alzheimer disease",
                    "execute_code": True,
                },
            }]
        }),
    )
    with pytest.raises(ToolSelectionError):
        select_tools("Find Alzheimer disease")


def test_rejects_missing_required_arguments(monkeypatch):
    mock_response(
        monkeypatch,
        '{"tool_calls":[{"tool":"search_disease","arguments":{}}]}',
    )
    with pytest.raises(ToolSelectionError, match="missing required"):
        select_tools("Find a disease")


def test_returns_empty_plan_when_model_selects_no_tools(monkeypatch):
    mock_response(monkeypatch, '{"tool_calls":[]}')
    assert select_tools("I have not decided what to research") == {
        "tool_calls": []
    }


def test_rejects_invalid_json(monkeypatch):
    mock_response(monkeypatch, "not JSON")
    with pytest.raises(ToolSelectionError, match="valid JSON"):
        select_tools("Find Alzheimer disease")


def test_rejects_empty_request():
    with pytest.raises(ValueError):
        select_tools("  ")
