"""Tests for the provider layer, with a fake Azure client (no network).

The Azure surface is where most of the risk sits: deployment-name routing, the
reasoning-vs-chat parameter split, and two failure modes (truncation, content filter) that
must not be mistaken for results.
"""

import json
import os
from types import SimpleNamespace

import httpx
import openai
import pytest

from panspatial.orchestrator.backends import (
    AzureOpenAIBackend,
    DeploymentCapabilities,
    LLMError,
    load_dotenv,
    parse_deployment_map,
    resolve_backend,
)
from panspatial.orchestrator.client import OrchestratorClient
from panspatial.orchestrator.registry import PromptRegistry

AZURE_ENV = {
    "AZURE_OPENAI_API_KEY": "test-key",
    "AZURE_OPENAI_ENDPOINT": "https://test.openai.azure.com",
    "AZURE_OPENAI_API_VERSION": "2024-10-21",
    "AZURE_OPENAI_DEPLOYMENTS": "claude-opus-5=gpt5-prod,fast=gpt4o-mini-dev",
}


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    for key in list(AZURE_ENV) + [
        "PANSPATIAL_BACKEND",
        "AZURE_OPENAI_REASONING_DEPLOYMENTS",
        "AZURE_OPENAI_CHAT_DEPLOYMENTS",
    ]:
        monkeypatch.delenv(key, raising=False)


@pytest.fixture(scope="module")
def prompt():
    return PromptRegistry("prompts").render(
        "phase3.nk_spatial", {"adata_path": "x.h5ad", "platform": "Xenium"}
    )


def _bad_request(message: str) -> openai.BadRequestError:
    response = httpx.Response(400, request=httpx.Request("POST", "https://test.openai.azure.com"))
    return openai.BadRequestError(message, response=response, body={"error": {"code": "badRequest"}})


def _completion(text="ok", finish="stop", prompt_tokens=1200, cached=1024, completion_tokens=300):
    return SimpleNamespace(
        id="chatcmpl-test",
        choices=[SimpleNamespace(finish_reason=finish, message=SimpleNamespace(content=text))],
        usage=SimpleNamespace(
            prompt_tokens=prompt_tokens,
            completion_tokens=completion_tokens,
            prompt_tokens_details=SimpleNamespace(cached_tokens=cached),
            completion_tokens_details=SimpleNamespace(reasoning_tokens=128),
        ),
    )


class FakeAzureClient:
    """Records calls and replays a scripted sequence of responses or exceptions."""

    def __init__(self, *responses):
        self._responses = list(responses) or [_completion()]
        self.calls = []
        self.chat = SimpleNamespace(completions=SimpleNamespace(create=self._create))

    def _create(self, **kwargs):
        self.calls.append(kwargs)
        item = self._responses.pop(0) if len(self._responses) > 1 else self._responses[0]
        if isinstance(item, Exception):
            raise item
        return item


def make_backend(client, **kwargs):
    kwargs.setdefault("deployments", AZURE_ENV["AZURE_OPENAI_DEPLOYMENTS"])
    return AzureOpenAIBackend(client=client, **kwargs)


# ------------------------------------------------------------------- config


def test_deployment_map_accepts_both_forms():
    assert parse_deployment_map("a=dep-a,b=dep-b") == {"a": "dep-a", "b": "dep-b"}
    assert parse_deployment_map("solo") == {"solo": "solo"}
    assert parse_deployment_map("  a=dep-a ,  bare  ") == {"a": "dep-a", "bare": "bare"}
    assert parse_deployment_map(None) == {}


def test_empty_deployment_spec_is_refused_with_an_actionable_message():
    with pytest.raises(LLMError, match="AZURE_OPENAI_DEPLOYMENTS"):
        AzureOpenAIBackend(client=FakeAzureClient(), deployments="")


def test_dotenv_loads_but_never_overrides_the_shell(tmp_path, monkeypatch):
    env_file = tmp_path / ".env"
    env_file.write_text(
        "# comment\n"
        'AZURE_OPENAI_API_KEY="from-file"\n'
        "export AZURE_OPENAI_ENDPOINT='https://from-file.example'\n"
        "AZURE_OPENAI_API_VERSION=2024-10-21\n"
        "\n"
        "malformed-line\n"
    )
    monkeypatch.setenv("AZURE_OPENAI_API_KEY", "from-shell")
    loaded = load_dotenv(env_file)

    assert os.environ["AZURE_OPENAI_API_KEY"] == "from-shell", "shell must win over .env"
    assert os.environ["AZURE_OPENAI_ENDPOINT"] == "https://from-file.example"
    assert "AZURE_OPENAI_API_KEY" not in loaded


def test_missing_dotenv_is_not_an_error(tmp_path):
    assert load_dotenv(tmp_path / "nope.env") == {}


def test_azure_is_auto_selected_when_its_variables_are_present(monkeypatch):
    for k, v in AZURE_ENV.items():
        monkeypatch.setenv(k, v)
    assert resolve_backend().name == "azure-openai"


def test_explicit_backend_name_wins_over_auto_detection(monkeypatch):
    for k, v in AZURE_ENV.items():
        monkeypatch.setenv(k, v)
    assert resolve_backend("azure").name == "azure-openai"
    with pytest.raises(LLMError, match="unknown backend"):
        resolve_backend("bedrock")


def test_missing_azure_credentials_names_what_is_missing(monkeypatch):
    monkeypatch.setenv("AZURE_OPENAI_DEPLOYMENTS", "dep")
    with pytest.raises(LLMError, match="AZURE_OPENAI_ENDPOINT"):
        AzureOpenAIBackend()


# ------------------------------------------------------- deployment resolution


def test_mapped_model_routes_to_its_deployment():
    backend = make_backend(FakeAzureClient())
    assert backend.resolve_model("claude-opus-5") == "gpt5-prod"
    assert backend.resolve_model("fast") == "gpt4o-mini-dev"


def test_unmapped_model_falls_back_to_the_default_deployment():
    """Prompt templates name Anthropic models; they must still run on Azure unedited."""
    backend = make_backend(FakeAzureClient())
    assert backend.resolve_model("claude-sonnet-5") == "gpt5-prod"


# --------------------------------------------------------------- happy path


def test_completion_sends_the_anchor_first_and_normalizes_usage(prompt):
    client = FakeAzureClient(_completion(text="result text"))
    result = make_backend(client).complete(prompt)

    sent = client.calls[0]
    assert sent["model"] == "gpt5-prod"
    assert sent["messages"][0]["role"] == "system"
    assert sent["messages"][0]["content"] == prompt.system
    assert sent["messages"][1]["content"] == prompt.user

    assert result.text == "result text"
    assert result.provider == "azure-openai"
    assert result.model == "gpt5-prod"
    assert result.cached_input_tokens == 1024
    assert result.input_tokens == 1200 - 1024, "cached tokens must not be double-counted"
    assert result.output_tokens == 300
    assert result.extra["reasoning_tokens"] == 128


def test_azure_cost_is_left_unpriced_rather_than_guessed(prompt):
    """Azure pricing varies by region and agreement; a list price here would be wrong."""
    result = make_backend(FakeAzureClient()).complete(prompt)
    assert result.cost_usd is None


@pytest.mark.parametrize(
    "effort,expected", [("low", "low"), ("high", "high"), ("xhigh", "high"), ("max", "high")]
)
def test_five_level_effort_maps_onto_azures_three(effort, expected):
    client = FakeAzureClient()
    make_backend(client).complete(
        PromptRegistry("prompts").render(
            "phase3.nk_spatial", {"adata_path": "x.h5ad", "platform": "Xenium"}, effort=effort
        )
    )
    assert client.calls[0]["reasoning_effort"] == expected


def test_structured_output_uses_strict_json_schema(prompt):
    schema = {"type": "object", "properties": {"a": {"type": "string"}}, "required": ["a"],
              "additionalProperties": False}
    client = FakeAzureClient(_completion(text='{"a": "b"}'))
    result = make_backend(client).complete(prompt, schema=schema)

    fmt = client.calls[0]["response_format"]
    assert fmt["type"] == "json_schema"
    assert fmt["json_schema"]["strict"] is True
    assert fmt["json_schema"]["schema"] == schema
    assert result.parsed == {"a": "b"}


def test_invalid_structured_output_is_an_error_not_a_string(prompt):
    client = FakeAzureClient(_completion(text="I cannot produce JSON."))
    with pytest.raises(LLMError, match="not valid JSON"):
        make_backend(client).complete(prompt, schema={"type": "object"})


# ------------------------------------------------------- capability probing


def test_deployment_wanting_max_tokens_is_detected_and_retried(prompt):
    client = FakeAzureClient(
        _bad_request("Unrecognized request argument supplied: max_completion_tokens"),
        _completion(),
    )
    backend = make_backend(client)
    backend.complete(prompt)

    assert "max_completion_tokens" in client.calls[0]
    assert "max_tokens" in client.calls[1]
    assert client.calls[1]["max_tokens"] == prompt.max_tokens


def test_learned_capabilities_are_reused_on_the_next_call(prompt):
    client = FakeAzureClient(
        _bad_request("Unrecognized request argument supplied: max_completion_tokens"),
        _completion(),
    )
    backend = make_backend(client)
    backend.complete(prompt)
    n_after_first = len(client.calls)
    backend.complete(prompt)
    assert len(client.calls) == n_after_first + 1, "the probe must not repeat per call"
    assert "max_tokens" in client.calls[-1]


def test_deployment_rejecting_reasoning_effort_is_retried_without_it(prompt):
    client = FakeAzureClient(
        _bad_request("Unsupported parameter: 'reasoning_effort' is not supported with this model."),
        _completion(),
    )
    make_backend(client).complete(prompt)
    assert "reasoning_effort" in client.calls[0]
    assert "reasoning_effort" not in client.calls[1]


def test_strict_schema_falls_back_to_json_object(prompt):
    client = FakeAzureClient(
        _bad_request("Invalid parameter: 'response_format.json_schema' is not supported."),
        _completion(text='{"a": "b"}'),
    )
    schema = {"type": "object", "properties": {"a": {"type": "string"}}}
    result = make_backend(client).complete(prompt, schema=schema)

    assert client.calls[1]["response_format"] == {"type": "json_object"}
    assert "conforming to this schema" in client.calls[1]["messages"][1]["content"]
    assert result.parsed == {"a": "b"}


def test_env_override_skips_probing_for_a_chat_deployment(prompt, monkeypatch):
    monkeypatch.setenv("AZURE_OPENAI_CHAT_DEPLOYMENTS", "gpt5-prod")
    client = FakeAzureClient()
    make_backend(client).complete(prompt)
    assert "max_tokens" in client.calls[0]
    assert "reasoning_effort" not in client.calls[0]
    assert len(client.calls) == 1, "an explicit override should need no probe"


def test_an_unadaptable_400_is_raised_not_retried_forever(prompt):
    client = FakeAzureClient(_bad_request("The prompt exceeds the model's context length."))
    with pytest.raises(LLMError, match="context length"):
        make_backend(client).complete(prompt)
    assert len(client.calls) == 1


# ------------------------------------------------------------ failure modes


def test_truncated_response_raises_instead_of_looking_complete(prompt):
    client = FakeAzureClient(_completion(text="Half an aud", finish="length"))
    with pytest.raises(LLMError, match="truncated"):
        make_backend(client).complete(prompt)


def test_content_filter_during_generation_raises(prompt):
    client = FakeAzureClient(_completion(text="partial", finish="content_filter"))
    with pytest.raises(LLMError, match="content filter"):
        make_backend(client).complete(prompt)


def test_content_filter_on_input_is_not_reported_as_an_empty_result(prompt):
    response = httpx.Response(400, request=httpx.Request("POST", "https://test.openai.azure.com"))
    client = FakeAzureClient(
        openai.BadRequestError(
            "The response was filtered due to the prompt triggering Azure OpenAI's content "
            "management policy.",
            response=response,
            body={"error": {"code": "content_filter"}},
        )
    )
    with pytest.raises(LLMError, match="content filter"):
        make_backend(client).complete(prompt)


def test_empty_completion_explains_the_reasoning_budget_cause(prompt):
    client = FakeAzureClient(_completion(text=""))
    with pytest.raises(LLMError, match="consumed by reasoning"):
        make_backend(client).complete(prompt)


def test_missing_deployment_names_the_actual_cause(prompt):
    response = httpx.Response(404, request=httpx.Request("POST", "https://test.openai.azure.com"))
    client = FakeAzureClient(openai.NotFoundError("DeploymentNotFound", response=response, body=None))
    with pytest.raises(LLMError, match="deployment name, not model name"):
        make_backend(client).complete(prompt)


# ------------------------------------------------------- orchestrator wiring


def test_client_records_the_provider_and_resolved_deployment(prompt, tmp_path):
    client = OrchestratorClient(tmp_path / "runs", backend=make_backend(FakeAzureClient()))
    result = client.run(prompt)

    assert result.provider == "azure-openai"
    assert result.model == "gpt5-prod"
    assert result.cost_display == "unpriced"

    records = list((tmp_path / "runs").glob("*.json"))
    assert len(records) == 1
    record = json.loads(records[0].read_text())
    assert record["prompt"]["requested_model"] == "claude-opus-5"
    assert record["prompt"]["resolved_model"] == "gpt5-prod"
    assert record["prompt"]["provider"] == "azure-openai"
    assert record["prompt"]["fingerprint"] == prompt.fingerprint


def test_the_verified_results_gate_is_provider_independent():
    """Swapping providers must not weaken the manuscript gate."""
    from panspatial.orchestrator.registry import PromptError

    registry = PromptRegistry("prompts")
    values = {
        "core_finding": "CAF-NK axis",
        "results_table": "see table",
        "n_tumor_types": "12",
        "limitations": "observational",
    }
    with pytest.raises(PromptError, match="gated"):
        registry.render("phase5.title_abstract", values)
