"""LLM wrapper tests with a fake SDK client (no network, no API key needed)."""

from types import SimpleNamespace

import anthropic
import pytest

from proteinlens.llm import AnthropicClient, CachedClient, LLMError, MockClient, get_client


class FakeMessages:
    def __init__(self, stop_reason="end_turn"):
        self.kwargs = None
        self.stop_reason = stop_reason

    def create(self, **kwargs):
        self.kwargs = kwargs
        return SimpleNamespace(
            stop_reason=self.stop_reason,
            content=[SimpleNamespace(type="text", text='{"summary": "x", "same_fold": true}')],
            usage=SimpleNamespace(input_tokens=1000, output_tokens=200),
        )


def make_client(monkeypatch, stop_reason="end_turn"):
    fake = FakeMessages(stop_reason)
    monkeypatch.setattr(anthropic, "Anthropic", lambda: SimpleNamespace(messages=fake))
    return AnthropicClient(), fake


def test_request_shape_and_cost(monkeypatch):
    monkeypatch.setenv("PROTEINLENS_MODEL", "claude-haiku-5-5")
    client, fake = make_client(monkeypatch)
    r = client.generate("sys", "user", {"type": "object"})
    assert fake.kwargs["model"] == "claude-haiku-5-5"
    assert fake.kwargs["output_config"]["format"]["type"] == "json_schema"
    assert fake.kwargs["system"] == "sys"
    assert r.cost_usd == pytest.approx((1000 * 0.10 + 200 * 0.50) / 1e6)


def test_default_model(monkeypatch):
    monkeypatch.delenv("PROTEINLENS_MODEL", raising=False)
    client, _ = make_client(monkeypatch)
    assert client.model == "claude-opus-5-5"


def test_refusal_raises(monkeypatch):
    client, _ = make_client(monkeypatch, stop_reason="refusal")
    with pytest.raises(LLMError, match="declined"):
        client.generate("s", "u")


def test_cache_avoids_second_call(tmp_path):
    inner = MockClient(["first", "second"])
    c = CachedClient(inner, folder=tmp_path)
    assert c.generate("s", "u").text == "first"
    again = c.generate("s", "u")
    assert again.text == "first" and again.cached and len(inner.calls) == 1


def test_no_key_means_no_client(monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    assert get_client() is None
