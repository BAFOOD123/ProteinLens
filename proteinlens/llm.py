"""A small interface to a language model, with a mock for tests and a response cache.

ProteinLens never needs an LLM: without an API key the template report is used.

* :class:`AnthropicClient` calls Claude through the official ``anthropic`` SDK.
  The API key comes from ``ANTHROPIC_API_KEY`` and the model from
  ``PROTEINLENS_MODEL`` (default ``claude-opus-5-5``).
* :class:`MockClient` returns scripted responses (used in tests and in the
  pipeline dry-run of the LLM experiment).
* :class:`CachedClient` stores every response on disk keyed by a hash of the
  request, so re-running an experiment costs nothing.
"""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

from proteinlens.data import cache_dir

DEFAULT_MODEL = "claude-opus-5-5"
# USD per million tokens (input, output), first-party API prices as of 2026-10.
PRICES = {
    "claude-opus-5-5": (4.00, 20.00),
    "claude-sonnet-5-5": (2.00, 10.00),
    "claude-haiku-5-5": (0.10, 0.50),
}


class LLMError(Exception):
    """The model could not produce a usable response."""


@dataclass
class LLMResponse:
    text: str
    model: str
    input_tokens: int = 0
    output_tokens: int = 0
    cached: bool = False

    @property
    def cost_usd(self) -> float:
        p_in, p_out = PRICES.get(self.model, (0.0, 0.0))
        return (self.input_tokens * p_in + self.output_tokens * p_out) / 1e6


class LLMClient(Protocol):
    model: str

    def generate(self, system: str, user: str, schema: dict[str, Any] | None = None) -> LLMResponse:
        """Return the model's reply. With ``schema`` the reply is JSON matching it."""
        ...


def api_key_available() -> bool:
    return bool(os.environ.get("ANTHROPIC_API_KEY"))


class AnthropicClient:
    def __init__(self, model: str | None = None, effort: str = "low", max_tokens: int = 2000):
        import anthropic  # imported here so the package works without network/API setup

        self.model = model or os.environ.get("PROTEINLENS_MODEL") or DEFAULT_MODEL
        self.effort = effort
        self.max_tokens = max_tokens
        self._client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY

    def generate(self, system: str, user: str, schema: dict[str, Any] | None = None) -> LLMResponse:
        import anthropic

        output_config: dict[str, Any] = {"effort": self.effort}
        if schema is not None:
            output_config["format"] = {"type": "json_schema", "schema": schema}
        try:
            resp = self._client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                system=system,
                messages=[{"role": "user", "content": user}],
                output_config=output_config,
            )
        except anthropic.APIStatusError as err:
            raise LLMError(f"API error {err.status_code}: {err.message}") from err
        except anthropic.APIConnectionError as err:
            raise LLMError(f"Could not reach the API: {err}") from err
        if resp.stop_reason == "refusal":
            raise LLMError("The model declined to answer (stop_reason=refusal).")
        if resp.stop_reason == "max_tokens":
            raise LLMError("The response was cut off (max_tokens).")
        text = next((b.text for b in resp.content if b.type == "text"), "")
        return LLMResponse(
            text=text,
            model=self.model,
            input_tokens=resp.usage.input_tokens,
            output_tokens=resp.usage.output_tokens,
        )


@dataclass
class MockClient:
    """Returns scripted responses in order (the last one repeats). Records every call."""

    responses: list[str]
    model: str = "mock"
    calls: list[dict[str, Any]] = field(default_factory=list)

    def generate(self, system: str, user: str, schema: dict[str, Any] | None = None) -> LLMResponse:
        self.calls.append({"system": system, "user": user, "schema": schema})
        text = self.responses[min(len(self.calls) - 1, len(self.responses) - 1)]
        return LLMResponse(text=text, model=self.model)


class CachedClient:
    """Wrap a client so identical requests are answered from disk."""

    def __init__(self, inner: LLMClient, folder: Path | None = None, tag: str = ""):
        self.inner = inner
        self.model = inner.model
        self.folder = folder or cache_dir() / "llm"
        self.folder.mkdir(parents=True, exist_ok=True)
        self.tag = tag  # e.g. "attempt2", so deliberate retries are cached separately

    def key(self, system: str, user: str, schema: dict[str, Any] | None) -> str:
        blob = json.dumps([self.model, self.tag, system, user, schema], sort_keys=True)
        return hashlib.sha256(blob.encode()).hexdigest()[:32]

    def generate(self, system: str, user: str, schema: dict[str, Any] | None = None) -> LLMResponse:
        path = self.folder / f"{self.key(system, user, schema)}.json"
        if path.exists():
            d = json.loads(path.read_text())
            return LLMResponse(d["text"], d["model"], d["input_tokens"], d["output_tokens"], True)
        r = self.inner.generate(system, user, schema)
        path.write_text(
            json.dumps(
                {
                    "text": r.text,
                    "model": r.model,
                    "input_tokens": r.input_tokens,
                    "output_tokens": r.output_tokens,
                }
            )
        )
        return r


def get_client(model: str | None = None, cached: bool = True) -> LLMClient | None:
    """The real client if an API key is set, otherwise None (callers use the template)."""
    if not api_key_available():
        return None
    client: LLMClient = AnthropicClient(model=model)
    return CachedClient(client) if cached else client
