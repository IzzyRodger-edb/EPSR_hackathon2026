"""Tests of ``mellor.llm.backend``: tool-call parsing, the scripted backend, the
backend factory, the Inference API mapping and the transformers backend with a fake
model (the real model only with ``MELLOR_GPU_TESTS=1``)."""

from __future__ import annotations

import json
import os
from types import SimpleNamespace

import pytest

from mellor.llm import backend as backend_module
from mellor.llm.backend import (
    DEFAULT_MODEL,
    AssistantMessage,
    HFInferenceBackend,
    ScriptedBackend,
    ToolCall,
    TransformersBackend,
    decode_arguments,
    make_backend,
    parse_tool_calls,
)

# --------------------------------------------------------------------------- parsing


def test_parse_single_tool_call_block():
    text = '<tool_call>\n{"name": "inspect_data", "arguments": {"files": null}}\n</tool_call>'
    prose, calls = parse_tool_calls(text)
    assert prose == ""
    assert len(calls) == 1
    assert calls[0].name == "inspect_data"
    assert calls[0].arguments == {"files": None}
    assert calls[0].id == "call_1"


def test_parse_several_blocks_with_prose_before():
    text = (
        "I will look at the data first.\n"
        '<tool_call>\n{"name": "inspect_data", "arguments": {}}\n</tool_call>\n'
        '<tool_call>\n{"name": "list_known_systems", "arguments": {}}\n</tool_call>'
    )
    prose, calls = parse_tool_calls(text)
    assert prose == "I will look at the data first."
    assert [call.name for call in calls] == ["inspect_data", "list_known_systems"]
    assert [call.id for call in calls] == ["call_1", "call_2"]


def test_parse_arguments_given_as_json_string():
    text = '<tool_call>{"name": "library_spec", "arguments": "{\\"system\\": \\"water\\"}"}</tool_call>'
    _, calls = parse_tool_calls(text)
    assert calls[0].arguments == {"system": "water"}


def test_parse_missing_arguments_means_empty_dict():
    _, calls = parse_tool_calls('<tool_call>{"name": "list_known_systems"}</tool_call>')
    assert calls[0].arguments == {}


def test_parse_parameters_alias():
    _, calls = parse_tool_calls(
        '<tool_call>{"name": "get_results", "parameters": {"section": "fit"}}</tool_call>'
    )
    assert calls[0].arguments == {"section": "fit"}


def test_malformed_block_stays_in_prose():
    text = 'Trying: <tool_call>\n{"name": "inspect_data", "arguments": {broken\n</tool_call> done'
    prose, calls = parse_tool_calls(text)
    assert calls == []
    assert "<tool_call>" in prose and "broken" in prose


def test_block_without_name_stays_in_prose():
    text = '<tool_call>{"arguments": {"a": 1}}</tool_call>'
    prose, calls = parse_tool_calls(text)
    assert calls == []
    assert prose == text


def test_non_object_arguments_is_malformed():
    text = '<tool_call>{"name": "x", "arguments": [1, 2]}</tool_call>'
    prose, calls = parse_tool_calls(text)
    assert calls == [] and prose == text


def test_unterminated_block_with_complete_json_is_accepted():
    text = (
        'Calling now.\n<tool_call>\n{"name": "finish", "arguments": {"message": "bye"}}'
    )
    prose, calls = parse_tool_calls(text)
    assert prose == "Calling now."
    assert calls[0].name == "finish" and calls[0].arguments == {"message": "bye"}


def test_fenced_json_fallback():
    text = 'Here is the call:\n```json\n{"name": "validate_spec", "arguments": {"spec": {"name": "w"}}}\n```\nThanks.'
    prose, calls = parse_tool_calls(text)
    assert calls[0].name == "validate_spec"
    assert calls[0].arguments == {"spec": {"name": "w"}}
    assert prose == "Here is the call:\n\nThanks."


def test_fenced_block_without_tool_call_is_left_alone():
    text = 'Example:\n```json\n{"density": 0.1}\n```'
    prose, calls = parse_tool_calls(text)
    assert calls == [] and prose == text


def test_bare_json_object_fallback():
    text = '{"name": "inspect_data", "arguments": {}}'
    prose, calls = parse_tool_calls(text)
    assert prose == "" and calls[0].name == "inspect_data"


def test_bare_json_inside_prose():
    text = (
        'Sure. {"name": "get_results", "arguments": {"section": "fit"}} Let me check.'
    )
    prose, calls = parse_tool_calls(text)
    assert calls[0].arguments == {"section": "fit"}
    assert prose == "Sure.  Let me check."


def test_openai_layout_inside_block():
    text = '<tool_call>{"type": "function", "function": {"name": "finish", "arguments": "{\\"message\\": \\"ok\\"}"}}</tool_call>'
    _, calls = parse_tool_calls(text)
    assert calls[0].name == "finish" and calls[0].arguments == {"message": "ok"}


def test_plain_text_has_no_calls():
    prose, calls = parse_tool_calls("The water model fits well (R-factor 0.02).")
    assert calls == [] and prose.startswith("The water model")
    assert parse_tool_calls("") == ("", [])


def test_think_blocks_are_stripped_from_prose():
    text = "<think>\nreasoning\n</think>\nFinal answer."
    prose, calls = parse_tool_calls(text)
    assert prose == "Final answer." and calls == []


def test_tool_calls_keep_the_block_order_and_prose_in_between():
    text = 'A<tool_call>{"name": "x", "arguments": {}}</tool_call>B<tool_call>{"name": "y", "arguments": {}}</tool_call>'
    prose, calls = parse_tool_calls(text)
    assert [c.name for c in calls] == ["x", "y"]
    assert prose == "AB"


def test_decode_arguments():
    assert decode_arguments(None) == {}
    assert decode_arguments("") == {}
    assert decode_arguments('{"a": 1}') == {"a": 1}
    assert decode_arguments({"a": 1}) == {"a": 1}
    with pytest.raises(ValueError):
        decode_arguments("[1]")
    with pytest.raises(ValueError):
        decode_arguments("{bad")


def test_assistant_message_to_message():
    message = AssistantMessage("hi", [ToolCall("x", {"a": 1}, "call_1")])
    entry = message.to_message()
    assert entry["role"] == "assistant" and entry["content"] == "hi"
    assert entry["tool_calls"][0]["function"] == {"name": "x", "arguments": {"a": 1}}
    assert "tool_calls" not in AssistantMessage("plain").to_message()


# --------------------------------------------------------------------------- scripted


def test_scripted_backend_sequence_and_exhaustion():
    backend = ScriptedBackend(
        [
            AssistantMessage("one"),
            "two",
            ("three", [ToolCall("inspect_data", {})]),
            {"name": "finish", "arguments": {"message": "done"}},
            '<tool_call>{"name": "get_results", "arguments": {"section": "fit"}}</tool_call>',
        ]
    )
    assert backend.name == "scripted"
    assert backend.chat([]).content == "one"
    assert backend.chat([]).content == "two"
    third = backend.chat([])
    assert third.content == "three" and third.tool_calls[0].name == "inspect_data"
    assert third.tool_calls[0].id == "call_1"
    assert backend.chat([]).tool_calls[0].arguments == {"message": "done"}
    assert not backend.exhausted
    fifth = backend.chat([])
    assert fifth.tool_calls[0].arguments == {"section": "fit"} and fifth.content == ""
    assert backend.exhausted
    empty = backend.chat([])
    assert empty.content == "" and empty.tool_calls == []
    assert backend.exhausted
    assert len(backend.requests) == 6


def test_scripted_backend_function():
    def reply(messages):
        return AssistantMessage(f"{len(messages)} messages")

    backend = ScriptedBackend(reply)
    assert backend.chat([{"role": "user", "content": "a"}]).content == "1 messages"
    assert not backend.exhausted


def test_scripted_backend_rejects_bad_reply():
    backend = ScriptedBackend([42])
    with pytest.raises(TypeError):
        backend.chat([])


# --------------------------------------------------------------------------- factory


def test_make_backend_scripted_and_unknown():
    backend = make_backend("scripted", replies=["hi"])
    assert isinstance(backend, ScriptedBackend)
    with pytest.raises(ValueError):
        make_backend("mystery")


def test_make_backend_auto_without_model_raises(monkeypatch):
    monkeypatch.setattr(backend_module, "_local_model_available", lambda: False)
    monkeypatch.delenv("HF_TOKEN", raising=False)
    monkeypatch.delenv("MELLOR_BACKEND", raising=False)
    with pytest.raises(RuntimeError, match="HF_TOKEN"):
        make_backend("auto")


def test_make_backend_auto_prefers_api_when_token_set(monkeypatch):
    monkeypatch.setattr(backend_module, "_local_model_available", lambda: False)
    monkeypatch.setenv("HF_TOKEN", "hf_dummy")
    monkeypatch.setenv("MELLOR_MODEL", "org/other-model")
    backend = make_backend(None)
    assert isinstance(backend, HFInferenceBackend)
    assert backend.name == "org/other-model"
    assert backend.token == "hf_dummy"


def test_make_backend_explicit_model_and_env(monkeypatch):
    monkeypatch.delenv("MELLOR_MODEL", raising=False)
    backend = make_backend("hf-api:org/model", max_new_tokens=100)
    assert isinstance(backend, HFInferenceBackend)
    assert backend.name == "org/model" and backend.max_tokens == 100
    monkeypatch.setenv("MELLOR_BACKEND", "hf-api")
    assert make_backend().name == DEFAULT_MODEL


def test_make_backend_transformers_is_lazy():
    backend = make_backend("transformers:org/tiny", max_tokens=64, temperature=0.0)
    assert isinstance(backend, TransformersBackend)
    assert backend.name == "org/tiny" and backend.max_new_tokens == 64
    assert not backend.loaded
    assert isinstance(backend.device_description, str)


# --------------------------------------------------------------------------- inference API


def _fake_response(content, tool_calls):
    return SimpleNamespace(
        choices=[
            SimpleNamespace(
                message=SimpleNamespace(content=content, tool_calls=tool_calls),
                finish_reason="stop",
            )
        ]
    )


def test_hf_backend_maps_tool_calls():
    response = _fake_response(
        None,
        [
            SimpleNamespace(
                id="abc",
                type="function",
                function=SimpleNamespace(
                    name="library_spec", arguments='{"system": "water"}'
                ),
            ),
            SimpleNamespace(
                id=None,
                type="function",
                function=SimpleNamespace(name="finish", arguments={"message": "x"}),
            ),
        ],
    )
    message = HFInferenceBackend.parse_response(response)
    assert message.content == ""
    assert [c.name for c in message.tool_calls] == ["library_spec", "finish"]
    assert (
        message.tool_calls[0].arguments == {"system": "water"}
        and message.tool_calls[0].id == "abc"
    )
    assert message.tool_calls[1].id == "call_2"


def test_hf_backend_dict_response_and_text_fallback():
    response = {
        "choices": [
            {
                "message": {
                    "content": '<tool_call>{"name": "finish", "arguments": {}}</tool_call>',
                    "tool_calls": None,
                }
            }
        ]
    }
    message = HFInferenceBackend.parse_response(response)
    assert message.tool_calls[0].name == "finish" and message.content == ""


def test_hf_backend_undecodable_arguments_stay_in_prose():
    response = _fake_response(
        "Calling.",
        [
            SimpleNamespace(
                id="1",
                type="function",
                function=SimpleNamespace(name="x", arguments="{bad"),
            )
        ],
    )
    message = HFInferenceBackend.parse_response(response)
    assert message.tool_calls == []
    assert "<tool_call>" in message.content and "{bad" in message.content


def test_hf_backend_chat_uses_client(monkeypatch):
    backend = HFInferenceBackend(
        "org/model", token="hf_x", max_tokens=50, temperature=0.1
    )
    seen = {}

    def chat_completion(**kwargs):
        seen.update(kwargs)
        return _fake_response("Answer text.", [])

    monkeypatch.setattr(
        backend, "_client", SimpleNamespace(chat_completion=chat_completion)
    )
    messages = [
        {"role": "system", "content": "s"},
        {
            "role": "assistant",
            "content": None,
            "tool_calls": [ToolCall("x", {"a": 1}, "call_1").to_message_entry()],
        },
        {"role": "tool", "name": "x", "tool_call_id": "call_1", "content": "{}"},
    ]
    tools = [
        {
            "type": "function",
            "function": {
                "name": "x",
                "description": "",
                "parameters": {"type": "object", "properties": {}},
            },
        }
    ]
    reply = backend.chat(messages, tools)
    assert reply.content == "Answer text." and reply.tool_calls == []
    assert seen["tools"] == tools and seen["tool_choice"] == "auto"
    assert seen["max_tokens"] == 50 and seen["temperature"] == 0.1
    sent = seen["messages"]
    assert sent[1]["content"] == ""
    assert sent[1]["tool_calls"][0]["function"]["arguments"] == '{"a": 1}'
    reply = backend.chat([{"role": "user", "content": "hi"}])
    assert "tools" not in seen or seen.get("tool_choice") == "auto"


# --------------------------------------------------------------------------- transformers (fake model)


def test_prepare_messages_makes_arguments_dicts():
    messages = [
        {"role": "user", "content": None},
        {
            "role": "assistant",
            "content": "",
            "tool_calls": [
                {
                    "id": "call_1",
                    "type": "function",
                    "function": {"name": "x", "arguments": '{"a": 1}'},
                },
                {
                    "id": "call_2",
                    "type": "function",
                    "function": {"name": "y", "arguments": "{bad"},
                },
            ],
        },
        {
            "role": "tool",
            "name": "x",
            "tool_call_id": "call_1",
            "content": {"ok": True},
        },
    ]
    prepared = TransformersBackend.prepare_messages(messages)
    assert prepared[0]["content"] == ""
    assert prepared[1]["tool_calls"][0]["function"]["arguments"] == {"a": 1}
    assert prepared[1]["tool_calls"][1]["function"]["arguments"] == {"raw": "{bad"}
    assert prepared[2]["content"] == '{"ok": true}'
    assert (
        messages[1]["tool_calls"][0]["function"]["arguments"] == '{"a": 1}'
    )  # untouched


torch = pytest.importorskip("torch")


class _FakeTokenizer:
    pad_token_id = None
    eos_token_id = 7

    def __init__(self, reply: str, n_prompt: int = 10):
        self.reply = reply
        self.n_prompt = n_prompt
        self.seen: dict = {}

    def apply_chat_template(
        self, messages, tools=None, add_generation_prompt=True, tokenize=False
    ):
        self.seen = {
            "messages": messages,
            "tools": tools,
            "add_generation_prompt": add_generation_prompt,
        }
        return "PROMPT"

    def __call__(self, prompt, return_tensors="pt", add_special_tokens=False):
        assert prompt == "PROMPT" and not add_special_tokens
        ids = torch.zeros((1, self.n_prompt), dtype=torch.long)
        return {"input_ids": ids, "attention_mask": torch.ones_like(ids)}

    def decode(self, ids, skip_special_tokens=True):
        assert skip_special_tokens and int(ids.shape[0]) == 5
        return self.reply


class _FakeModel:
    device = torch.device("cpu")

    def __init__(self):
        self.kwargs: dict = {}

    def generate(self, input_ids, attention_mask=None, **kwargs):
        self.kwargs = kwargs
        return torch.cat([input_ids, torch.full((1, 5), 3, dtype=torch.long)], dim=1)


def _fake_backend(
    reply: str, **options
) -> tuple[TransformersBackend, _FakeTokenizer, _FakeModel]:
    backend = TransformersBackend("org/fake", **options)
    tokenizer, model = _FakeTokenizer(reply), _FakeModel()
    backend._tokenizer, backend._model = tokenizer, model
    return backend, tokenizer, model


def test_transformers_chat_with_fake_model_parses_tool_call():
    reply = 'Let me check.\n<tool_call>\n{"name": "add", "arguments": {"a": 2, "b": 3}}\n</tool_call>'
    backend, tokenizer, model = _fake_backend(
        reply, temperature=0.7, top_p=0.9, top_k=5, max_new_tokens=33, seed=1
    )
    messages = [
        {"role": "user", "content": "add"},
        {
            "role": "assistant",
            "content": "",
            "tool_calls": [
                {
                    "id": "c1",
                    "type": "function",
                    "function": {"name": "add", "arguments": '{"a": 1}'},
                }
            ],
        },
        {"role": "tool", "name": "add", "tool_call_id": "c1", "content": "1"},
    ]
    tools = [
        {
            "type": "function",
            "function": {"name": "add", "description": "", "parameters": {}},
        }
    ]
    message = backend.chat(messages, tools)
    assert message.content == "Let me check." and message.raw == reply
    assert message.tool_calls[0].name == "add" and message.tool_calls[0].arguments == {
        "a": 2,
        "b": 3,
    }
    assert tokenizer.seen["tools"] == tools and tokenizer.seen["add_generation_prompt"]
    assert tokenizer.seen["messages"][1]["tool_calls"][0]["function"]["arguments"] == {
        "a": 1
    }
    assert model.kwargs["do_sample"] is True and model.kwargs["temperature"] == 0.7
    assert model.kwargs["top_p"] == 0.9 and model.kwargs["top_k"] == 5
    assert model.kwargs["max_new_tokens"] == 33 and model.kwargs["pad_token_id"] == 7
    assert backend.last_prompt_tokens == 10 and backend.last_new_tokens == 5


def test_transformers_greedy_when_temperature_zero():
    backend, _, model = _fake_backend("plain answer", temperature=0.0)
    message = backend.chat([{"role": "user", "content": "hi"}])
    assert message.content == "plain answer" and message.tool_calls == []
    assert model.kwargs["do_sample"] is False and "temperature" not in model.kwargs


def test_transformers_context_limit():
    backend, _, _ = _fake_backend("x", max_context_tokens=5)
    with pytest.raises(RuntimeError, match="10 tokens"):
        backend.chat([{"role": "user", "content": "hi"}])


def test_transformers_bad_dtype():
    backend = TransformersBackend("org/fake", dtype="float99")
    with pytest.raises(ValueError, match="dtype"):
        backend._resolve_dtype("cpu")


@pytest.mark.gpu
def test_transformers_real_model_calls_add_tool():
    if os.environ.get("MELLOR_GPU_TESTS") != "1":
        pytest.skip("set MELLOR_GPU_TESTS=1 to run the local model")
    backend = TransformersBackend(temperature=0.0, max_new_tokens=200)
    tools = [
        {
            "type": "function",
            "function": {
                "name": "add",
                "description": "Add two numbers and return the sum.",
                "parameters": {
                    "type": "object",
                    "properties": {"a": {"type": "number"}, "b": {"type": "number"}},
                    "required": ["a", "b"],
                },
            },
        }
    ]
    messages = [
        {
            "role": "system",
            "content": "You are a calculator. Always use the add tool for additions.",
        },
        {"role": "user", "content": "Use the add tool to add 2 and 3."},
    ]
    reply = backend.chat(messages, tools)
    assert reply.tool_calls, f"no tool call parsed from: {reply.raw!r}"
    call = reply.tool_calls[0]
    assert call.name == "add"
    assert {float(call.arguments["a"]), float(call.arguments["b"])} == {2.0, 3.0}
    assert "cuda" in backend.device_description or backend.device_description == "cpu"
    assert json.dumps(reply.to_message())  # serialisable for the transcript
