"""Chat backends: where Mellor's words come from.

Messages use the OpenAI-style layout that the Hugging Face chat templates and the
Inference API both understand::

    {"role": "system" | "user" | "assistant" | "tool", "content": "..."}
    assistant messages may carry "tool_calls": [{"id": ..., "type": "function",
        "function": {"name": ..., "arguments": {...}}}]
    tool results: {"role": "tool", "name": ..., "tool_call_id": ..., "content": "<json>"}

Tools are OpenAI-style function schemas::

    {"type": "function", "function": {"name": ..., "description": ..., "parameters": {json schema}}}

Backends:

* :class:`TransformersBackend`: a Hugging Face model run locally with ``transformers``
  (default ``Qwen/Qwen3-4B-Instruct-2507``; bf16 on the GPU, float32 on the CPU).
  Tool calls are rendered by the model's chat template (``tools=``) and parsed from
  ``<tool_call>{json}</tool_call>`` blocks (Qwen/Hermes format), with a fallback for a
  fenced JSON block containing ``{"name": ..., "arguments": ...}``.
* :class:`HFInferenceBackend`: the Hugging Face Inference API
  (``huggingface_hub.InferenceClient.chat_completion`` with ``tools``), token from
  ``HF_TOKEN``.
* :class:`ScriptedBackend`: replies from a list (or a function of the messages); for
  tests and for running the pipeline without a model.
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any, Protocol

__all__ = [
    "DEFAULT_MODEL",
    "ToolCall",
    "AssistantMessage",
    "ChatBackend",
    "TransformersBackend",
    "HFInferenceBackend",
    "ScriptedBackend",
    "make_backend",
    "parse_tool_calls",
    "decode_arguments",
]

DEFAULT_MODEL = "Qwen/Qwen3-4B-Instruct-2507"

_TOOL_CALL_BLOCK = re.compile(r"<tool_call>\s*(.*?)\s*</tool_call>", re.DOTALL)
_OPEN_TOOL_CALL = re.compile(r"<tool_call>\s*(.*)\Z", re.DOTALL)
_FENCED_BLOCK = re.compile(r"```(?:json|JSON)?[ \t]*\n(.*?)```", re.DOTALL)
_THINK_BLOCK = re.compile(r"<think>.*?</think>\s*", re.DOTALL)


@dataclass
class ToolCall:
    name: str
    arguments: dict[str, Any]
    id: str = ""

    def to_message_entry(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "type": "function",
            "function": {"name": self.name, "arguments": self.arguments},
        }


@dataclass
class AssistantMessage:
    content: str = ""
    tool_calls: list[ToolCall] = field(default_factory=list)
    raw: str = ""  # the model's full text before parsing (for the transcript)

    def to_message(self) -> dict[str, Any]:
        """The assistant message to append to the conversation."""
        message: dict[str, Any] = {"role": "assistant", "content": self.content}
        if self.tool_calls:
            message["tool_calls"] = [
                call.to_message_entry() for call in self.tool_calls
            ]
        return message


class ChatBackend(Protocol):
    name: str

    def chat(
        self, messages: list[dict[str, Any]], tools: list[dict[str, Any]] | None = None
    ) -> AssistantMessage: ...


# --------------------------------------------------------------------------- parsing
def decode_arguments(value: Any) -> dict[str, Any]:
    """The ``arguments`` of a tool call as a dict.

    ``None`` means no arguments; a JSON string is decoded. Raises ``ValueError`` when
    the value is not (or does not decode to) a JSON object.
    """
    if value is None:
        return {}
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return {}
        try:
            value = json.loads(text)
        except json.JSONDecodeError as error:
            raise ValueError(f"arguments are not valid JSON: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"arguments must be a JSON object, not {type(value).__name__}")
    return value


def _tool_call_from_json(text: str) -> ToolCall | None:
    """A :class:`ToolCall` from the JSON text of one block, or None if malformed."""
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return None
    return _tool_call_from_object(data)


def _tool_call_from_object(data: Any) -> ToolCall | None:
    if not isinstance(data, dict):
        return None
    if "function" in data and isinstance(data["function"], dict):  # OpenAI layout
        data = data["function"]
    name = data.get("name")
    if not isinstance(name, str) or not name.strip():
        return None
    arguments = data.get("arguments", data.get("parameters"))
    try:
        return ToolCall(name=name.strip(), arguments=decode_arguments(arguments))
    except ValueError:
        return None


def _scan_bare_objects(text: str) -> list[tuple[int, int, ToolCall]]:
    """``(start, end, call)`` of every bare JSON object with name/arguments in ``text``."""
    decoder = json.JSONDecoder()
    found: list[tuple[int, int, ToolCall]] = []
    position = 0
    while True:
        start = text.find("{", position)
        if start < 0:
            break
        try:
            data, end = decoder.raw_decode(text, start)
        except json.JSONDecodeError:
            position = start + 1
            continue
        call = (
            _tool_call_from_object(data)
            if isinstance(data, dict) and "name" in data
            else None
        )
        if call is not None and ("arguments" in data or "parameters" in data):
            found.append((start, end, call))
            position = end
        else:
            position = start + 1
    return found


def parse_tool_calls(text: str) -> tuple[str, list[ToolCall]]:
    """Split a model reply into its prose and its tool calls.

    Recognises ``<tool_call>{...}</tool_call>`` blocks and, failing that, one fenced
    ```` ```json ```` block or bare JSON object with ``name`` and ``arguments`` keys.
    ``arguments`` given as a JSON string are decoded. Unparsable blocks are left in
    the prose. Ids are ``call_1``, ``call_2``...
    """
    if not text:
        return "", []
    prose = _THINK_BLOCK.sub("", text)
    calls: list[ToolCall] = []

    def take(match: re.Match[str]) -> str:
        call = _tool_call_from_json(match.group(1))
        if call is None:
            return match.group(0)
        calls.append(call)
        return ""

    prose = _TOOL_CALL_BLOCK.sub(take, prose)
    if "<tool_call>" in prose and "</tool_call>" not in prose:
        # the model stopped before closing the block (end of generation): accept
        # the block when its JSON is complete
        prose = _OPEN_TOOL_CALL.sub(take, prose)
    if not calls:
        prose = _FENCED_BLOCK.sub(take, prose)
    if not calls:
        spans = _scan_bare_objects(prose)
        if spans:
            pieces = []
            last = 0
            for start, end, call in spans:
                pieces.append(prose[last:start])
                calls.append(call)
                last = end
            pieces.append(prose[last:])
            prose = "".join(pieces)
    for index, call in enumerate(calls, start=1):
        call.id = f"call_{index}"
    return prose.strip(), calls


# --------------------------------------------------------------------------- local model
class TransformersBackend:
    """A local Hugging Face model (loaded lazily on the first call)."""

    name: str

    def __init__(
        self,
        model_id: str = DEFAULT_MODEL,
        *,
        device: str = "auto",
        dtype: str = "auto",
        max_new_tokens: int = 2048,
        temperature: float = 0.7,
        top_p: float = 0.8,
        top_k: int = 20,
        max_context_tokens: int = 24_000,
        seed: int | None = None,
    ) -> None:
        self.name = model_id
        self.model_id = model_id
        self.device = device
        self.dtype = dtype
        self.max_new_tokens = int(max_new_tokens)
        self.temperature = float(temperature)
        self.top_p = float(top_p)
        self.top_k = int(top_k)
        self.max_context_tokens = int(max_context_tokens)
        self.seed = seed
        self._tokenizer: Any = None
        self._model: Any = None
        self.last_prompt_tokens: int = 0
        self.last_new_tokens: int = 0

    # ---- device handling
    def _resolve_device(self) -> str:
        if self.device != "auto":
            return self.device
        try:
            import torch
        except ImportError:
            return "cpu"
        return "cuda" if torch.cuda.is_available() else "cpu"

    def _resolve_dtype(self, device: str) -> Any:
        import torch

        if self.dtype == "auto":
            return torch.bfloat16 if device.startswith("cuda") else torch.float32
        try:
            return getattr(torch, self.dtype)
        except AttributeError as error:
            raise ValueError(
                f"unknown dtype {self.dtype!r}; use 'auto', 'bfloat16', 'float16' or 'float32'"
            ) from error

    @property
    def device_description(self) -> str:
        """e.g. 'cuda:0 (NVIDIA GeForce RTX 5070 Ti Laptop GPU, 12.2 GB)' or 'cpu'."""
        try:
            import torch
        except ImportError:
            return "unavailable (torch is not installed)"
        device = self._resolve_device()
        if device.startswith("cuda") and torch.cuda.is_available():
            index = int(device.split(":", 1)[1]) if ":" in device else 0
            properties = torch.cuda.get_device_properties(index)
            memory_gb = properties.total_memory / 1e9
            return f"cuda:{index} ({properties.name}, {memory_gb:.1f} GB)"
        return "cpu"

    @property
    def loaded(self) -> bool:
        return self._model is not None and self._tokenizer is not None

    def load(self) -> None:
        """Download (Hugging Face cache) and load the model and tokenizer."""
        if self.loaded:
            return
        try:
            import torch
            from transformers import AutoModelForCausalLM, AutoTokenizer
        except ImportError as error:
            raise RuntimeError(
                "the transformers backend needs torch and transformers: install them "
                "with `uv sync --extra local` (or use the hf-api backend with HF_TOKEN)"
            ) from error
        device = self._resolve_device()
        dtype = self._resolve_dtype(device)
        self._tokenizer = AutoTokenizer.from_pretrained(self.model_id)
        kwargs: dict[str, Any] = {"dtype": dtype}
        if device.startswith("cuda"):
            if not torch.cuda.is_available():
                raise RuntimeError(
                    f"device {device!r} requested but CUDA is not available"
                )
            kwargs["device_map"] = device
        model: Any = AutoModelForCausalLM.from_pretrained(self.model_id, **kwargs)
        if not device.startswith("cuda"):
            model = model.to(device)
        model.eval()
        self._model = model

    # ---- messages
    @staticmethod
    def prepare_messages(messages: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """Messages as the chat template wants them: string contents, tool-call
        arguments as dicts (the Qwen3 template renders them with ``tojson``)."""
        prepared = []
        for message in messages:
            copy = dict(message)
            content = copy.get("content")
            if content is None:
                copy["content"] = ""
            elif not isinstance(content, str):
                copy["content"] = json.dumps(content, ensure_ascii=False)
            if copy.get("role") == "assistant" and copy.get("tool_calls"):
                calls = []
                for entry in copy["tool_calls"]:
                    entry = json.loads(json.dumps(entry))  # deep copy of plain data
                    function = entry.get("function", entry)
                    try:
                        function["arguments"] = decode_arguments(
                            function.get("arguments")
                        )
                    except ValueError:
                        function["arguments"] = {"raw": str(function.get("arguments"))}
                    calls.append(entry)
                copy["tool_calls"] = calls
            prepared.append(copy)
        return prepared

    def chat(
        self, messages: list[dict[str, Any]], tools: list[dict[str, Any]] | None = None
    ) -> AssistantMessage:
        """``apply_chat_template(messages, tools=tools, add_generation_prompt=True)``,
        ``generate`` (sampling with the given parameters; greedy when temperature is 0),
        decode only the new tokens, :func:`parse_tool_calls`. Raises ``RuntimeError``
        with a clear message if the prompt exceeds ``max_context_tokens``."""
        import torch

        self.load()
        tokenizer, model = self._tokenizer, self._model
        prompt = tokenizer.apply_chat_template(
            self.prepare_messages(messages),
            tools=tools or None,
            add_generation_prompt=True,
            tokenize=False,
        )
        encoded = tokenizer(prompt, return_tensors="pt", add_special_tokens=False)
        input_ids = encoded["input_ids"]
        n_prompt = int(input_ids.shape[1])
        self.last_prompt_tokens = n_prompt
        if n_prompt > self.max_context_tokens:
            raise RuntimeError(
                f"the prompt has {n_prompt} tokens, more than the limit of "
                f"{self.max_context_tokens}: shorten the conversation (start a new "
                "chat, or ask for smaller tool results)"
            )
        device = model.device
        inputs = {key: value.to(device) for key, value in encoded.items()}
        generate_kwargs: dict[str, Any] = {"max_new_tokens": self.max_new_tokens}
        pad_token_id = getattr(tokenizer, "pad_token_id", None)
        if pad_token_id is None:
            pad_token_id = getattr(tokenizer, "eos_token_id", None)
        if pad_token_id is not None:
            generate_kwargs["pad_token_id"] = pad_token_id
        if self.temperature > 0:
            generate_kwargs.update(
                do_sample=True,
                temperature=self.temperature,
                top_p=self.top_p,
                top_k=self.top_k,
            )
        else:
            generate_kwargs["do_sample"] = False
        if self.seed is not None:
            torch.manual_seed(self.seed)
        with torch.no_grad():
            output = model.generate(**inputs, **generate_kwargs)
        new_tokens = output[0, n_prompt:]
        self.last_new_tokens = int(new_tokens.shape[0])
        text = tokenizer.decode(new_tokens, skip_special_tokens=True)
        content, calls = parse_tool_calls(text)
        return AssistantMessage(content=content, tool_calls=calls, raw=text)


# --------------------------------------------------------------------------- inference API
class HFInferenceBackend:
    """The Hugging Face Inference API (serverless or dedicated endpoint)."""

    name: str

    def __init__(
        self,
        model_id: str = DEFAULT_MODEL,
        *,
        token: str | None = None,
        max_tokens: int = 2048,
        temperature: float = 0.7,
        timeout: float = 120.0,
    ) -> None:
        from huggingface_hub import InferenceClient

        self.name = model_id
        self.model_id = model_id
        self.max_tokens = int(max_tokens)
        self.temperature = float(temperature)
        self.timeout = float(timeout)
        self.token = token or os.environ.get("HF_TOKEN") or None
        self._client = InferenceClient(
            model=model_id, token=self.token, timeout=self.timeout
        )

    @staticmethod
    def prepare_messages(messages: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """OpenAI layout for the API: tool-call arguments as JSON strings."""
        prepared = []
        for message in messages:
            copy = dict(message)
            if copy.get("content") is None:
                copy["content"] = ""
            if copy.get("role") == "assistant" and copy.get("tool_calls"):
                calls = []
                for entry in copy["tool_calls"]:
                    entry = json.loads(json.dumps(entry))
                    function = entry.get("function", entry)
                    arguments = function.get("arguments")
                    if not isinstance(arguments, str):
                        function["arguments"] = json.dumps(arguments or {})
                    entry.setdefault("type", "function")
                    calls.append(entry)
                copy["tool_calls"] = calls
            prepared.append(copy)
        return prepared

    def chat(
        self, messages: list[dict[str, Any]], tools: list[dict[str, Any]] | None = None
    ) -> AssistantMessage:
        kwargs: dict[str, Any] = {
            "messages": self.prepare_messages(messages),
            "max_tokens": self.max_tokens,
            "temperature": self.temperature,
        }
        if tools:
            kwargs["tools"] = tools
            kwargs["tool_choice"] = "auto"
        response = self._client.chat_completion(**kwargs)
        return self.parse_response(response)

    @staticmethod
    def parse_response(response: Any) -> AssistantMessage:
        """Map a ``chat_completion`` response (object or dict) to an :class:`AssistantMessage`."""
        choices = _get(response, "choices") or []
        if not choices:
            raise RuntimeError("the Inference API returned no choices")
        message = _get(choices[0], "message")
        content = _get(message, "content") or ""
        calls: list[ToolCall] = []
        unparsed: list[str] = []
        for index, entry in enumerate(_get(message, "tool_calls") or [], start=1):
            function = _get(entry, "function")
            name = _get(function, "name") or ""
            raw_arguments = _get(function, "arguments")
            try:
                arguments = decode_arguments(raw_arguments)
            except ValueError:
                unparsed.append(json.dumps({"name": name, "arguments": raw_arguments}))
                continue
            call_id = _get(entry, "id") or f"call_{index}"
            calls.append(ToolCall(name=str(name), arguments=arguments, id=str(call_id)))
        raw = content
        if unparsed:  # malformed calls stay in the prose, as parse_tool_calls does
            raw = "\n".join(
                [content] + [f"<tool_call>\n{u}\n</tool_call>" for u in unparsed]
            )
            content = raw
        if not calls and content:
            # some providers leave Qwen-style <tool_call> blocks in the text
            content, extra = parse_tool_calls(content)
            calls.extend(extra)
        return AssistantMessage(content=content.strip(), tool_calls=calls, raw=raw)


def _get(obj: Any, key: str) -> Any:
    if obj is None:
        return None
    if isinstance(obj, dict):
        return obj.get(key)
    return getattr(obj, key, None)


# --------------------------------------------------------------------------- scripted
class ScriptedBackend:
    """Replies from a script: a sequence of :class:`AssistantMessage` (or strings, or
    ``(content, [ToolCall...])`` tuples), or a function ``messages -> AssistantMessage``.
    When the sequence is exhausted it returns an empty assistant message."""

    name: str = "scripted"

    def __init__(
        self,
        replies: Sequence[Any] | Callable[[list[dict[str, Any]]], AssistantMessage],
    ) -> None:
        self._function: Callable[[list[dict[str, Any]]], Any] | None
        self._replies: list[Any]
        if callable(replies):
            self._function = replies
            self._replies = []
        else:
            self._function = None
            self._replies = list(replies)
        self._index = 0
        self.requests: list[list[dict[str, Any]]] = []  # the messages seen at each call

    @property
    def exhausted(self) -> bool:
        return self._function is None and self._index >= len(self._replies)

    @staticmethod
    def coerce(reply: Any) -> AssistantMessage:
        """Any accepted reply form as an :class:`AssistantMessage`."""
        if isinstance(reply, AssistantMessage):
            for index, existing in enumerate(reply.tool_calls, start=1):
                if not existing.id:
                    existing.id = f"call_{index}"
            if not reply.raw:
                reply.raw = reply.content
            return reply
        if reply is None:
            return AssistantMessage()
        if isinstance(reply, str):
            content, calls = parse_tool_calls(reply)
            return AssistantMessage(content=content, tool_calls=calls, raw=reply)
        if isinstance(reply, ToolCall):
            reply = ("", [reply])
        if isinstance(reply, dict):
            call = _tool_call_from_object(reply)
            if call is None:
                raise TypeError(f"a dict reply needs 'name' and 'arguments': {reply!r}")
            reply = ("", [call])
        if isinstance(reply, tuple) and len(reply) == 2:
            content, calls = reply
            tool_calls = []
            for index, item in enumerate(calls or [], start=1):
                if isinstance(item, ToolCall):
                    call = item
                else:
                    call = _tool_call_from_object(item)
                    if call is None:
                        raise TypeError(f"bad tool call in scripted reply: {item!r}")
                if not call.id:
                    call.id = f"call_{index}"
                tool_calls.append(call)
            message = AssistantMessage(
                content=str(content or ""), tool_calls=tool_calls
            )
            message.raw = message.content
            return message
        raise TypeError(f"unsupported scripted reply: {reply!r}")

    def chat(
        self, messages: list[dict[str, Any]], tools: list[dict[str, Any]] | None = None
    ) -> AssistantMessage:
        self.requests.append([dict(message) for message in messages])
        if self._function is not None:
            return self.coerce(self._function(messages))
        if self._index >= len(self._replies):
            return AssistantMessage()
        reply = self._replies[self._index]
        self._index += 1
        return self.coerce(reply)


# --------------------------------------------------------------------------- factory
def _local_model_available() -> bool:
    """True when torch and transformers can be imported (without importing them)."""
    return all(
        importlib.util.find_spec(name) is not None for name in ("torch", "transformers")
    )


def make_backend(spec: str | None = None, **options: Any) -> ChatBackend:
    """A backend from a short string: ``"transformers[:model]"``, ``"hf-api[:model]"``,
    ``"scripted"`` or ``"auto"`` (transformers when torch is importable, else hf-api
    when ``HF_TOKEN`` is set, else an error). ``None`` reads ``MELLOR_BACKEND`` (default
    ``auto``); ``MELLOR_MODEL`` overrides the default model id."""
    spec = (spec or os.environ.get("MELLOR_BACKEND") or "auto").strip()
    kind, _, model = spec.partition(":")
    kind = kind.strip().lower()
    model = (
        options.pop("model", None)
        or options.pop("model_id", None)
        or model.strip()
        or os.environ.get("MELLOR_MODEL")
        or DEFAULT_MODEL
    )
    if kind == "auto":
        if _local_model_available():
            kind = "transformers"
        elif os.environ.get("HF_TOKEN"):
            kind = "hf-api"
        else:
            raise RuntimeError(
                "no language-model backend is available: install torch and transformers "
                "(`uv sync --extra local`) to run the model locally, or set HF_TOKEN to "
                "use the Hugging Face Inference API (--backend hf-api). Without a model, "
                "`mellor run --system KEY` builds the specification from the library."
            )
    if kind == "transformers":
        if "max_tokens" in options:
            options["max_new_tokens"] = options.pop("max_tokens")
        return TransformersBackend(model, **options)
    if kind in ("hf-api", "hf_api", "hfapi", "api", "inference"):
        if "max_new_tokens" in options:
            options["max_tokens"] = options.pop("max_new_tokens")
        return HFInferenceBackend(model, **options)
    if kind == "scripted":
        return ScriptedBackend(options.pop("replies", []))
    raise ValueError(
        f"unknown backend {spec!r}; use 'transformers[:model]', 'hf-api[:model]', "
        "'scripted' or 'auto'"
    )
