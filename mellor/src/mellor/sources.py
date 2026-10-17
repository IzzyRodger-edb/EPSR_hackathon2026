"""Finding comparable experimental data.

Three sources, all without credentials:

* **local**: ``.mint01`` files in this repository (``experimental_data/``,
  ``validation/``, ``EPSR-test/``) and in folders the user names;
* **dissolve**: the examples of the Dissolve project on GitHub
  (``disorderedmaterials/dissolve``, ``examples/<case>/data/*.mint01`` with the
  ``input.txt`` that describes the system) through the GitHub contents API;
* **zenodo**: records on Zenodo through its public search API
  (``https://zenodo.org/api/records?q=...``).

Only plain ``urllib`` is used, with a User-Agent, timeouts and a size limit on
downloads. Network errors become empty results with a note, never exceptions, so the
model can report them.

Scoring: the query is split into lowercase words (stop words dropped); a candidate
scores 1 for every query word found in its title, description, file names or key
extras, 0.5 for a word found through a chemical alias (``water`` matches ``D2O``),
plus 1 when a multi-word query appears as a whole phrase. Candidates are returned best
first; local candidates with a zero score are dropped, the other sources keep their
top results (the remote search already ranked them).

The Dissolve listing is cached in memory for the process: the unauthenticated GitHub
API allows 60 requests per hour. One Git Trees request lists the whole ``examples``
folder; if GitHub refuses that request (an HTTP error other than the rate limit) the
contents API is used instead, with one request per example folder.
A ``GITHUB_TOKEN`` environment variable, when set, is sent to raise the limit.
"""

from __future__ import annotations

import email.message
import hashlib
import html
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import warnings
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "DataCandidate",
    "SourceError",
    "HttpStatusError",
    "RateLimitError",
    "SOURCES",
    "find_similar_data",
    "search_local",
    "search_dissolve_examples",
    "search_zenodo",
    "download",
    "repository_root",
    "default_local_roots",
    "query_words",
    "score_text",
    "clear_caches",
]

USER_AGENT = "mellor/0.1 (EPSR agent; +https://github.com/)"
DISSOLVE_REPO = "disorderedmaterials/dissolve"
DISSOLVE_BRANCH = "develop"
MAX_DOWNLOAD_BYTES = 50_000_000

GITHUB_API = "https://api.github.com"
RAW_GITHUB = "https://raw.githubusercontent.com"
ZENODO_API = "https://zenodo.org/api/records"

#: The sources in the order used to break score ties.
SOURCES = ("local", "dissolve", "zenodo")
#: Folders of the repository searched by default (relative to :func:`repository_root`).
LOCAL_DATA_FOLDERS = ("experimental_data", "validation", "EPSR-test")
#: Instrument from the prefix of a ``.mint01`` source name (``sls18498`` -> SANDALS).
INSTRUMENTS = {
    "SLS": "SANDALS (ISIS)",
    "NIMROD": "NIMROD (ISIS)",
    "GEM": "GEM (ISIS)",
    "POLARIS": "POLARIS (ISIS)",
    "LAD": "LAD (ISIS)",
    "D4": "D4 (ILL)",
}
#: Query words that also match these (half-weight) chemical names.
ALIASES: dict[str, tuple[str, ...]] = {
    "water": ("h2o", "d2o", "hdo"),
    "silica": ("sio2",),
    "benzene": ("c6h6", "c6d6"),
    "mof": ("cubtc", "cu-btc", "hkust"),
    "hkust": ("cubtc", "cu-btc"),
}
#: Words that carry no information for the search.
STOP_WORDS = frozenset(
    (
        "a", "an", "the", "of", "in", "on", "at", "to", "for", "and", "or", "with",
        "from", "by", "is", "are", "data", "dataset", "datasets", "file", "files",
        "mint01",
    )
)  # fmt: skip
#: ``extra`` entries that take part in the scoring.
SCORED_EXTRA_KEYS = (
    "case",
    "instrument",
    "source",
    "user",
    "keywords",
    "creators",
    "doi",
)
#: What the Dissolve examples contain (species, isotopes, instrument), for the scoring
#: and the description; anything else is described from its file names only.
DISSOLVE_EXAMPLE_HINTS = {
    "water": (
        "Liquid water at 25 C measured on SANDALS (ISIS): H2O, D2O and a 50:50 "
        "H2O/D2O mixture (the EPSR validation water case)."
    ),
    "benzene": "Liquid benzene: C6H6, C6D6 and a 50:50 mixture, neutron total scattering.",
    "silica_bulk": "Bulk amorphous silica (SiO2 glass) measured on NIMROD (ISIS).",
    "cu-btc": (
        "Cu-BTC (HKUST-1) metal-organic framework with empty pores, measured on "
        "NIMROD (ISIS)."
    ),
    "emim-oac": (
        "1-ethyl-3-methylimidazolium acetate ([EMIM][OAc]) ionic liquid measured on "
        "SANDALS (ISIS), several H/D isotopologues of cation and anion."
    ),
    "argon": "Liquid argon (no reference data in the example).",
}

_DOWNLOAD_CHUNK = 64 * 1024
_MAX_HEADER_LINES = 64
_DESCRIPTION_LENGTH = 400
_WORD_SPLIT = re.compile(r"[\s,;]+")
_SOURCE_NAME = re.compile(r"^([A-Za-z]+)0*(\d+)")
_HTML_TAG = re.compile(r"<[^>]+>")
_GITHUB_BLOB = re.compile(r"^https?://github\.com/([^/]+)/([^/]+)/blob/(.+)$")
_CACHE: dict[str, Any] = {}


class SourceError(RuntimeError):
    """A data source could not be searched or a file could not be downloaded."""


class HttpStatusError(SourceError):
    """The server answered with an error status (``code``), e.g. 404."""

    def __init__(self, message: str, code: int):
        super().__init__(message)
        self.code = code


class RateLimitError(HttpStatusError):
    """The remote API refused the request because its rate limit is exhausted."""

    def __init__(self, message: str, code: int = 429):
        super().__init__(message, code)


@dataclass
class DataCandidate:
    """One data set (or record) that may be comparable to the user's."""

    source: str  # "local" | "dissolve" | "zenodo"
    title: str
    description: str
    url: str  # download URL of the file, or the record page
    files: list[str] = field(default_factory=list)  # file names available
    score: float = 0.0  # keyword match score, higher is better
    #: instrument, sample title, DOI, download links...
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """A JSON-serialisable dict with every field."""
        return asdict(self)


# --------------------------------------------------------------------------- scoring
def query_words(query: str) -> list[str]:
    """The distinct lowercase words of a query, without stop words and punctuation."""
    words: list[str] = []
    for raw in _WORD_SPLIT.split(query.lower()):
        word = raw.strip("\"'()[]{}.:!?")
        if len(word) < 2 or word in STOP_WORDS or word in words:
            continue
        words.append(word)
    return words


def _compact(text: str) -> str:
    return text.replace("-", "").replace("_", "")


def score_text(query: str, text: str) -> float:
    """The match score of ``text`` for ``query`` (see the module docstring)."""
    words = query_words(query)
    if not words:
        return 0.0
    haystack = " ".join(text.lower().split())
    compact = _compact(haystack)
    score = 0.0
    for word in words:
        if word in haystack or _compact(word) in compact:
            score += 1.0
        elif any(alias in haystack for alias in ALIASES.get(word, ())):
            score += 0.5
    phrase = " ".join(query.lower().split())
    if len(words) > 1 and phrase in haystack:
        score += 1.0
    return score


def _haystack(candidate: DataCandidate) -> str:
    parts = [candidate.title, candidate.description, *candidate.files]
    for key in SCORED_EXTRA_KEYS:
        value = candidate.extra.get(key)
        if isinstance(value, str):
            parts.append(value)
        elif isinstance(value, (list, tuple)):
            parts.extend(str(item) for item in value)
    return " ".join(parts)


def _rank(
    candidates: list[DataCandidate], query: str, *, limit: int, keep_zero: bool
) -> list[DataCandidate]:
    for candidate in candidates:
        candidate.score = score_text(query, _haystack(candidate))
    ranked = sorted(candidates, key=lambda candidate: -candidate.score)
    if not keep_zero:
        ranked = [candidate for candidate in ranked if candidate.score > 0]
    return ranked[: max(0, limit)]


def _check_query(query: str) -> str:
    if not isinstance(query, str) or not query.strip():
        raise ValueError("the search query must not be empty")
    return query.strip()


# --------------------------------------------------------------------------- local
def repository_root() -> Path:
    """The EPSR_hackathon repository: the parent of ``mellor/`` that holds the data
    folders (computed from this file, so it also works for an editable install)."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "mellor").is_dir() and any(
            (parent / name).is_dir() for name in LOCAL_DATA_FOLDERS
        ):
            return parent
    return here.parents[3] if len(here.parents) > 3 else here.parent


def default_local_roots() -> list[Path]:
    """The repository data folders that exist."""
    root = repository_root()
    return [root / name for name in LOCAL_DATA_FOLDERS if (root / name).is_dir()]


def _read_header(path: Path) -> dict[str, str]:
    """``source``, ``title``, ``user`` and ``date`` from the ``#`` header of a ``.mint01``."""
    lines: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("#") or len(lines) >= _MAX_HEADER_LINES:
                break
            lines.append(line[1:].strip())
    header = {"source": "", "title": "", "user": "", "date": ""}
    if lines:
        header["source"] = lines[0]
    if len(lines) > 1:
        header["title"] = lines[1]
    if len(lines) > 4 and not re.fullmatch(r"[\d\s.eE+-]*", lines[4]):
        header["user"] = lines[4]
    if len(lines) > 5 and re.search(r"\d{4}", lines[5]):
        header["date"] = lines[5]
    return header


def _instrument(source: str) -> tuple[str, str]:
    """Instrument name and run number from a source name such as ``sls18498.mdcs01``."""
    match = _SOURCE_NAME.match(source.strip())
    if not match:
        return "", ""
    prefix, run = match.group(1).upper(), match.group(2)
    return INSTRUMENTS.get(prefix, ""), run


def _is_hidden(path: Path, root: Path) -> bool:
    try:
        parts = path.relative_to(root).parts
    except ValueError:
        parts = path.parts
    return any(part.startswith(".") for part in parts)


def _relative_path(path: Path, root: Path) -> str:
    for base in (repository_root(), root):
        try:
            return path.relative_to(base).as_posix()
        except ValueError:
            continue
    return path.as_posix()


def _local_candidate(path: Path, root: Path) -> DataCandidate:
    header = _read_header(path)
    instrument, run = _instrument(header["source"])
    relative = _relative_path(path, root)
    parts = [
        instrument + (f" run {run}" if run else "") if instrument else "",
        relative,
    ]
    if header["user"]:
        parts.append(f"user {header['user']}")
    if header["date"]:
        parts.append(header["date"])
    return DataCandidate(
        source="local",
        title=header["title"] or path.stem,
        description="; ".join(part for part in parts if part),
        url=str(path),
        files=[path.name],
        extra={
            "path": str(path),
            "folder": str(path.parent),
            "source": header["source"],
            "instrument": instrument,
            "run_number": run,
            "user": header["user"],
            "date": header["date"],
            "duplicates": [],
        },
    )


def search_local(
    query: str, roots: Sequence[str | Path] | None = None, *, limit: int = 10
) -> list[DataCandidate]:
    """``.mint01`` files under ``roots`` (default: the repository data folders) whose
    header title or file name matches the query words.

    The folder path, the instrument (from the source name) and the user are searched
    too, so ``water SANDALS`` finds ``experimental_data/water/data/*.mint01``. Files
    with identical content are reported once, the other copies in
    ``extra["duplicates"]``. Named roots must exist (``FileNotFoundError``); files
    that cannot be read are skipped with a warning.
    """
    query = _check_query(query)
    if roots is None:
        search_roots = default_local_roots()
    else:
        search_roots = [Path(root).expanduser() for root in roots]
        for root in search_roots:
            if not root.is_dir():
                raise FileNotFoundError(f"folder not found: {root}")
    seen: dict[str, DataCandidate] = {}
    for root in search_roots:
        root = root.resolve()
        for path in sorted(root.rglob("*")):
            if path.suffix.lower() != ".mint01" or not path.is_file():
                continue
            if _is_hidden(path, root):
                continue
            try:
                digest = hashlib.sha1(path.read_bytes()).hexdigest()
                candidate = _local_candidate(path, root)
            except OSError as error:
                warnings.warn(f"skipping {path}: {error}", stacklevel=2)
                continue
            if digest in seen:
                seen[digest].extra["duplicates"].append(str(path))
            else:
                seen[digest] = candidate
    return _rank(list(seen.values()), query, limit=limit, keep_zero=False)


# --------------------------------------------------------------------------- HTTP
def _request(url: str, headers: dict[str, str] | None = None) -> urllib.request.Request:
    merged = {"User-Agent": USER_AGENT}
    merged.update(headers or {})
    return urllib.request.Request(url, headers=merged)


def _rate_limit_error(
    service: str, error: urllib.error.HTTPError
) -> RateLimitError | None:
    if error.code not in (403, 429):
        return None
    remaining = error.headers.get("X-RateLimit-Remaining") if error.headers else None
    if remaining is None or remaining.strip() != "0":
        return None
    message = f"{service} rate limit exceeded"
    reset = error.headers.get("X-RateLimit-Reset")
    if reset and reset.strip().isdigit():
        minutes = max(1, round((int(reset) - time.time()) / 60))
        message += f"; it resets in about {minutes} min"
    if service == "GitHub":
        message += (
            " (the GitHub API allows 60 requests per hour without a GITHUB_TOKEN)"
        )
    return RateLimitError(message, error.code)


def _get_json(
    url: str, *, service: str, timeout: float, headers: dict[str, str] | None = None
) -> Any:
    """GET ``url`` and parse the JSON body; every failure becomes a :class:`SourceError`."""
    request = _request(url, {"Accept": "application/json", **(headers or {})})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read()
    except urllib.error.HTTPError as error:
        limited = _rate_limit_error(service, error)
        if limited is not None:
            raise limited from error
        raise HttpStatusError(
            f"{service} request failed: HTTP {error.code} {error.reason} ({url})",
            error.code,
        ) from error
    except urllib.error.URLError as error:
        raise SourceError(f"{service} unreachable: {error.reason}") from error
    except (OSError, TimeoutError) as error:
        raise SourceError(f"{service} request failed: {error}") from error
    try:
        return json.loads(payload.decode("utf-8"))
    except ValueError as error:
        raise SourceError(f"{service} returned invalid JSON ({url})") from error


# --------------------------------------------------------------------------- dissolve
def _github_headers() -> dict[str, str]:
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    token = os.environ.get("GITHUB_TOKEN", "").strip()
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers


def _github_json(url: str, timeout: float) -> Any:
    return _get_json(url, service="GitHub", timeout=timeout, headers=_github_headers())


def _dissolve_listing_from_tree(timeout: float) -> dict[str, list[str]]:
    """``case -> [paths inside the example]`` from one Git Trees request."""
    url = (
        f"{GITHUB_API}/repos/{DISSOLVE_REPO}/git/trees/"
        f"{DISSOLVE_BRANCH}:examples?recursive=1"
    )
    tree = _github_json(url, timeout)
    listing: dict[str, list[str]] = {}
    for entry in tree.get("tree", []) if isinstance(tree, dict) else []:
        if entry.get("type") != "blob":
            continue
        case, separator, rest = str(entry.get("path", "")).partition("/")
        if separator and rest:
            listing.setdefault(case, []).append(rest)
    if not listing:
        raise SourceError("GitHub returned no files for the Dissolve examples tree")
    return listing


def _dissolve_listing_from_contents(timeout: float) -> dict[str, list[str]]:
    """The same listing through the contents API: one request per example folder."""
    base = f"{GITHUB_API}/repos/{DISSOLVE_REPO}/contents/examples"
    entries = _github_json(f"{base}?ref={DISSOLVE_BRANCH}", timeout)
    listing: dict[str, list[str]] = {}
    for entry in entries if isinstance(entries, list) else []:
        if entry.get("type") != "dir":
            continue
        case = str(entry.get("name", ""))
        listing[case] = ["input.txt"]
        try:
            files = _github_json(f"{base}/{case}/data?ref={DISSOLVE_BRANCH}", timeout)
        except RateLimitError:
            raise
        except HttpStatusError:  # no data folder (404)
            continue
        for item in files if isinstance(files, list) else []:
            if item.get("type") == "file":
                listing[case].append(f"data/{item.get('name', '')}")
    if not listing:
        raise SourceError("GitHub returned no example folders for Dissolve")
    return listing


def _dissolve_listing(timeout: float) -> dict[str, list[str]]:
    """The cached ``case -> [paths]`` listing of the Dissolve examples."""
    listing = _CACHE.get("dissolve")
    if listing is None:
        try:
            listing = _dissolve_listing_from_tree(timeout)
        except RateLimitError:
            raise
        except HttpStatusError:  # the trees endpoint refused: use the contents API
            listing = _dissolve_listing_from_contents(timeout)
        _CACHE["dissolve"] = listing
    return listing


def _dissolve_candidate(case: str, paths: list[str]) -> DataCandidate | None:
    raw = f"{RAW_GITHUB}/{DISSOLVE_REPO}/{DISSOLVE_BRANCH}/examples/{case}"
    data_files = sorted(
        path for path in paths if path.lower().endswith(".mint01") and "/" in path
    )
    if not data_files:
        return None
    names = [path.rsplit("/", 1)[-1] for path in data_files]
    hint = DISSOLVE_EXAMPLE_HINTS.get(case, "")
    summary = (
        f"Dissolve example '{case}' with {len(names)} .mint01 data file(s): "
        + ", ".join(names)
    )
    return DataCandidate(
        source="dissolve",
        title=f"Dissolve example: {case}",
        description=(hint + " " if hint else "") + summary,
        url=f"https://github.com/{DISSOLVE_REPO}/tree/{DISSOLVE_BRANCH}/examples/{case}",
        files=names,
        extra={
            "case": case,
            "downloads": {
                name: f"{raw}/{urllib.parse.quote(path)}"
                for name, path in zip(names, data_files, strict=True)
            },
            "input_file": f"{raw}/input.txt",
            "repository": DISSOLVE_REPO,
        },
    )


def search_dissolve_examples(
    query: str, *, limit: int = 10, timeout: float = 20.0
) -> list[DataCandidate]:
    """Dissolve example data sets matching the query words (case name, species, data names).

    One candidate per example folder that has ``data/*.mint01`` files; the raw
    download URL of every file is in ``extra["downloads"]`` and the example's
    ``input.txt`` in ``extra["input_file"]``. Raises :class:`SourceError` (or
    :class:`RateLimitError`) when GitHub cannot be queried.
    """
    query = _check_query(query)
    listing = _dissolve_listing(timeout)
    candidates = []
    for case in sorted(listing):
        candidate = _dissolve_candidate(case, listing[case])
        if candidate is not None:
            candidates.append(candidate)
    return _rank(candidates, query, limit=limit, keep_zero=True)


# --------------------------------------------------------------------------- zenodo
def _clean_html(text: str, length: int = _DESCRIPTION_LENGTH) -> str:
    plain = html.unescape(_HTML_TAG.sub(" ", text or ""))
    plain = " ".join(plain.split())
    if len(plain) > length:
        plain = plain[: length - 1].rstrip() + "…"
    return plain


def _zenodo_candidate(hit: dict[str, Any]) -> DataCandidate:
    metadata = hit.get("metadata") or {}
    links = hit.get("links") or {}
    record_id = hit.get("id") or hit.get("recid") or ""
    title = metadata.get("title") or hit.get("title") or f"Zenodo record {record_id}"
    files = [item for item in hit.get("files") or [] if isinstance(item, dict)]
    names = [str(item.get("key", "")) for item in files if item.get("key")]
    downloads = {
        str(item["key"]): str((item.get("links") or {}).get("self", ""))
        for item in files
        if item.get("key")
    }
    creators = [
        str(creator.get("name", ""))
        for creator in metadata.get("creators") or []
        if isinstance(creator, dict) and creator.get("name")
    ]
    resource = metadata.get("resource_type") or {}
    keywords = metadata.get("keywords") or []
    return DataCandidate(
        source="zenodo",
        title=str(title),
        description=_clean_html(str(metadata.get("description") or "")),
        url=str(
            links.get("self_html")
            or hit.get("doi_url")
            or (f"https://zenodo.org/records/{record_id}" if record_id else "")
        ),
        files=names,
        extra={
            "doi": str(metadata.get("doi") or hit.get("doi") or ""),
            "record_id": str(record_id),
            "creators": creators,
            "publication_date": str(metadata.get("publication_date") or ""),
            "resource_type": str(resource.get("title") or resource.get("type") or ""),
            "access": str(metadata.get("access_right") or ""),
            "keywords": [
                str(keyword) for keyword in keywords if isinstance(keyword, str)
            ],
            "downloads": downloads,
        },
    )


def search_zenodo(
    query: str, *, limit: int = 10, timeout: float = 20.0
) -> list[DataCandidate]:
    """Zenodo records for the query (title, DOI, file names, record URL).

    The record page is ``url``; the download link of every file is in
    ``extra["downloads"]``. Raises :class:`SourceError` when Zenodo cannot be queried.
    """
    query = _check_query(query)
    parameters = {"q": query, "size": max(1, min(int(limit), 100)), "sort": "bestmatch"}
    url = (
        ZENODO_API
        + "?"
        + urllib.parse.urlencode(parameters, quote_via=urllib.parse.quote)
    )
    payload = _get_json(url, service="Zenodo", timeout=timeout)
    hits = (
        (payload.get("hits") or {}).get("hits") if isinstance(payload, dict) else None
    )
    if not isinstance(hits, list):
        raise SourceError("Zenodo returned an unexpected response (no hits)")
    candidates = [_zenodo_candidate(hit) for hit in hits if isinstance(hit, dict)]
    return _rank(candidates, query, limit=limit, keep_zero=True)


# --------------------------------------------------------------------------- combined
def find_similar_data(
    query: str,
    *,
    sources: Sequence[str] | str | None = ("local", "dissolve", "zenodo"),
    limit: int = 10,
) -> dict[str, Any]:
    """``{"query": ..., "candidates": [DataCandidate dicts, best first], "notes": [...]}``.

    At most ``limit`` candidates per source, sorted by score and then by source
    (local, dissolve, zenodo). A source that fails (network, rate limit, unknown
    name) adds a note such as ``"zenodo: Zenodo unreachable: ..."`` instead of
    raising; ``sources`` may also be a comma-separated string, and ``None`` means
    all of them.
    """
    query = _check_query(query)
    if sources is None:
        sources = SOURCES
    elif isinstance(sources, str):
        sources = tuple(part for part in sources.split(",") if part.strip())
    searchers = {
        "local": lambda: search_local(query, limit=limit),
        "dissolve": lambda: search_dissolve_examples(query, limit=limit),
        "zenodo": lambda: search_zenodo(query, limit=limit),
    }
    notes: list[str] = []
    candidates: list[DataCandidate] = []
    for name in sources:
        source = str(name).strip().lower()
        searcher = searchers.get(source)
        if searcher is None:
            notes.append(f"{source}: unknown source (choose from {', '.join(SOURCES)})")
            continue
        try:
            found = searcher()
        except Exception as error:  # every failure must become a note, not a crash
            notes.append(f"{source}: {error}")
            continue
        if not found:
            notes.append(f"{source}: no matches")
        candidates.extend(found)
    order = {name: index for index, name in enumerate(SOURCES)}
    candidates.sort(
        key=lambda candidate: (-candidate.score, order.get(candidate.source, 99))
    )
    return {
        "query": query,
        "candidates": [candidate.to_dict() for candidate in candidates],
        "notes": notes,
    }


# --------------------------------------------------------------------------- download
def _download_url(url: str) -> str:
    """Validate the URL; GitHub ``blob`` page links become raw file links."""
    if not isinstance(url, str) or not url.strip():
        raise ValueError("the download URL must not be empty")
    url = url.strip()
    match = _GITHUB_BLOB.match(url)
    if match:
        owner, repo, rest = match.groups()
        url = f"{RAW_GITHUB}/{owner}/{repo}/{rest}"
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme not in ("http", "https") or not parsed.netloc:
        raise ValueError(f"only http(s) URLs can be downloaded, not {url!r}")
    return url


def _download_name(url: str, content_disposition: str | None) -> str:
    """The file name: from ``Content-Disposition`` if given, else from the URL path."""
    name = ""
    if content_disposition:
        message = email.message.Message()
        message["Content-Disposition"] = content_disposition
        name = message.get_filename() or ""
    if not name:
        segments = [
            urllib.parse.unquote(segment)
            for segment in urllib.parse.urlparse(url).path.split("/")
            if segment
        ]
        # Zenodo file links end in ".../files/<name>/content".
        if len(segments) > 1 and segments[-1] in ("content", "download"):
            segments.pop()
        name = segments[-1] if segments else ""
    name = Path(name.replace("\\", "/")).name.strip()
    if not name or name in (".", ".."):
        name = "download"
    return name


def download(url: str, folder: str | Path, *, timeout: float = 60.0) -> Path:
    """Download an ``http(s)`` URL into ``folder`` (file name from the URL) and return the path.

    Refuses other schemes and files over :data:`MAX_DOWNLOAD_BYTES`.

    The folder is created if needed and an existing file of the same name is
    replaced. The body is streamed through a ``.part`` file so a failed download
    leaves nothing behind. A web page (``text/html``) is refused with a hint, since
    record pages are not data files. ``ValueError`` for a bad URL,
    :class:`SourceError` for everything the network or the server does wrong.
    """
    url = _download_url(url)
    folder = Path(folder).expanduser()
    folder.mkdir(parents=True, exist_ok=True)
    request = _request(url)
    temp: Path | None = None
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            headers = response.headers
            length = (headers.get("Content-Length") or "").strip()
            if length.isdigit() and int(length) > MAX_DOWNLOAD_BYTES:
                raise SourceError(
                    f"refusing to download {url}: {int(length)} bytes is over the "
                    f"limit of {MAX_DOWNLOAD_BYTES} bytes"
                )
            name = _download_name(url, headers.get("Content-Disposition"))
            content_type = (
                (headers.get("Content-Type") or "").split(";")[0].strip().lower()
            )
            if content_type == "text/html" and not name.lower().endswith(
                (".html", ".htm")
            ):
                raise SourceError(
                    f"{url} is a web page, not a data file; use the file's download "
                    "link (for example one of extra['downloads'] of a candidate)"
                )
            target = folder / name
            temp = target.with_name(target.name + ".part")
            total = 0
            with open(temp, "wb") as handle:
                while True:
                    chunk = response.read(_DOWNLOAD_CHUNK)
                    if not chunk:
                        break
                    total += len(chunk)
                    if total > MAX_DOWNLOAD_BYTES:
                        raise SourceError(
                            f"download of {url} aborted: over the limit of "
                            f"{MAX_DOWNLOAD_BYTES} bytes"
                        )
                    handle.write(chunk)
            temp.replace(target)
            temp = None
    except urllib.error.HTTPError as error:
        limited = _rate_limit_error("the server", error)
        if limited is not None:
            raise limited from error
        raise HttpStatusError(
            f"download failed: HTTP {error.code} {error.reason} ({url})", error.code
        ) from error
    except urllib.error.URLError as error:
        raise SourceError(f"download failed: {error.reason} ({url})") from error
    except (OSError, TimeoutError) as error:
        raise SourceError(f"download of {url} failed: {error}") from error
    finally:
        if temp is not None and temp.exists():
            temp.unlink()
    return target


def clear_caches() -> None:
    """Forget the cached remote listings (tests, or to see new Dissolve examples)."""
    _CACHE.clear()
