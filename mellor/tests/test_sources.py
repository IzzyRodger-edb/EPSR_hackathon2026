"""Tests for ``mellor.sources``.

The offline tests replace ``urllib.request.urlopen`` with fake responses; the
``network`` tests hit GitHub, Zenodo and raw.githubusercontent.com for real.
"""

from __future__ import annotations

import email.message
import io
import json
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from mellor import sources
from mellor.sources import (
    DISSOLVE_REPO,
    GITHUB_API,
    RAW_GITHUB,
    DataCandidate,
    RateLimitError,
    SourceError,
    download,
    find_similar_data,
    query_words,
    score_text,
    search_dissolve_examples,
    search_local,
    search_zenodo,
)

REPO = sources.repository_root()
WATER_FILES = ["SLS18498-H2O.mint01", "SLS18500-HDO5050.mint01", "SLS18502-D2O.mint01"]
TREE_URL = f"{GITHUB_API}/repos/{DISSOLVE_REPO}/git/trees/develop:examples?recursive=1"
CONTENTS_URL = f"{GITHUB_API}/repos/{DISSOLVE_REPO}/contents/examples"
RAW_WATER = f"{RAW_GITHUB}/{DISSOLVE_REPO}/develop/examples/water/data"

TREE = {
    "sha": "abc",
    "truncated": False,
    "tree": [
        {"path": "CMakeLists.txt", "type": "blob"},
        {"path": "water", "type": "tree"},
        {"path": "water/input.txt", "type": "blob"},
        {"path": "water/data", "type": "tree"},
        {"path": "water/data/SLS18502-D2O.mint01", "type": "blob"},
        {"path": "water/data/SLS18498-H2O.mint01", "type": "blob"},
        {"path": "water/data/SLS18500-HDO5050.mint01", "type": "blob"},
        {"path": "benzene/input.txt", "type": "blob"},
        {"path": "benzene/data/C6H6.mint01", "type": "blob"},
        {"path": "benzene/data/C6D6.mint01", "type": "blob"},
        {"path": "benzene/data/5050.mint01", "type": "blob"},
        {"path": "argon/input.txt", "type": "blob"},
    ],
}

ZENODO_HITS = {
    "hits": {
        "total": 2,
        "hits": [
            {
                "id": 123,
                "metadata": {
                    "title": "Neutron diffraction of liquid water",
                    "description": "<p>Total scattering &amp; EPSR of <b>water</b> at 298 K.</p>",
                    "doi": "10.5281/zenodo.123",
                    "creators": [{"name": "Doe, Jane"}],
                    "publication_date": "2024-01-01",
                    "resource_type": {"title": "Dataset", "type": "dataset"},
                    "keywords": ["water", "SANDALS"],
                },
                "links": {"self_html": "https://zenodo.org/records/123"},
                "files": [
                    {
                        "key": "H2O.mint01",
                        "links": {
                            "self": "https://zenodo.org/api/records/123/files/H2O.mint01/content"
                        },
                    }
                ],
            },
            {
                "id": 456,
                "metadata": {
                    "title": "Rietveld refinement of a zeolite",
                    "description": "",
                    "doi": "10.5281/zenodo.456",
                    "creators": [],
                },
                "links": {"self_html": "https://zenodo.org/records/456"},
            },
        ],
    }
}


# --------------------------------------------------------------------------- fakes
class FakeResponse:
    """What ``urlopen`` returns: a context manager with headers and a streaming read."""

    def __init__(self, body: bytes, headers: dict[str, str] | None = None):
        self._buffer = io.BytesIO(body)
        self.headers = email.message.Message()
        for key, value in (headers or {}).items():
            self.headers[key] = value
        self.status = 200

    def read(self, size: int = -1) -> bytes:
        return self._buffer.read(size)

    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *exc: object) -> bool:
        return False


def json_response(payload: Any) -> FakeResponse:
    body = json.dumps(payload).encode()
    return FakeResponse(
        body, {"Content-Type": "application/json", "Content-Length": str(len(body))}
    )


def http_error(
    url: str, code: int, reason: str, headers: dict[str, str] | None = None
) -> urllib.error.HTTPError:
    message = email.message.Message()
    for key, value in (headers or {}).items():
        message[key] = value
    return urllib.error.HTTPError(url, code, reason, message, io.BytesIO(b""))


class FakeNetwork:
    """Routes exact URLs to response factories (or exceptions) and records requests."""

    def __init__(self) -> None:
        self.routes: dict[str, Callable[[], Any]] = {}
        self.calls: list[urllib.request.Request] = []

    def add(self, url: str, factory: Callable[[], Any]) -> None:
        self.routes[url] = factory

    def fail(self, url: str, error: BaseException) -> None:
        def raise_it() -> Any:
            raise error

        self.routes[url] = raise_it

    def urlopen(self, request: Any, timeout: float | None = None) -> Any:
        assert isinstance(request, urllib.request.Request)
        assert request.get_header("User-agent") == sources.USER_AGENT
        assert timeout is not None and timeout > 0
        self.calls.append(request)
        factory = self.routes.get(request.full_url)
        if factory is None:
            raise AssertionError(f"unexpected request: {request.full_url}")
        return factory()

    @property
    def urls(self) -> list[str]:
        return [request.full_url for request in self.calls]


@pytest.fixture
def network(monkeypatch: pytest.MonkeyPatch) -> FakeNetwork:
    fake = FakeNetwork()
    monkeypatch.setattr(urllib.request, "urlopen", fake.urlopen)
    sources.clear_caches()
    yield fake
    sources.clear_caches()


def write_mint01(path: Path, header: list[str], rows: int = 5) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {line}" for line in header]
    lines += [f"  {0.5 + 0.05 * i:.6E}  {-0.3:.6E}  {0.01:.6E}" for i in range(rows)]
    path.write_text("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------- scoring
def test_query_words_drops_stop_words_and_punctuation() -> None:
    assert query_words("Water SANDALS, the data (H2O).") == ["water", "sandals", "h2o"]
    assert query_words("in of the") == []
    assert query_words("") == []
    assert query_words("water water") == ["water"]


def test_score_text_counts_words_aliases_and_phrase() -> None:
    assert score_text("water sandals", "H2O in a can; SANDALS (ISIS)") == 1.5
    assert score_text("cu-btc", "CuBTC pellets") == 1.0
    assert score_text("cu btc", "experimental_data/cu-btc/data") == 2.0
    assert score_text("H2O in 1mm can", "H2O in 1mm can 5 in position 5 : 25C") == 4.0
    assert score_text("argon", "water") == 0.0
    assert score_text("", "anything") == 0.0


def test_candidate_to_dict_is_json_serialisable() -> None:
    candidate = DataCandidate(
        source="zenodo",
        title="t",
        description="d",
        url="https://zenodo.org/records/1",
        files=["a.mint01"],
        score=1.5,
        extra={"doi": "10.1/x", "downloads": {"a.mint01": "https://x/a"}},
    )
    text = json.dumps(candidate.to_dict())
    assert json.loads(text)["extra"]["downloads"]["a.mint01"] == "https://x/a"
    assert set(json.loads(text)) == {
        "source", "title", "description", "url", "files", "score", "extra"
    }  # fmt: skip


# --------------------------------------------------------------------------- local
def _needs_repository_data() -> None:
    if not (REPO / "experimental_data" / "water" / "data").is_dir():
        pytest.skip("repository data folders not found")


def test_search_local_water_sandals_finds_the_three_water_files() -> None:
    _needs_repository_data()
    results = search_local("water SANDALS")
    names = [candidate.files[0] for candidate in results]
    assert names[:3] == WATER_FILES
    for candidate in results[:3]:
        assert candidate.source == "local"
        assert candidate.score == 2.0
        assert Path(candidate.url).is_file()
        assert candidate.extra["instrument"] == "SANDALS (ISIS)"
        assert candidate.extra["user"] == "Bychkov"
        # The identical copy in validation/water is reported as a duplicate.
        assert len(candidate.extra["duplicates"]) == 1
        assert "validation/water" in candidate.extra["duplicates"][0].replace("\\", "/")
    first = results[0]
    assert first.title == "H2O in 1mm can 5 in position 5 : 25C"
    assert first.extra["run_number"] == "18498"
    assert "experimental_data/water/data/SLS18498-H2O.mint01" in first.description
    assert "SiO2_NIMROD.mint01" not in names
    assert "Empty-CuBTC.mint01" not in names
    assert [c.score for c in results] == sorted(
        (c.score for c in results), reverse=True
    )
    json.dumps([candidate.to_dict() for candidate in results])


def test_search_local_silica_finds_sio2_nimrod() -> None:
    _needs_repository_data()
    results = search_local("silica")
    assert len(results) == 1
    (candidate,) = results
    assert candidate.files == ["SiO2_NIMROD.mint01"]
    assert candidate.score == 1.0
    assert candidate.title == "3mm SiO2 flat plate in pos 9"
    assert candidate.extra["instrument"] == "NIMROD (ISIS)"
    assert candidate.extra["run_number"] == "1122"
    assert len(candidate.extra["duplicates"]) == 1


def test_search_local_limit_roots_and_errors() -> None:
    _needs_repository_data()
    assert len(search_local("water SANDALS", limit=2)) == 2
    results = search_local("water", roots=[REPO / "validation" / "water"])
    assert sorted(candidate.files[0] for candidate in results) == WATER_FILES
    assert all(candidate.extra["duplicates"] == [] for candidate in results)
    assert all("validation/water" in candidate.description for candidate in results)
    assert search_local("argon") == []
    with pytest.raises(FileNotFoundError):
        search_local("water", roots=[REPO / "no-such-folder"])
    with pytest.raises(ValueError):
        search_local("   ")


def test_search_local_custom_folder_reads_headers(tmp_path: Path) -> None:
    header = [
        "lad12345.raw",
        " Benzene C6D6 in a flat cell",
        "    1   100",
        "   10",
        "Someone",
        "01-JAN-2020 10:00:00",
    ]
    sample = write_mint01(tmp_path / "sample.mint01", header)
    copy = write_mint01(tmp_path / "sub" / "copy.mint01", header)
    write_mint01(tmp_path / ".hidden" / "secret.mint01", ["x", " Benzene hidden"])
    bare = tmp_path / "bare.MINT01"
    bare.write_text("  0.5  -0.3  0.01\n")

    results = search_local("benzene", roots=[str(tmp_path)])
    assert len(results) == 1
    (candidate,) = results
    assert candidate.url == str(sample)
    assert candidate.title == "Benzene C6D6 in a flat cell"
    assert candidate.extra["source"] == "lad12345.raw"
    assert candidate.extra["instrument"] == "LAD (ISIS)"
    assert candidate.extra["run_number"] == "12345"
    assert candidate.extra["user"] == "Someone"
    assert candidate.extra["date"] == "01-JAN-2020 10:00:00"
    assert candidate.extra["duplicates"] == [str(copy)]

    (candidate,) = search_local("bare", roots=[tmp_path])
    assert candidate.title == "bare"
    assert candidate.extra["source"] == ""
    assert candidate.extra["instrument"] == ""


# --------------------------------------------------------------------------- dissolve
def test_search_dissolve_examples_offline_and_cached(network: FakeNetwork) -> None:
    network.add(TREE_URL, lambda: json_response(TREE))
    results = search_dissolve_examples("water")
    assert [candidate.extra["case"] for candidate in results] == ["water", "benzene"]
    water, benzene = results
    assert water.source == "dissolve"
    assert water.score >= 1.0
    assert benzene.score == 0.0  # kept: the remote source keeps its top results
    assert water.files == WATER_FILES
    assert (
        water.url == f"https://github.com/{DISSOLVE_REPO}/tree/develop/examples/water"
    )
    assert water.extra["downloads"] == {
        name: f"{RAW_WATER}/{name}" for name in WATER_FILES
    }
    assert water.extra["input_file"] == (
        f"{RAW_GITHUB}/{DISSOLVE_REPO}/develop/examples/water/input.txt"
    )
    assert "SANDALS" in water.description
    assert benzene.files == ["5050.mint01", "C6D6.mint01", "C6H6.mint01"]
    request = network.calls[0]
    assert request.get_header("Accept") == "application/vnd.github+json"

    # A second search uses the cached listing: no new request.
    results = search_dissolve_examples("benzene C6D6", limit=1)
    assert [candidate.extra["case"] for candidate in results] == ["benzene"]
    assert results[0].score == 2.0
    assert len(network.calls) == 1


def test_search_dissolve_rate_limit_is_reported(network: FakeNetwork) -> None:
    reset = str(int(time.time()) + 600)
    network.fail(
        TREE_URL,
        http_error(
            TREE_URL,
            403,
            "rate limit exceeded",
            {"X-RateLimit-Remaining": "0", "X-RateLimit-Reset": reset},
        ),
    )
    with pytest.raises(RateLimitError) as info:
        search_dissolve_examples("water")
    assert "rate limit exceeded" in str(info.value)
    assert "10 min" in str(info.value)
    assert len(network.calls) == 1  # no fallback to the contents API

    result = find_similar_data("water", sources=("dissolve",))
    assert result["candidates"] == []
    assert result["notes"][0].startswith("dissolve: GitHub rate limit exceeded")
    assert len(network.calls) == 2  # failures are not cached


def test_search_dissolve_falls_back_to_contents_api(network: FakeNetwork) -> None:
    network.fail(TREE_URL, http_error(TREE_URL, 404, "Not Found"))
    network.add(
        f"{CONTENTS_URL}?ref=develop",
        lambda: json_response(
            [
                {"name": "CMakeLists.txt", "type": "file"},
                {"name": "argon", "type": "dir"},
                {"name": "water", "type": "dir"},
            ]
        ),
    )
    network.fail(
        f"{CONTENTS_URL}/argon/data?ref=develop",
        http_error(f"{CONTENTS_URL}/argon/data?ref=develop", 404, "Not Found"),
    )
    network.add(
        f"{CONTENTS_URL}/water/data?ref=develop",
        lambda: json_response([{"name": "SLS18498-H2O.mint01", "type": "file"}]),
    )
    results = search_dissolve_examples("water")
    assert [candidate.extra["case"] for candidate in results] == ["water"]
    assert results[0].files == ["SLS18498-H2O.mint01"]
    assert len(network.calls) == 4
    search_dissolve_examples("argon")
    assert len(network.calls) == 4


def test_search_dissolve_unreachable(network: FakeNetwork) -> None:
    network.fail(TREE_URL, urllib.error.URLError("Name or service not known"))
    with pytest.raises(SourceError, match="GitHub unreachable"):
        search_dissolve_examples("water")
    assert len(network.calls) == 1  # no point in trying the contents API


# --------------------------------------------------------------------------- zenodo
def test_search_zenodo_offline(network: FakeNetwork) -> None:
    url = f"{sources.ZENODO_API}?q=water%20neutron&size=5&sort=bestmatch"
    network.add(url, lambda: json_response(ZENODO_HITS))
    results = search_zenodo("water neutron", limit=5)
    assert [candidate.title for candidate in results] == [
        "Neutron diffraction of liquid water",
        "Rietveld refinement of a zeolite",
    ]
    water, zeolite = results
    assert water.source == "zenodo"
    assert water.score == 2.0
    assert water.description == "Total scattering & EPSR of water at 298 K."
    assert water.url == "https://zenodo.org/records/123"
    assert water.files == ["H2O.mint01"]
    assert water.extra["downloads"] == {
        "H2O.mint01": "https://zenodo.org/api/records/123/files/H2O.mint01/content"
    }
    assert water.extra["doi"] == "10.5281/zenodo.123"
    assert water.extra["creators"] == ["Doe, Jane"]
    assert water.extra["keywords"] == ["water", "SANDALS"]
    assert water.extra["resource_type"] == "Dataset"
    assert zeolite.score == 0.0
    assert zeolite.files == []
    assert zeolite.extra["downloads"] == {}
    json.dumps([candidate.to_dict() for candidate in results])


def test_search_zenodo_errors(network: FakeNetwork) -> None:
    url = f"{sources.ZENODO_API}?q=water&size=10&sort=bestmatch"
    network.fail(url, urllib.error.URLError("timed out"))
    with pytest.raises(SourceError, match="Zenodo unreachable"):
        search_zenodo("water")
    network.add(url, lambda: FakeResponse(b"<html>not json</html>"))
    with pytest.raises(SourceError, match="invalid JSON"):
        search_zenodo("water")
    network.add(url, lambda: json_response({"message": "no hits here"}))
    with pytest.raises(SourceError, match="unexpected response"):
        search_zenodo("water")
    network.fail(url, http_error(url, 500, "Internal Server Error"))
    with pytest.raises(SourceError, match="HTTP 500"):
        search_zenodo("water")
    with pytest.raises(ValueError):
        search_zenodo("")


# --------------------------------------------------------------------------- combined
def test_find_similar_data_merges_sources_and_notes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def make(source: str, title: str, score: float) -> DataCandidate:
        return DataCandidate(source, title, "", f"https://x/{title}", score=score)

    def fake_local(
        query: str, roots: Any = None, *, limit: int = 10
    ) -> list[DataCandidate]:
        return [make("local", "local-2", 2.0), make("local", "local-1", 1.0)][:limit]

    def fake_zenodo(
        query: str, *, limit: int = 10, timeout: float = 20.0
    ) -> list[DataCandidate]:
        return [make("zenodo", "zenodo-3", 3.0), make("zenodo", "zenodo-2", 2.0)][
            :limit
        ]

    def fake_dissolve(
        query: str, *, limit: int = 10, timeout: float = 20.0
    ) -> list[DataCandidate]:
        raise SourceError("GitHub unreachable: boom")

    monkeypatch.setattr(sources, "search_local", fake_local)
    monkeypatch.setattr(sources, "search_zenodo", fake_zenodo)
    monkeypatch.setattr(sources, "search_dissolve_examples", fake_dissolve)

    result = find_similar_data("water SANDALS")
    assert result["query"] == "water SANDALS"
    # Sorted by score, ties broken by the source order local, dissolve, zenodo.
    assert [c["title"] for c in result["candidates"]] == [
        "zenodo-3", "local-2", "zenodo-2", "local-1"
    ]  # fmt: skip
    assert result["notes"] == ["dissolve: GitHub unreachable: boom"]
    json.dumps(result)

    result = find_similar_data("water", sources="local, nope", limit=1)
    assert [c["title"] for c in result["candidates"]] == ["local-2"]
    assert result["notes"] == [
        "nope: unknown source (choose from local, dissolve, zenodo)"
    ]

    # None means every source; a list (JSON array from the model) works too.
    result = find_similar_data("water", sources=None, limit=1)
    assert [c["source"] for c in result["candidates"]] == ["zenodo", "local"]
    result = find_similar_data("water", sources=["zenodo"], limit=1)
    assert [c["title"] for c in result["candidates"]] == ["zenodo-3"]
    assert result["notes"] == []

    with pytest.raises(ValueError):
        find_similar_data("")


def test_find_similar_data_never_raises_on_source_failures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def broken(query: str, **kwargs: Any) -> list[DataCandidate]:
        raise RuntimeError("unexpected failure")

    def empty(query: str, **kwargs: Any) -> list[DataCandidate]:
        return []

    monkeypatch.setattr(sources, "search_local", empty)
    monkeypatch.setattr(sources, "search_zenodo", broken)
    monkeypatch.setattr(sources, "search_dissolve_examples", broken)
    result = find_similar_data("water")
    assert result["candidates"] == []
    assert result["notes"] == [
        "local: no matches",
        "dissolve: unexpected failure",
        "zenodo: unexpected failure",
    ]


# --------------------------------------------------------------------------- download
BODY = b"# sls18498.mdcs01\n#  H2O in 1mm can\n  0.175  -0.32  0.006\n"


def text_response(body: bytes = BODY, **headers: str) -> FakeResponse:
    return FakeResponse(
        body,
        {"Content-Type": "text/plain", "Content-Length": str(len(body)), **headers},
    )


def test_download_writes_the_file_and_creates_the_folder(
    network: FakeNetwork, tmp_path: Path
) -> None:
    url = "https://example.org/data/SLS18498-H2O.mint01"
    network.add(url, text_response)
    path = download(url, tmp_path / "new" / "sub")
    assert path == tmp_path / "new" / "sub" / "SLS18498-H2O.mint01"
    assert path.read_bytes() == BODY
    assert sorted(p.name for p in path.parent.iterdir()) == ["SLS18498-H2O.mint01"]
    # Downloading again replaces the file.
    network.add(url, lambda: text_response(b"# new\n"))
    assert download(url, path.parent).read_bytes() == b"# new\n"


def test_download_file_names(network: FakeNetwork, tmp_path: Path) -> None:
    zenodo = "https://zenodo.org/api/records/1/files/H2O%20data.mint01/content"
    network.add(zenodo, text_response)
    assert download(zenodo, tmp_path).name == "H2O data.mint01"

    url = "https://example.org/get?id=7"
    network.add(
        url,
        lambda: text_response(
            **{"Content-Disposition": 'attachment; filename="renamed.dat"'}
        ),
    )
    assert download(url, tmp_path).name == "renamed.dat"

    blob = f"https://github.com/{DISSOLVE_REPO}/blob/develop/examples/water/data/SLS18498-H2O.mint01"
    raw = f"{RAW_WATER}/SLS18498-H2O.mint01"
    network.add(raw, text_response)
    assert download(blob, tmp_path).name == "SLS18498-H2O.mint01"
    assert network.urls[-1] == raw


def test_download_refuses_bad_urls(tmp_path: Path) -> None:
    for url in ("ftp://example.org/a.mint01", "file:///etc/passwd", "", "not a url"):
        with pytest.raises(ValueError):
            download(url, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_download_size_limit(
    network: FakeNetwork, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    url = "https://example.org/big.mint01"
    network.add(
        url,
        lambda: text_response(
            **{"Content-Length": str(sources.MAX_DOWNLOAD_BYTES + 1)}
        ),
    )
    with pytest.raises(SourceError, match="over the limit"):
        download(url, tmp_path)
    assert list(tmp_path.iterdir()) == []

    # No Content-Length: the limit is enforced while streaming, nothing is left behind.
    monkeypatch.setattr(sources, "MAX_DOWNLOAD_BYTES", 10)
    monkeypatch.setattr(sources, "_DOWNLOAD_CHUNK", 4)
    network.add(url, lambda: FakeResponse(b"x" * 20, {"Content-Type": "text/plain"}))
    with pytest.raises(SourceError, match="over the limit"):
        download(url, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_download_refuses_web_pages_and_reports_http_errors(
    network: FakeNetwork, tmp_path: Path
) -> None:
    page = "https://zenodo.org/records/123"
    network.add(
        page,
        lambda: FakeResponse(
            b"<html></html>", {"Content-Type": "text/html; charset=utf-8"}
        ),
    )
    with pytest.raises(SourceError, match="web page"):
        download(page, tmp_path)

    missing = "https://example.org/missing.mint01"
    network.fail(missing, http_error(missing, 404, "Not Found"))
    with pytest.raises(SourceError, match="HTTP 404"):
        download(missing, tmp_path)

    down = "https://example.org/down.mint01"
    network.fail(down, urllib.error.URLError("Connection refused"))
    with pytest.raises(SourceError, match="Connection refused"):
        download(down, tmp_path)

    limited = "https://example.org/limited.mint01"
    network.fail(
        limited,
        http_error(limited, 429, "Too Many Requests", {"X-RateLimit-Remaining": "0"}),
    )
    with pytest.raises(RateLimitError):
        download(limited, tmp_path)
    assert list(tmp_path.iterdir()) == []


# --------------------------------------------------------------------------- network
_OFFLINE_MARKERS = (
    "unreachable",
    "timed out",
    "Temporary failure",
    "Name or service",
    "Network is unreachable",
    "Connection refused",
)


def _skip_if_offline(error: SourceError) -> None:
    if isinstance(error, RateLimitError) or any(
        m in str(error) for m in _OFFLINE_MARKERS
    ):
        pytest.skip(str(error))
    raise error


@pytest.mark.network
def test_dissolve_examples_online() -> None:
    sources.clear_caches()
    try:
        results = search_dissolve_examples("water", timeout=20)
    except SourceError as error:
        _skip_if_offline(error)
    cases = {candidate.extra["case"]: candidate for candidate in results}
    assert results[0].extra["case"] == "water"
    assert cases["water"].files == WATER_FILES
    assert cases["cu-btc"].files == ["Empty-CuBTC.mint01"]
    assert cases["silica_bulk"].files == ["SiO2_NIMROD.mint01"]
    assert set(cases["benzene"].files) == {"C6H6.mint01", "C6D6.mint01", "5050.mint01"}
    assert "argon" not in cases  # no data folder


@pytest.mark.network
def test_zenodo_online() -> None:
    try:
        results = search_zenodo("neutron diffraction liquid water", limit=3, timeout=20)
    except SourceError as error:
        _skip_if_offline(error)
    assert 0 < len(results) <= 3
    for candidate in results:
        assert candidate.source == "zenodo"
        assert candidate.title
        assert candidate.url.startswith("https://zenodo.org/")
    assert [c.score for c in results] == sorted(
        (c.score for c in results), reverse=True
    )


@pytest.mark.network
def test_download_online(tmp_path: Path) -> None:
    url = f"{RAW_WATER}/SLS18498-H2O.mint01"
    try:
        path = download(url, tmp_path, timeout=30)
    except SourceError as error:
        _skip_if_offline(error)
    assert path.name == "SLS18498-H2O.mint01"
    lines = path.read_text().splitlines()
    assert lines[0].startswith("#")
    assert lines[1].strip() == "#  H2O in 1mm can 5 in position 5 : 25C"
    assert path.stat().st_size > 10_000
