"""Test fixtures: a throwaway data root and a client with the job manager running."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
#: Real Gudrun files from the repository, so the upload path is exercised against
#: the format it will actually see.
WATER_DATA = sorted((REPO / "experimental_data" / "water" / "data").glob("*.mint01"))


@pytest.fixture()
def settings(tmp_path, monkeypatch):
    from epsrweb.config import load_settings

    monkeypatch.setenv("EPSRWEB_DATA_ROOT", str(tmp_path / "var"))
    monkeypatch.setenv("EPSRWEB_SECRET_KEY", "test-secret-key-not-for-production")
    monkeypatch.setenv("EPSRWEB_WORKERS", "1")
    monkeypatch.delenv("HF_TOKEN", raising=False)
    return load_settings()


@pytest.fixture()
def client(settings):
    from fastapi.testclient import TestClient

    from epsrweb.app import create_app

    with TestClient(create_app(settings)) as test_client:
        yield test_client


@pytest.fixture()
def signed_in(client):
    """A registered, signed-in client. Returns (client, csrf token)."""
    response = client.post(
        "/register",
        data={"email": "chemist@example.org", "password": "a-long-password",
              "confirm": "a-long-password"},
        follow_redirects=False,
    )
    assert response.status_code == 303, response.text
    page = client.get("/new")
    assert page.status_code == 200
    return client, _csrf_from(page.text)


def _csrf_from(html: str) -> str:
    marker = 'name="csrf" value="'
    start = html.index(marker) + len(marker)
    return html[start:html.index('"', start)]


@pytest.fixture()
def csrf_from():
    return _csrf_from


@pytest.fixture()
def water_files():
    if not WATER_DATA:
        import pytest as _pytest

        _pytest.skip("the repository's water .mint01 files are not present")
    return WATER_DATA[:2]
