"""End-to-end tests through the HTTP layer, including a complete job."""

from __future__ import annotations

import io
import time
import zipfile

import pytest

TERMINAL = {"done", "failed", "cancelled", "interrupted"}
PNG_MAGIC = bytes([0x89]) + b"PNG"


def _csrf(client, path="/new"):
    return client.get(path).text.split('name="csrf" value="')[1].split('"')[0]


def _upload(client, csrf, files, **extra):
    payload = {"csrf": csrf, "backend": "auto", "quick": "1", **extra}
    uploads = [("files", (path.name, path.read_bytes(), "text/plain")) for path in files]
    return client.post("/runs", data=payload, files=uploads, follow_redirects=False)


def _run_id(response):
    assert response.status_code == 303, response.text
    location = response.headers["location"]
    assert "error" not in location, location
    return location.split("/runs/")[1].split("?")[0]


def _await_run(client, run_id, timeout=300.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = client.get(f"/runs/{run_id}/status").json()
        if body["status"] in TERMINAL:
            return body
        time.sleep(0.5)
    pytest.fail(f"run {run_id} did not finish within {timeout} s")


# ------------------------------------------------------------------- pages


def test_public_pages_render(client):
    for path in ("/", "/about", "/instructions", "/login", "/register", "/healthz"):
        assert client.get(path).status_code == 200, path
    assert "Empirical Potential Structure Refinement" in client.get("/about").text
    assert "mint01" in client.get("/instructions").text


def test_pages_that_need_an_account_redirect(client):
    for path in ("/new", "/settings"):
        response = client.get(path, follow_redirects=False)
        assert response.status_code == 303
        assert "/login" in response.headers["location"]


def test_register_then_sign_in_and_out(client):
    created = client.post(
        "/register",
        data={"email": "a@b.org", "password": "a-long-password", "confirm": "a-long-password"},
        follow_redirects=False,
    )
    assert created.status_code == 303
    assert client.get("/new").status_code == 200

    client.post("/logout", data={"csrf": _csrf(client, "/")}, follow_redirects=False)
    assert client.get("/new", follow_redirects=False).status_code == 303

    denied = client.post("/login", data={"email": "a@b.org", "password": "wrong-password"},
                         follow_redirects=False)
    assert "error" in denied.headers["location"]
    allowed = client.post("/login", data={"email": "a@b.org", "password": "a-long-password"},
                          follow_redirects=False)
    assert allowed.status_code == 303
    assert client.get("/new").status_code == 200


def test_registration_rejects_short_passwords_and_duplicates(client):
    short = client.post("/register", data={"email": "c@d.org", "password": "short",
                                           "confirm": "short"}, follow_redirects=False)
    assert "at+least+10" in short.headers["location"].replace("%20", "+")
    client.post("/register", data={"email": "e@f.org", "password": "a-long-password",
                                   "confirm": "a-long-password"}, follow_redirects=False)
    again = client.post("/register", data={"email": "E@F.org", "password": "a-long-password",
                                           "confirm": "a-long-password"},
                        follow_redirects=False)
    assert "already+exists" in again.headers["location"].replace("%20", "+")


# ------------------------------------------------------------------ uploads


def test_upload_rejects_a_file_that_is_not_a_table(signed_in):
    client, csrf = signed_in
    response = client.post(
        "/runs",
        data={"csrf": csrf, "backend": "auto", "quick": "1"},
        files=[("files", ("notes.txt", b"this is not a diffraction pattern\n", "text/plain"))],
        follow_redirects=False,
    )
    assert response.status_code == 303
    assert "error" in response.headers["location"]


def test_upload_without_a_csrf_token_is_refused(signed_in, water_files):
    client, _ = signed_in
    response = client.post(
        "/runs",
        data={"csrf": "forged", "backend": "auto", "quick": "1"},
        files=[("files", (water_files[0].name, water_files[0].read_bytes(), "text/plain"))],
        follow_redirects=False,
    )
    assert response.status_code == 403


# ------------------------------------------------------------- a whole run


def test_a_run_goes_from_upload_to_figures_and_a_report(signed_in, water_files, settings):
    client, csrf = signed_in
    run_id = _run_id(_upload(client, csrf, water_files, title="Liquid water"))

    status = _await_run(client, run_id)
    assert status["status"] == "done", status
    assert status["figures"], "the run produced no figures"
    # The app's own seaborn figures, drawn from the tidy CSV tables.
    assert any(key.startswith("sb_") for key in status["figures"])

    page = client.get(f"/runs/{run_id}")
    assert page.status_code == 200
    assert "Liquid water" in page.text

    run_dir = settings.user_root(1) / run_id
    for name in ("report.pdf", "analysis.json", "tables/fits.csv", "tables/history.csv"):
        assert (run_dir / name).is_file(), f"{name} is missing"

    report = client.get(f"/runs/{run_id}/files/report.pdf")
    assert report.status_code == 200
    assert report.content.startswith(b"%PDF")

    figure = client.get(f"/runs/{run_id}/files/webfigures/sb_fit_overlay.png")
    assert figure.status_code == 200
    assert figure.content.startswith(PNG_MAGIC)
    assert client.get(f"/runs/{run_id}/files/webfigures/sb_fit_overlay.pdf").content.startswith(
        b"%PDF"
    )


def test_the_trajectory_and_other_paths_are_not_served(signed_in, water_files, settings):
    client, csrf = signed_in
    run_id = _run_id(_upload(client, csrf, water_files))
    _await_run(client, run_id)
    run_dir = settings.user_root(1) / run_id
    (run_dir / "trajectory").mkdir(exist_ok=True)
    (run_dir / "trajectory" / "history.tsv").write_text("secret", encoding="utf-8")
    (settings.data_root / "secret.key").write_text("key", encoding="utf-8")

    assert client.get(f"/runs/{run_id}/files/trajectory/history.tsv").status_code == 404
    assert client.get(f"/runs/{run_id}/files/events.jsonl").status_code == 404
    for attempt in ("../secret.key", "..%2Fsecret.key", "data/../../../secret.key"):
        assert client.get(f"/runs/{run_id}/files/{attempt}").status_code in (400, 404)

    # The user's own uploads are theirs to download again.
    page = client.get(f"/runs/{run_id}").text
    uploaded = page.split("/files/uploads/")[1].split('"')[0]
    assert client.get(f"/runs/{run_id}/files/uploads/{uploaded}").status_code == 200


def test_the_event_stream_replays_and_then_ends(signed_in, water_files):
    client, csrf = signed_in
    run_id = _run_id(_upload(client, csrf, water_files))
    _await_run(client, run_id)
    with client.stream("GET", f"/runs/{run_id}/events") as stream:
        body = "".join(chunk for chunk in stream.iter_text())
    assert "event: end" in body
    assert "read " in body  # the data files were announced in the log


def test_figures_zip_contains_latex_blocks(signed_in, water_files):
    client, csrf = signed_in
    run_id = _run_id(_upload(client, csrf, water_files))
    status = _await_run(client, run_id)
    assert status["status"] == "done"

    archive = client.post(
        f"/runs/{run_id}/figures.zip",
        data={"csrf": csrf, "figures": status["figures"][:2], "figure_format": "pdf"},
    )
    assert archive.status_code == 200
    with zipfile.ZipFile(io.BytesIO(archive.content)) as zipped:
        names = zipped.namelist()
        assert "figures.tex" in names
        assert any(name.startswith("figures/") for name in names)
        tex = zipped.read("figures.tex").decode("utf-8")
    assert chr(92) + "includegraphics" in tex
    # Names are prefixed with the run id so two runs cannot collide in Overleaf.
    assert run_id in tex


def test_zip_needs_a_selection(signed_in, water_files):
    client, csrf = signed_in
    run_id = _run_id(_upload(client, csrf, water_files))
    _await_run(client, run_id)
    assert client.post(f"/runs/{run_id}/figures.zip",
                       data={"csrf": csrf, "figure_format": "pdf"}).status_code == 400


def test_overleaf_push_without_a_token_is_refused(signed_in, water_files):
    client, csrf = signed_in
    run_id = _run_id(_upload(client, csrf, water_files))
    status = _await_run(client, run_id)
    refused = client.post(f"/runs/{run_id}/overleaf",
                          json={"figures": status["figures"][:1]},
                          headers={"X-CSRF-Token": csrf})
    assert refused.status_code == 400
    assert "token" in refused.json()["detail"].lower()


def test_chat_answers_after_the_run(signed_in, water_files):
    client, csrf = signed_in
    run_id = _run_id(_upload(client, csrf, water_files))
    _await_run(client, run_id)
    answered = client.post(f"/runs/{run_id}/chat", json={"message": "How good is the fit?"},
                           headers={"X-CSRF-Token": csrf})
    assert answered.status_code == 200
    assert answered.json()["answer"].strip()
    assert client.post(f"/runs/{run_id}/chat", json={"message": "hi"},
                       headers={"X-CSRF-Token": "forged"}).status_code == 403


def test_a_run_can_be_deleted_with_its_folder(signed_in, water_files, settings):
    client, csrf = signed_in
    run_id = _run_id(_upload(client, csrf, water_files))
    _await_run(client, run_id)
    run_dir = settings.user_root(1) / run_id
    assert run_dir.is_dir()
    deleted = client.post(f"/runs/{run_id}/delete", data={"csrf": csrf}, follow_redirects=False)
    assert deleted.status_code == 303
    assert not run_dir.exists()
    assert client.get(f"/runs/{run_id}", follow_redirects=False).status_code == 404


# ------------------------------------------------------------- separation


def test_one_user_cannot_reach_another_users_run(client, water_files):
    client.post("/register", data={"email": "one@x.org", "password": "a-long-password",
                                   "confirm": "a-long-password"}, follow_redirects=False)
    run_id = _run_id(_upload(client, _csrf(client), water_files))
    _await_run(client, run_id)

    client.post("/logout", data={"csrf": _csrf(client, "/")}, follow_redirects=False)
    client.post("/register", data={"email": "two@x.org", "password": "a-long-password",
                                   "confirm": "a-long-password"}, follow_redirects=False)

    assert client.get(f"/runs/{run_id}", follow_redirects=False).status_code == 404
    assert client.get(f"/runs/{run_id}/status").status_code == 404
    assert client.get(f"/runs/{run_id}/files/report.pdf").status_code == 404
