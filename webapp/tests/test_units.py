"""Unit tests for the parts that do not need a server."""

from __future__ import annotations

import zipfile
from io import BytesIO

import pytest

from epsrweb import overleaf, security, storage

GOOD = """# NIMROD00024401.mint01
# Water H2O on SANDALS
#     1  2172
   0.1000000E+00  0.1000000E+00  0.1000000E-01
   0.2000000E+00  0.1500000E+00  0.1000000E-01
   0.3000000E+00  0.2000000E+00  0.1000000E-01
   0.4000000E+00  0.2500000E+00  0.1000000E-01
   0.5000000E+00  0.3000000E+00  0.1000000E-01
   0.6000000E+00  0.3500000E+00  0.1000000E-01
   0.7000000E+00  0.4000000E+00  0.1000000E-01
   0.8000000E+00  0.4500000E+00  0.1000000E-01
   0.9000000E+00  0.5000000E+00  0.1000000E-01
   0.1000000E+01  0.5500000E+00  0.1000000E-01
   0.1100000E+01  0.0000000E+00  0.0000000E+00
"""


# ----------------------------------------------------------------- passwords


def test_password_round_trip():
    stored = security.hash_password("a-long-password")
    assert stored.startswith("scrypt$")
    assert security.verify_password("a-long-password", stored)
    assert not security.verify_password("a-long-passwore", stored)
    assert not security.verify_password("", stored)


def test_verify_password_rejects_rubbish():
    assert not security.verify_password("x", "not-a-hash")
    assert not security.verify_password("x", "scrypt$bad")


def test_session_token_round_trip():
    secret = b"secret"
    token = security.make_session_token(secret, 42, 3600)
    assert security.read_session_token(secret, token) == 42
    assert security.read_session_token(b"other", token) is None
    assert security.read_session_token(secret, token + "x") is None
    assert security.read_session_token(secret, None) is None


def test_expired_session_is_rejected():
    secret = b"secret"
    assert security.read_session_token(secret, security.make_session_token(secret, 1, -10)) is None


def test_csrf_token_is_bound_to_the_session():
    secret = b"secret"
    session = security.make_session_token(secret, 7, 3600)
    token = security.make_csrf_token(secret, session)
    assert security.check_csrf_token(secret, session, token)
    assert not security.check_csrf_token(secret, "another-session", token)
    assert not security.check_csrf_token(secret, session, "forged")
    assert not security.check_csrf_token(secret, None, token)


def test_stored_secrets_are_encrypted_and_recoverable():
    secret = b"secret-key-material"
    box = security.encrypt_secret(secret, "olp_sensitive_token")
    assert "olp_sensitive_token" not in box
    assert security.decrypt_secret(secret, box) == "olp_sensitive_token"
    assert security.decrypt_secret(b"wrong-key", box) is None


# ------------------------------------------------------------ upload checks


def test_validate_mint01_reads_the_title_from_the_second_header():
    info = storage.validate_mint01(GOOD, name="water.mint01")
    assert info.title == "Water H2O on SANDALS"
    assert info.points == 11
    assert info.q_min == pytest.approx(0.1)
    assert info.warnings == ()


def test_validate_mint01_rejects_non_tables():
    with pytest.raises(ValueError, match="three numeric columns"):
        storage.validate_mint01("# header\n1.0 2.0\nnonsense here now\n", name="bad")
    with pytest.raises(ValueError, match="not a .mint01"):
        storage.validate_mint01("# header\n1.0 2.0 3.0\n", name="short")


def test_validate_mint01_warns_about_silently_discarded_q():
    rows = "".join(f"{0.05 * i:.4f} 1.0 0.1\n" for i in range(1, 800))
    info = storage.validate_mint01("# f\n# title\n" + rows, name="wide")
    assert any("discarded" in warning for warning in info.warnings)
    assert any("sentinel" in warning for warning in info.warnings)


def test_data_file_names_are_blank_free():
    # EPSR truncates a data path at the first blank, so this matters.
    assert " " not in storage.safe_data_filename("my data file.mint01")
    assert storage.safe_data_filename("../../etc/passwd") == "passwd.mint01"
    assert storage.safe_data_filename("x/y/SLS18498-H2O.mint01") == "SLS18498-H2O.mint01"


def test_resolve_within_refuses_escapes(tmp_path):
    (tmp_path / "figures").mkdir()
    assert storage.resolve_within(tmp_path, "figures") == (tmp_path / "figures").resolve()
    for attempt in ("../secrets", "figures/../../secrets", "/etc/passwd"):
        with pytest.raises(ValueError):
            storage.resolve_within(tmp_path, attempt)


# ----------------------------------------------------------------- Overleaf


@pytest.mark.parametrize("given, expected", [
    ("https://www.overleaf.com/project/65f1c0ffee1234567890abcd", "65f1c0ffee1234567890abcd"),
    ("https://git.overleaf.com/65f1c0ffee1234567890abcd", "65f1c0ffee1234567890abcd"),
    ("65f1c0ffee1234567890abcd", "65f1c0ffee1234567890abcd"),
    ("https://www.overleaf.com/project/65f1c0ffee1234567890abcd/", "65f1c0ffee1234567890abcd"),
])
def test_parse_project_id(given, expected):
    assert overleaf.parse_project_id(given) == expected


@pytest.mark.parametrize("given", ["", "   ", "https://example.com/", "a b c", "--upload-pack=x"])
def test_parse_project_id_rejects_junk(given):
    with pytest.raises(overleaf.OverleafError):
        overleaf.parse_project_id(given)


def test_redact_removes_tokens_and_credentials():
    assert "olp_s3cret" not in overleaf.redact("fatal: olp_s3cret rejected", "olp_s3cret")
    assert overleaf.redact("https://git:olp_x@git.overleaf.com/abc", None).startswith(
        "https://***@"
    )


def test_latex_snippet_escapes_the_caption():
    snippet = overleaf.latex_snippet("run_01.pdf", "R&D at 50% with an_underscore", "run_01")
    backslash = chr(92)
    assert snippet.startswith(backslash + "begin{figure}")
    assert "figures/run_01.pdf" in snippet
    assert backslash + "&" in snippet and backslash + "%" in snippet
    assert backslash + "label{fig:run_01}" in snippet


def test_build_zip_carries_the_figures_and_the_snippets(tmp_path):
    figure = tmp_path / "sb_fit_overlay.pdf"
    figure.write_bytes(b"%PDF-1.4 fake")
    payload = overleaf.build_zip([(figure, "run42_sb_fit_overlay.pdf", "The fit")])
    with zipfile.ZipFile(BytesIO(payload)) as archive:
        names = archive.namelist()
        assert "figures/run42_sb_fit_overlay.pdf" in names
        assert "figures.tex" in names
        tex = archive.read("figures.tex").decode("utf-8")
    assert "figures/run42_sb_fit_overlay.pdf" in tex
    assert "The fit" in tex


def test_push_requires_a_token(tmp_path):
    with pytest.raises(overleaf.OverleafError, match="No Overleaf token"):
        overleaf.push_figures("65f1c0ffee1234567890abcd", "", [(tmp_path, "a.pdf", "c")],
                              commit_message="x")
