"""The page font: always the bundled one, never the government's typeface.

The app looks like a GOV.UK service but is not one, so it must never render in
GDS Transport, the typeface reserved for government services. The stylesheet
bundles its own font under a private family name; these tests pin that down.
"""

from __future__ import annotations

import re
from pathlib import Path

import epsrweb

PACKAGE = Path(epsrweb.__file__).resolve().parent
CSS = (PACKAGE / "static" / "app.css").read_text(encoding="utf-8")
CSS_CODE = re.sub(r"/\*.*?\*/", "", CSS, flags=re.S)       # declarations only
FONT_FILE = "fonts/UbuntuSans-wdth-wght.ttf"
FORBIDDEN = ("transport", "gds")


def _font_values() -> list[str]:
    """Every value that can choose a font: font, font-family, the font tokens."""
    pattern = r"(?<![\w-])(font-family|font|--font|--mono)\s*:\s*([^;}]+)"
    return [value.strip() for _name, value in re.findall(pattern, CSS_CODE)]


def _font_faces() -> list[str]:
    return re.findall(r"@font-face\s*{([^}]*)}", CSS_CODE)


def test_no_font_declaration_names_the_government_typeface():
    values = _font_values()
    assert values, "found no font declarations at all; is the pattern stale?"
    for value in values + _font_faces():
        assert not any(word in value.lower() for word in FORBIDDEN), value


def test_the_page_font_is_the_bundled_one():
    token = re.search(r"--font\s*:\s*([^;]+);", CSS_CODE).group(1)
    assert token.split(",")[0].strip().strip('"') == "EPSR Sans"

    faces = [face for face in _font_faces() if '"EPSR Sans"' in face]
    assert len(faces) == 1
    # Only the bundled file: a local() source would let an installed font in.
    assert "local(" not in faces[0]
    assert FONT_FILE in faces[0]
    assert (PACKAGE / "static" / FONT_FILE).is_file()


def test_templates_leave_fonts_to_the_stylesheet():
    for template in (PACKAGE / "templates").glob("*.html"):
        text = template.read_text(encoding="utf-8").lower()
        assert "font-family" not in text, template.name
        assert "govuk-frontend" not in text, template.name
        assert not any(word in text for word in FORBIDDEN), template.name


def test_the_font_and_its_licence_are_served(client):
    font = client.get(f"/static/{FONT_FILE}")
    assert font.status_code == 200
    assert font.content[:4] == b"\x00\x01\x00\x00"           # a TrueType file
    assert client.get("/static/fonts/LICENCE-UbuntuFont.txt").status_code == 200
