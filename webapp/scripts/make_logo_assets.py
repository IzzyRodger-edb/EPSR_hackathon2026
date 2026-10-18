"""Build the web app's logo images from docs/EPSR_Mellor_logo.png.

Run from webapp/:  uv run python scripts/make_logo_assets.py

The source is black and blue on white. The app needs it on light pages, dark
pages and a coloured header, so this writes into src/epsrweb/static/:

  logo.png        the full logo, background removed (light pages)
  logo-dark.png   the same with black turned near-white and the blue lightened
  logo-mark.png   the cube alone, in white with light-blue atoms (the header)
  favicon.png     the cube on a white tile, 64 x 64

The white background is removed exactly ("colour to alpha"), so logo.png over
white is the original image.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "docs" / "EPSR_Mellor_logo.png"
STATIC = HERE.parent / "src" / "epsrweb" / "static"


def unmix_white(rgb: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Colour and alpha such that colour over white, at that alpha, is `rgb`."""
    alpha = (1.0 - rgb).max(axis=2)
    with np.errstate(divide="ignore", invalid="ignore"):
        colour = 1.0 - (1.0 - rgb) / alpha[..., None]
    colour[alpha == 0] = 0.0
    return np.clip(colour, 0, 1), alpha


def is_blue(colour: np.ndarray) -> np.ndarray:
    return colour[..., 2] - colour[..., 0] > 0.25


def save(colour: np.ndarray, alpha: np.ndarray, path: Path, width: int | None = None) -> None:
    rgba = np.dstack([colour, alpha])
    image = Image.fromarray((rgba * 255 + 0.5).astype(np.uint8), "RGBA")
    if width:
        image = image.resize((width, round(image.height * width / image.width)), Image.LANCZOS)
    image.save(path, optimize=True)
    print(f"{path.name}: {image.width} x {image.height}, {path.stat().st_size // 1024} kB")


def main() -> None:
    rgb = np.asarray(Image.open(SOURCE).convert("RGB"), dtype=float) / 255.0
    colour, alpha = unmix_white(rgb)
    blue = is_blue(colour)

    # Trim to the drawing, with a small margin.
    rows, cols = np.nonzero(alpha > 0.02)
    top, bottom = max(rows.min() - 20, 0), min(rows.max() + 20, alpha.shape[0])
    left, right = max(cols.min() - 20, 0), min(cols.max() + 20, alpha.shape[1])
    crop = (slice(top, bottom), slice(left, right))
    save(colour[crop], alpha[crop], STATIC / "logo.png", width=1240)

    # Dark pages: neutrals flipped (black -> near-white), blue lightened.
    dark = colour.copy()
    grey = dark.mean(axis=2, keepdims=True)
    dark = np.where(blue[..., None], dark * 0.6 + 0.4, np.repeat(0.95 - 0.6 * grey, 3, axis=2))
    save(dark[crop], alpha[crop], STATIC / "logo-dark.png", width=1240)

    # The cube: everything left of the lettering.
    columns = np.nonzero((alpha > 0.02).any(axis=0))[0]
    gaps = np.nonzero(np.diff(columns) > 60)[0]
    cube_right = columns[gaps[0]] if len(gaps) else columns.max()
    cube_cols = slice(max(left, columns.min() - 12), cube_right + 12)
    cube_rows_any = np.nonzero((alpha[:, cube_cols] > 0.02).any(axis=1))[0]
    cube = (slice(cube_rows_any.min() - 12, cube_rows_any.max() + 12), cube_cols)

    mark = np.ones_like(colour)
    mark[blue] = (0.73, 0.83, 0.92)                     # light blue atoms
    save(mark[cube], alpha[cube], STATIC / "logo-mark.png", width=120)

    # Favicon: the original cube on a white tile.
    tile_size = 64
    cube_img = Image.fromarray((np.dstack([colour[cube], alpha[cube]]) * 255 + 0.5).astype(np.uint8), "RGBA")
    scale = (tile_size - 8) / max(cube_img.size)
    cube_img = cube_img.resize((round(cube_img.width * scale), round(cube_img.height * scale)), Image.LANCZOS)
    tile = Image.new("RGBA", (tile_size, tile_size), (255, 255, 255, 255))
    tile.alpha_composite(cube_img, ((tile_size - cube_img.width) // 2, (tile_size - cube_img.height) // 2))
    tile.save(STATIC / "favicon.png", optimize=True)
    print(f"favicon.png: {tile_size} x {tile_size}")


if __name__ == "__main__":
    main()
