"""A simulation folder and its EPSR naming conventions."""

from __future__ import annotations

import shutil
from os import PathLike
from pathlib import Path

from .io import KeywordFile, Table, read_keyword_file, read_table_series
from .simulation import Overrides, Simulation, find_inp

__all__ = ["Project", "copy_folder"]


def copy_folder(src: str | PathLike[str], dst: str | PathLike[str], *, overwrite: bool = False) -> Path:
    """Copy the regular files of ``src`` (not subfolders) into a new folder ``dst``.

    Refuses to touch an existing non-empty ``dst`` unless ``overwrite`` is true,
    in which case ``dst`` is deleted first.
    """
    src, dst = Path(src), Path(dst)
    if not src.is_dir():
        raise FileNotFoundError(f"template folder {src} does not exist")
    if dst.exists() and any(dst.iterdir()):
        if not overwrite:
            raise FileExistsError(f"{dst} exists and is not empty; pass overwrite=True to replace it")
        shutil.rmtree(dst)
    dst.mkdir(parents=True, exist_ok=True)
    for p in src.iterdir():
        if p.is_file():
            shutil.copy2(p, dst / p.name)
    return dst


class Project:
    """An EPSR simulation folder.

    >>> proj = Project.from_template("EPSR-test", "work/dmso")   # copy, then work on the copy
    >>> sim = proj.simulation(overrides={"ntimes": 5})
    >>> sim.run()

    Every EPSR run rewrites ``.ato``, ``.pcof``, ``.EPSR.inp`` and ``.EPSR.inpa``
    in the folder and appends to ``.erg``. Always work on a copy of your
    original folder; :meth:`snapshot` keeps intermediate states.
    """

    def __init__(self, folder: str | PathLike[str], inp_name: str | None = None) -> None:
        self.folder = Path(folder)
        if not self.folder.is_dir():
            raise FileNotFoundError(f"simulation folder {self.folder} does not exist")
        self.inp_name = inp_name if inp_name is not None else find_inp(self.folder)

    def __repr__(self) -> str:
        return f"<epsr.Project {self.folder / self.inp_name}>"

    @classmethod
    def from_template(
        cls,
        template: str | PathLike[str],
        work: str | PathLike[str],
        *,
        inp_name: str | None = None,
        overwrite: bool = False,
    ) -> Project:
        """Copy the files of ``template`` into ``work`` and open the copy."""
        return cls(copy_folder(template, work, overwrite=overwrite), inp_name)

    @property
    def run_name(self) -> str:
        """Prefix of the output files, e.g. ``200DMSO400Water.EPSR``."""
        return self.inp_name[: -len(".inp")] if self.inp_name.endswith(".inp") else self.inp_name

    def simulation(self, overrides: Overrides = None, *, verbose: bool = False) -> Simulation:
        """Read the folder into a new :class:`Simulation` (one EPSR run)."""
        return Simulation(self.folder, self.inp_name, overrides=overrides, verbose=verbose)

    def snapshot(self, label: str, *, overwrite: bool = False) -> Path:
        """Copy the folder's files to ``<folder>.snapshots/<label>/`` and return that path."""
        dst = self.folder.parent / f"{self.folder.name}.snapshots" / label
        return copy_folder(self.folder, dst, overwrite=overwrite)

    def output(self, ext: str) -> Path:
        """Path of an output file: ``output("g01")`` → ``<folder>/<run_name>.g01``."""
        return self.folder / f"{self.run_name}.{ext.lstrip('.')}"

    def read_table(self, ext: str) -> Table:
        """Read a result table (``"g01"``, ``"f01"``, ``"u01"``...), merging ``02``, ``03``... files."""
        return read_table_series(self.output(ext))

    def history(self):  # type: ignore[no-untyped-def]  # -> History (import cycle)
        """The ``.erg`` energy history as a :class:`epsr.History`."""
        from .refine import History

        return History.from_erg(self.output("erg"))

    def input_file(self) -> KeywordFile:
        """The parsed ``.EPSR.inp`` file."""
        return read_keyword_file(self.folder / self.inp_name)

    def files(self) -> dict[str, Path | list[Path]]:
        """The main input files named by the ``.EPSR.inp`` file."""
        inp = self.input_file()
        out: dict[str, Path | list[Path]] = {"inp": self.folder / self.inp_name}
        for key, name in (("ato", "fnameato"), ("pcof", "fnamepcof")):
            words = inp.words(name)
            if words:
                out[key] = self.folder / words[0]
        out["data"] = [self.folder / b["datafile"].split()[0] for b in inp.blocks if b.get("datafile")]
        out["weights"] = [self.folder / b["wtsfile"].split()[0] for b in inp.blocks if b.get("wtsfile")]
        out["erg"] = self.output("erg")
        out["out"] = self.output("out")
        return out
