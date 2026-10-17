"""A simulation folder and its EPSR naming conventions."""

from __future__ import annotations

import shutil
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING

from .io import KeywordFile, Table, read_keyword_file, read_table_series
from .simulation import Overrides, Simulation, find_inp

if TYPE_CHECKING:
    from .refine import History

__all__ = ["Project", "copy_folder"]


def copy_folder(
    source: str | PathLike[str],
    destination: str | PathLike[str],
    *,
    overwrite: bool = False,
) -> Path:
    """Copy the regular files of ``source`` (not subfolders) into a new folder ``destination``.

    Refuses to touch an existing non-empty ``destination`` unless ``overwrite``
    is true, in which case ``destination`` is deleted first.
    """
    source, destination = Path(source), Path(destination)
    if not source.is_dir():
        raise FileNotFoundError(f"template folder {source} does not exist")
    if destination.exists() and any(destination.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"{destination} exists and is not empty; pass overwrite=True to replace it"
            )
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    for source_file in source.iterdir():
        if source_file.is_file():
            shutil.copy2(source_file, destination / source_file.name)
    return destination


class Project:
    """An EPSR simulation folder.

    >>> project = Project.from_template("EPSR-test", "work/dmso")  # copy, then use the copy
    >>> simulation = project.simulation(overrides={"ntimes": 5})
    >>> simulation.run()

    Every EPSR run rewrites ``.ato``, ``.pcof``, ``.EPSR.inp`` and ``.EPSR.inpa``
    in the folder and appends to ``.erg``. Always work on a copy of your
    original folder; :meth:`snapshot` keeps intermediate states.
    """

    def __init__(
        self, folder: str | PathLike[str], inp_name: str | None = None
    ) -> None:
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
        return self.inp_name.removesuffix(".inp")

    def simulation(
        self, overrides: Overrides = None, *, verbose: bool = False
    ) -> Simulation:
        """Read the folder into a new :class:`Simulation` (one EPSR run)."""
        return Simulation(
            self.folder, self.inp_name, overrides=overrides, verbose=verbose
        )

    def snapshot(self, label: str, *, overwrite: bool = False) -> Path:
        """Copy the folder's files to ``<folder>.snapshots/<label>/`` and return that path."""
        snapshot_folder = self.folder.parent / f"{self.folder.name}.snapshots" / label
        return copy_folder(self.folder, snapshot_folder, overwrite=overwrite)

    def output(self, extension: str) -> Path:
        """Path of an output file: ``output("g01")`` is ``<folder>/<run_name>.g01``."""
        return self.folder / f"{self.run_name}.{extension.lstrip('.')}"

    def read_table(self, extension: str) -> Table:
        """Read a result table (``"g01"``, ``"u01"``...).

        Continuation files (``02``, ``03``...) are merged in.
        """
        return read_table_series(self.output(extension))

    def history(self) -> History:
        """The ``.erg`` energy history as a :class:`epsr.History`."""
        from .refine import History

        return History.from_erg(self.output("erg"))

    def input_file(self) -> KeywordFile:
        """The parsed ``.EPSR.inp`` file."""
        return read_keyword_file(self.folder / self.inp_name)

    def files(self) -> dict[str, Path | list[Path]]:
        """The main input files named by the ``.EPSR.inp`` file."""
        input_file = self.input_file()
        paths: dict[str, Path | list[Path]] = {"inp": self.folder / self.inp_name}
        for key, keyword in (("ato", "fnameato"), ("pcof", "fnamepcof")):
            file_name = input_file.words(keyword)
            if file_name:
                paths[key] = self.folder / file_name[0]
        paths["data"] = [
            self.folder / block["datafile"].split()[0]
            for block in input_file.blocks
            if block.get("datafile")
        ]
        paths["weights"] = [
            self.folder / block["wtsfile"].split()[0]
            for block in input_file.blocks
            if block.get("wtsfile")
        ]
        paths["erg"] = self.output("erg")
        paths["out"] = self.output("out")
        return paths
