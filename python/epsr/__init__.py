"""Python interface to the Rust port of EPSR (Empirical Potential Structure Refinement).

EPSR (A. K. Soper) builds an atomistic model of a liquid or amorphous material that
is consistent with measured diffraction data: a Monte Carlo simulation with a
reference potential is driven towards the data by an empirical potential refined
from the difference between measured and simulated structure factors.

Quick start::

    import epsr
    proj = epsr.Project.from_template("EPSR-test", "work/dmso")   # copy the folder
    sim = proj.simulation(overrides={"ntimes": 5})
    sim.run()                                  # one EPSR run (files updated in place)
    sim.summary.r_factor, sim.gr_inter(), sim.potential.empirical()
    hist = epsr.refine(proj, n_runs=20)        # repeated runs

See ``python/docs/`` for the user guide and API reference.
"""

from __future__ import annotations

from . import io, plot, tools
from ._core import VERSION_DATE
from ._core import __version__ as __version__
from .errors import EpsrError, EpsrInputError, EpsrIOError, EpsrUnsupportedError
from .io import read_erg, read_keyword_file, read_table
from .potential import Potential, PotParams, TypeParams
from .project import Project, copy_folder
from .refine import History, refine, refine_many
from .simulation import Control, Dataset, Fit, Simulation, Summary, render_overrides

__all__ = [
    "VERSION_DATE",
    "Control",
    "Dataset",
    "EpsrError",
    "EpsrIOError",
    "EpsrInputError",
    "EpsrUnsupportedError",
    "Fit",
    "History",
    "PotParams",
    "Potential",
    "Project",
    "Simulation",
    "Summary",
    "TypeParams",
    "copy_folder",
    "io",
    "plot",
    "read_erg",
    "read_keyword_file",
    "read_table",
    "refine",
    "refine_many",
    "render_overrides",
    "tools",
]
