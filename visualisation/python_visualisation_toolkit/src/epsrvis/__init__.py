"""epsrvis — analysis and visualisation of EPSR refinements.

Written to be embedded in the Python wrapper of the Rust EPSR port (e.g. as
``epsr.vis``); it only uses relative imports and numpy/scipy/matplotlib.
The 3D viewer (:mod:`epsrvis.viewer`) additionally needs pyvista and is not
imported here.

Main entry points
-----------------
:class:`Trajectory`          per-iteration record (load from disk / build in memory)
:class:`TrajectoryRecorder`  record a running wrapper ``Simulation``
:func:`record_run`           one recorded ``epsr`` invocation (== ``sim.run()``)
:func:`make_report`          all figures + summary tables
:mod:`epsrvis.metrics`       R-factor, R², χ², acceptance, coordination, angles
:mod:`epsrvis.plots`         individual matplotlib figures
:mod:`epsrvis.outputs`       readers for the EPSR output files (.?01, .erg)
:mod:`epsrvis.viewer`        pyvista 3D configuration + partial g(r)
"""

from .outputs import OutputFolder, read_epsr_table, read_erg
from .recorder import TrajectoryRecorder, record_run
from .report import make_report
from .trajectory import Trajectory

__version__ = "0.1.0"

__all__ = ["Trajectory", "TrajectoryRecorder", "record_run", "make_report", "OutputFolder",
           "read_epsr_table", "read_erg"]
