"""EPSR input files generated from a :class:`SystemSpec`.

Everything EPSR reads is produced here, as text, so it can be inspected, embedded in
the standalone script and written to the simulation folder:

* ``<molecule>.mol`` for ``makemole`` (gmol format, see ``validation/water/Water.mol``);
* the answers to ``mixato`` (number of types; per type the template name and count;
  density; box name);
* ``<label>.NWTS.dat`` for ``nwts`` (keyword file; use ``epsr.tools.keywords("nwts")``
  for the keyword table and descriptions, and ``write_keyword_file``'s layout);
* ``<name>.EPSR.inp`` from the template ``validation/water/waterbox.EPSR.inp`` with the
  run name, ``fnameato``/``fnamepcof``, ``rho``, ``ndata`` and one ``data n`` block
  per data set, ``qmin`` and the other ``epsr_keywords`` substituted;
* for ``pdb`` systems the answers to ``readpdb`` and the ``changeato`` commands
  (see ``validation/cu-btc/readpdb.stdin`` and ``changeato.stdin``), and the
  ``vanderWaalsRadii.txt`` table ``readpdb`` reads from the folder to find the bonds.

The ``.mol`` text uses the layout ``makemole`` itself writes when it rewrites the
file (Fortran formats, reals rounded to single precision), so the generated file is
the one EPSR would produce from the same numbers; ``validation/water/Water.mol`` is
reproduced byte for byte.
"""

from __future__ import annotations

import math
import shutil
import struct
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from epsr import tools as epsr_tools

from .spec import DATASET_KEYWORDS, Dataset, Molecule, SpecError, SystemSpec

__all__ = [
    "EpsrInputs",
    "generate",
    "write",
    "render_mol",
    "render_nwts",
    "render_inp",
    "mixato_answers",
    "readpdb_answers",
    "changeato_commands",
    "auto_rshmin",
    "EPSR_INP_TEMPLATE",
    "EPSR_ELEMENT_MASSES",
    "VAN_DER_WAALS_RADII",
    "VDW_RADII_FILE",
    "render_vdw_radii",
    "pdb_atom_labels",
    "format_value",
]

#: The ``.EPSR.inp`` template: every global keyword with the test-case default and its
#: description, in file order, plus the per-dataset block layout. Copied from
#: ``validation/water/waterbox.EPSR.inp`` (one ``data`` block kept as the layout of
#: every data set). ``render_inp`` substitutes the title, ``fnameato``,
#: ``fnamepcof``, ``rho``, ``ntimes``, ``ndata``, the data blocks and the
#: ``epsr_keywords`` of the specification; ``ireset 1``, ``iinit 1``, ``nsumt -1``
#: and ``niter 1`` are the scratch start.
EPSR_INP_TEMPLATE: str = """\
waterbox.EPSR               Title of this file
feedback    0.9               Confidence factor - should be < 1. [0.8]
potfac      1.000000   1.000000               >0.0 to enable potential refinement, 0.0 to inhibit
rspcrmin    1.00               Minimum distance for calculating the R-space coefficient [1.00]
rspcfrac    0.5               Fraction of R-space coefficient in control level [0.2]
num_threds  0               No. of parallel threads to be used, (0 to let program choose) [0]
porodpwr    0 0               Q-value and power for Porod scattering [0 0]
nmolcell    10               Average number of molecules in a cell [10]
control-p   none 0 0               Control pressure [a, b and/or c axes control plus value plus step, none 0 0]
ref_intra   0.000    0.000               Weighting on EP and Coulomb terms for intra-molecular structure [0 0]
sizefactor  1.00000    0.90000  0.00000E+00               Multiplying factor for box dimension, decline rate, and threshold. [1.0 0.9 0.0]
nq          600               Number of Q values. [600]
qstep       0.05               Size of Q step [1/A]. [0.05]
qinscalefc  1.0 0.0               Scale factor on input Q values and broadening power for self terms. [1.0 0.0]
ireset      1               1: complete reset; 2: sets the Empirical Potential to zero
iinit       1               Sets accumulators to zero. Recalculates r and Q. [1]
ntimes      5               Number of MC cycles between potential refinements. [5]
niter       1               Number of potential refinements before exitting. [1]
nsumt       -1               Number of iterations already accumulated. [-1 with reset]
intra       100               Number of iterations between molecule shakes. [100]
rotfreq     5               Number of iterations between internal rotation moves. [5]
inter       5               Number of iterations in running averages. [5]
rho         8.80699679E-02               Atomic number density - will be derived from .ato file
cellst      0.03               Size of r step [A]. [0.03]
rmaxgr      0.000000E+00               Range of g(r) and F.T. (0.0 will use half the cell box) [0.0]
ngrsamples  0               Requested no. of origin molecules to sample g(r). (0 will use 1000 molecules) [0]
fwhm        0.0               Resolution width - Q independent term. [0.0]
fwhmq       0.02               Resolution width - Q dependent term. [0.02 for SLS]
nsmoop      1               1 means background subtraction is ON, 0 means OFF
fnameato    waterbox.ato               Name of .ato file
fnamepcof   waterbox.pcof               Name of potential coefficients file.
revlorch    0.0 0.0               Broadening factor in Q space. [0.0 0.0]
qwidthqmax  0.01 0               Broadening and maximum Q for Bragg peak calculation
mplicities  1 1 1               No. of unit cells along a, b and c for Bragg peak calculation
hklqmin     0.0 0.0 0.0               Minimum value of qhkl to be used, minimum radius for Bragg g(r), and Debye-Waller factor
diffuse     0               No. of unit cells along a, b and c for diffuse scattering calculation, maximum l, steps in l and m [0]
nlatrefine  0               No. of lattice refinement steps [0]
latrefine   0 0 0 0 0 0               Refinement steps for each of a, b, c, alpha, beta, gamma [0 0 0 0 0 0]
rejrate     0.75               Rejection rate [0.75]
qmin        0.05000    0.00000               Qmin for Fourier transforms and for potential fits. [0.05 0.0]
ndata       3               Number of data files to be fit by EPSR

data   1

datafile    SLS18502-D2O.mint01               Name of data file to be fit
wtsfile     D2O.NWTStot.wts               Name of weights file for this data set
nrtype      5               Data type - see User Manual for more details
rshmin      0.7               Minimum radius [A] - used for background subtraction
szeros      0.0               Zero limit - 0 means use first data point for Q=0
tweak       1.0               Scaling factor for this data set. [1.0]
efilereq    1.0               Requested energy amplitude for this data set [1.0]
porodadd    0.0               Additional amplitude of Porod scattering for this dataset [0.0]

q
"""

#: EPSR's own atomic weights (``chemical_symbols.inc``), which ``readpdb`` assigns
#: from the element of every atom. ``changeato`` only needs ``atmass`` when the
#: specification's mass differs from this value.
EPSR_ELEMENT_MASSES: dict[str, float] = {
    "H": 1.00797, "He": 4.00260, "Li": 6.94090, "Be": 9.01218, "B": 10.81003,
    "C": 12.01110, "N": 14.00670, "O": 15.99940, "F": 18.99840, "Ne": 20.17900,
    "Na": 22.98977, "Mg": 24.30500, "Al": 26.98154, "Si": 28.08550, "P": 30.97376,
    "S": 32.06600, "Cl": 35.45300, "Ar": 39.94800, "K": 39.09830, "Ca": 40.07800,
    "Sc": 44.95591, "Ti": 47.88000, "V": 50.94150, "Cr": 51.99600, "Mn": 54.93800,
    "Fe": 55.84700, "Co": 58.93320, "Ni": 58.67000, "Cu": 63.54600, "Zn": 65.39000,
    "Ga": 69.72300, "Ge": 72.59000, "As": 74.92159, "Se": 78.96000, "Br": 79.90400,
    "Kr": 83.80000, "Rb": 85.46780, "Sr": 87.62000, "Y": 88.90585, "Zr": 91.22400,
    "Nb": 92.90638, "Mo": 95.94000, "Tc": 98.90625, "Ru": 101.07000, "Rh": 102.90550,
    "Pd": 106.42000, "Ag": 107.86820, "Cd": 112.41000, "In": 114.82000,
    "Sn": 118.71000, "Sb": 121.75000, "Te": 127.60000, "I": 126.90447,
    "Xe": 131.28999, "Cs": 132.90543, "Ba": 137.33000, "La": 138.90550,
    "Ce": 140.12000, "Pr": 140.90765, "Nd": 144.24001, "Pm": 145.00000,
    "Sm": 150.36000, "Eu": 151.96001, "Gd": 157.25000, "Tb": 158.92534,
    "Dy": 162.50000, "Ho": 164.93031, "Er": 167.25999, "Tm": 168.93420,
    "Yb": 173.03999, "Lu": 174.96700, "Hf": 178.49001, "Ta": 180.94791,
    "W": 183.95093, "Re": 186.20700, "Os": 190.20000, "Ir": 192.22000,
    "Pt": 195.08000, "Au": 196.96654, "Hg": 200.59000, "Tl": 204.38300,
    "Pb": 207.20000, "Bi": 208.98038, "Po": 209.00000, "At": 210.00000,
    "Rn": 222.00000, "Fr": 223.00000, "Ra": 226.02541, "Ac": 227.00000,
    "Th": 232.03806, "Pa": 231.03587, "U": 238.02890, "Np": 237.04817,
    "Pu": 244.06274, "Am": 243.00000, "Cm": 247.00000, "Bk": 247.00000,
    "Cf": 251.00000, "Es": 254.00000, "Fm": 257.00000, "Md": 256.00000,
    "No": 254.00000, "Lr": 258.00000,
}  # fmt: skip

#: Van der Waals radii (Å) by element, the ``vanderWaalsRadii.txt`` table ``readpdb``
#: reads from the simulation folder: two atoms closer than the van der Waals factor
#: times the sum of their radii are bonded (``validation/cu-btc/vanderWaalsRadii.txt``).
VAN_DER_WAALS_RADII: dict[str, float] = {
    "H": 0.0, "He": 2, "Li": 2.2, "Be": 1.9, "B": 1.8, "C": 1.7, "N": 1.6, "O": 1.55,
    "F": 1.5, "Ne": 1.5, "Na": 2.4, "Mg": 2.2, "Al": 2.1, "Si": 2.1, "P": 1.95,
    "S": 1.8, "Cl": 1.8, "Ar": 1.8, "K": 2.8, "Ca": 2.4, "Sc": 2.3, "Ti": 2.15,
    "V": 2.05, "Cr": 2.05, "Mn": 2.05, "Fe": 2.05, "Co": 2, "Ni": 2, "Cu": 2,
    "Zn": 2.1, "Ga": 2.1, "Ge": 2.1, "As": 2.05, "Se": 1.9, "Br": 1.9, "Kr": 1.9,
    "Rb": 2.9, "Sr": 2.55, "Y": 2.4, "Zr": 2.3, "Nb": 2.15, "Mo": 2.1, "Tc": 2.05,
    "Ru": 2.05, "Rh": 2, "Pd": 2.05, "Ag": 2.1, "Cd": 2.2, "In": 2.2, "Sn": 2.25,
    "Sb": 2.2, "Te": 2.1, "I": 2.1, "Xe": 2.1, "Cs": 3, "Ba": 2.7, "La": 2.5,
    "Ce": 2.5, "Pr": 2.5, "Nd": 2.5, "Pm": 2.5, "Sm": 2.5, "Eu": 2.5, "Gd": 2.5,
    "Tb": 2.5, "Dy": 2.5, "Ho": 2.5, "Er": 2.5, "Tm": 2.5, "Yb": 2.5, "Lu": 2.5,
    "Hf": 2.25, "Ta": 2.2, "W": 2.1, "Re": 2.05, "Os": 2, "Ir": 2, "Pt": 2.05,
    "Au": 2.1, "Hg": 2.05, "Tl": 2.2, "Pb": 2.3, "Bi": 2.3, "Po": 2.3, "At": 2.3,
    "Rn": 2.3, "Fr": 2.3, "Ra": 2.3, "Ac": 2.3, "Th": 2.4, "Pa": 2.4, "U": 2.3,
    "Np": 2.3, "Pu": 2.3, "Am": 2.3, "Cm": 2.3, "Bk": 2.3, "Cf": 2.3, "Es": 2.3,
    "Fm": 2.3, "Md": 2.3, "No": 2.3, "Lr": 2.3,
}  # fmt: skip
#: The file name ``readpdb`` looks for (in the current directory).
VDW_RADII_FILE = "vanderWaalsRadii.txt"

#: ``rshmin`` defaults (Å): the EPSR test-case value for molecular systems, and a
#: value below the shortest expected distance for atomic systems (silica Si-O 1.6 Å).
MOLECULAR_RSHMIN = 0.7
ATOMIC_RSHMIN = 1.2
#: Fixed ``.mol`` keywords (EPSR test-case values).
MOL_VIBTEMP = 65.0
MOL_ANGTEMP = 1.0
MOL_DIHTEMP = 10.0
MOL_DENSITY = 0.001
#: ``makemole`` labels the atoms of a template ``001``..``999``.
MAX_ATOMS_PER_MOLECULE = 999
#: The ``changeato`` step-size commands, in the order of ``validation/cu-btc``.
CHANGEATO_STEP_KEYS = ("stepmi", "stepri", "stepr", "stepm")
#: Keyword files: keyword padded to 11, a blank, the value, 15 blanks, description.
KEYWORD_WIDTH = epsr_tools.KEYWORD_COLUMN_WIDTH
DESCRIPTION_GAP = epsr_tools.DESCRIPTION_GAP
TITLE_MARKER = "Title of this file"


# --------------------------------------------------------------------------- numbers
def _f32(value: float) -> float:
    """``value`` rounded to single precision, as EPSR stores every real."""
    return struct.unpack("f", struct.pack("f", float(value)))[0]


def format_value(value: Any) -> str:
    """A keyword or answer value as EPSR reads it: lists joined with blanks,
    booleans ``T``/``F``, floats in their shortest exact form (``0.1``, not
    ``0.10000``; exponents written as ``E``)."""
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, (list, tuple)):
        return " ".join(format_value(item) for item in value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"cannot write {value!r} to an EPSR file")
        return repr(value).replace("e", "E")
    return str(value)


def _fixed(value: float, width: int, decimals: int) -> str:
    """Fortran ``Fw.d`` of a single-precision value."""
    return f"{_f32(value):{width}.{decimals}f}"


def _exponent(value: float, width: int = 13, decimals: int = 6) -> str:
    """Fortran ``Ew.d``: ``0.100000E+01`` style (mantissa below 1), right-justified."""
    value = _f32(value)
    if value == 0.0:
        mantissa, exponent = 0.0, 0
    else:
        exponent = math.floor(math.log10(abs(value))) + 1
        mantissa = value / 10.0**exponent
        if abs(round(mantissa, decimals)) >= 1.0:  # rounding carried into a new digit
            mantissa /= 10.0
            exponent += 1
    text = f"{mantissa:.{decimals}f}E{exponent:+03d}"
    return f"{text:>{width}}"


def _list_real(value: float) -> str:
    """A single-precision real as Fortran list-directed output prints it
    (``0.976000011``, ``15.9994001``, ``4.99999989E-02``)."""
    value = _f32(value)
    if value == 0.0:
        return "0.00000000"
    exponent = math.floor(math.log10(abs(value)))
    if -1 <= exponent <= 8:
        return f"{value:.{8 - exponent}f}"
    return f"{value:.8E}"


def _keyword_line(keyword: str, value: str, description: str) -> str:
    return f"{keyword:<{KEYWORD_WIDTH}} {value}{DESCRIPTION_GAP}{description}"


# --------------------------------------------------------------------------- the inputs
@dataclass
class EpsrInputs:
    """The generated inputs of one specification, as texts keyed by file name."""

    spec: SystemSpec
    mol_files: dict[str, str] = field(default_factory=dict)  # 'Water.mol' -> text
    mixato_answers: list[str] = field(default_factory=list)
    nwts_files: dict[str, str] = field(default_factory=dict)  # 'D2O.NWTS.dat' -> text
    inp_name: str = ""  # '<name>.EPSR.inp'
    inp_text: str = ""
    readpdb_answers: list[str] = field(default_factory=list)  # pdb kind only
    changeato_commands: list[str] = field(default_factory=list)  # pdb kind only
    data_files: list[Path] = field(default_factory=list)  # the .mint01 files to copy
    pdb_file: Path | None = None
    #: Other files a setup program reads: ``vanderWaalsRadii.txt`` for ``readpdb``.
    aux_files: dict[str, str] = field(default_factory=dict)

    @property
    def box_name(self) -> str:
        return self.spec.box_name

    @property
    def pdb_name(self) -> str:
        """The PDB file name inside the simulation folder (``<name>.pdb``), so that
        ``readpdb`` writes ``<name>.ato``."""
        return f"{self.spec.name}.pdb"

    @property
    def setup_steps(self) -> list[str]:
        """Human-readable list of the setup programs and arguments, in order."""
        spec = self.spec
        steps: list[str] = []
        if spec.kind == "pdb":
            steps.append(f"readpdb {self.pdb_name} < readpdb.stdin")
            steps.append(f"changeato {self.box_name}.ato < changeato.stdin")
        else:
            steps += [f"makemole {molecule.name}" for molecule in spec.molecules]
            steps.append("mixato < mixato.stdin")
            if spec.schedule.randomise:
                steps.append(f"randomise {self.box_name}")
            steps += [f"delete {molecule.name}.mol" for molecule in spec.molecules]
            if spec.kind == "molecular" and spec.schedule.fmole_iterations > 0:
                steps.append(f"fmole {self.box_name} {spec.schedule.fmole_iterations}")
        steps += [f"nwts {dataset.nwts_file}" for dataset in spec.datasets]
        return steps

    def files(self) -> dict[str, str]:
        """Every generated text keyed by file name (``.mol``, ``.NWTS.dat``, ``.EPSR.inp``,
        and ``vanderWaalsRadii.txt`` for a pdb system)."""
        files = dict(self.mol_files)
        files.update(self.nwts_files)
        files[self.inp_name] = self.inp_text
        files.update(self.aux_files)
        return files

    def answer_files(self) -> dict[str, str]:
        """The standard input of the interactive programs as files (``mixato.stdin``,
        or ``readpdb.stdin`` and ``changeato.stdin``), as in ``validation/``."""
        if self.spec.kind == "pdb":
            return {
                "readpdb.stdin": _lines_text(self.readpdb_answers),
                "changeato.stdin": _lines_text(self.changeato_commands),
            }
        return {"mixato.stdin": _lines_text(self.mixato_answers)}


def _lines_text(lines: list[str]) -> str:
    return "".join(f"{line}\n" for line in lines)


def generate(spec: SystemSpec) -> EpsrInputs:
    """All inputs of ``spec``. Validates the specification first (raises ``SpecError``)."""
    spec.validate()
    _check_names(spec)
    inputs = EpsrInputs(spec=spec)
    if spec.kind == "pdb":
        assert spec.pdb is not None  # validate() guarantees it
        inputs.pdb_file = Path(spec.pdb.file).expanduser().resolve()
        inputs.readpdb_answers = readpdb_answers(spec)
        inputs.changeato_commands = changeato_commands(spec)
        inputs.aux_files[VDW_RADII_FILE] = render_vdw_radii()
        if inputs.pdb_file.is_file():
            _check_pdb_atom_types(inputs.pdb_file, spec)
    else:
        for molecule in spec.molecules:
            inputs.mol_files[f"{molecule.name}.mol"] = render_mol(molecule, spec)
        inputs.mixato_answers = mixato_answers(spec)
    for dataset in spec.datasets:
        inputs.nwts_files[dataset.nwts_file] = render_nwts(dataset, spec)
    inputs.inp_name = spec.inp_name
    inputs.inp_text = render_inp(spec)
    inputs.data_files = [
        Path(dataset.file).expanduser().resolve() for dataset in spec.datasets
    ]
    return inputs


def _check_names(spec: SystemSpec) -> None:
    """File-name clashes inside the simulation folder that ``spec.validate`` does not see."""
    errors = []
    for molecule in spec.molecules:
        if molecule.name == spec.box_name:
            errors.append(
                f"molecule {molecule.name!r} has the run name: mixato would overwrite "
                f"the template {molecule.name}.ato with the box"
            )
        if molecule.n_atoms > MAX_ATOMS_PER_MOLECULE:
            errors.append(
                f"molecule {molecule.name}: {molecule.n_atoms} atoms; makemole labels "
                f"atoms with three digits (at most {MAX_ATOMS_PER_MOLECULE})"
            )
    file_names: dict[str, str] = {}
    for dataset in spec.datasets:
        other = file_names.setdefault(dataset.file_name, dataset.label)
        if other != dataset.label:
            errors.append(
                f"data sets {other} and {dataset.label} use the same file name "
                f"{dataset.file_name}; they would overwrite each other"
            )
    if errors:
        raise SpecError(errors)


def write(inputs: EpsrInputs, folder: str | Path) -> list[Path]:
    """Write the generated files and copy the data (and PDB) files into ``folder``."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    texts = dict(inputs.files())
    texts.update(inputs.answer_files())
    for name, text in texts.items():
        path = folder / name
        path.write_text(text, encoding="latin-1")
        written.append(path)
    for data_file in inputs.data_files:
        if not data_file.is_file():
            raise FileNotFoundError(f"data file {data_file} not found")
        written.append(Path(shutil.copy2(data_file, folder / data_file.name)))
    if inputs.pdb_file is not None:
        if not inputs.pdb_file.is_file():
            raise FileNotFoundError(f"PDB file {inputs.pdb_file} not found")
        written.append(Path(shutil.copy2(inputs.pdb_file, folder / inputs.pdb_name)))
    return written


# --------------------------------------------------------------------------- .mol
def render_mol(molecule: Molecule, spec: SystemSpec) -> str:
    """The ``.mol`` (gmol) text of a molecule template.

    Layout (list-directed, so spacing is free but keep it readable):
    `` .gmol           0``; one ``atom`` line per atom (index, 3-digit label, x, y, z,
    number of neighbours, neighbour indices); ``bond`` lines with the length;
    ``angle`` lines for every :meth:`Molecule.angle_triples`; ``rot`` lines;
    ``changelabel NNN <label>``; ``pseudofac 1.00``; one ``potential`` line per atom
    type used (label, epsilon, sigma, mass, charge, element, nsph, rmin, rmax, cg);
    ``temperature``; ``vibtemp 65``; ``angtemp 1``; ``dihtemp 10``;
    ``ecoredcore 0 0``; ``density 0.001``.
    """
    if molecule.n_atoms > MAX_ATOMS_PER_MOLECULE:
        raise SpecError(
            [
                f"molecule {molecule.name}: {molecule.n_atoms} atoms; makemole labels "
                f"atoms with three digits (at most {MAX_ATOMS_PER_MOLECULE})"
            ]
        )
    neighbours = molecule.neighbours()
    lines = [f" .gmol{0:12d}"]
    for index, atom in enumerate(molecule.atoms, start=1):
        around = neighbours[index]
        line = (
            f" atom{index:13d}  {index:03d}"
            f"{_list_real(atom.x):>14}{_list_real(atom.y):>17}{_list_real(atom.z):>17}"
            f"{len(around):16d}"
        )
        lines.append(line + "".join(f"{other:12d}" for other in around))
    for first, second in molecule.bonds:
        length = molecule.bond_length(first, second)
        lines.append(f"bond {first:03d} {second:03d} {_fixed(length, 10, 5)}")
    for first, vertex, second in molecule.angle_triples():
        angle = molecule.angle_value(first, vertex, second)
        lines.append(
            f"angle {first:03d} {vertex:03d} {second:03d} {_fixed(angle, 10, 5)}"
        )
    for first, second in molecule.rotations:
        lines.append(f"rot {first:3d} {second:3d}")
    for index, atom in enumerate(molecule.atoms, start=1):
        lines.append(f"changelabel {index:03d} {atom.label:<3}")
    lines.append(f"pseudofac {_fixed(1.0, 8, 2)}")
    for label in dict.fromkeys(molecule.labels):
        atom_type = spec.atom_type(label)
        lines.append(
            f"potential {atom_type.label:<3} {_fixed(atom_type.epsilon, 12, 6)} "
            f"{_fixed(atom_type.sigma, 12, 6)} {_fixed(atom_type.mass, 12, 6)} "
            f"{_fixed(atom_type.charge, 12, 6)} {atom_type.element:<3} "
            f"{_exponent(atom_type.nsph)} {_fixed(atom_type.rmin, 12, 6)} "
            f"{_fixed(atom_type.rmax_value, 12, 6)} {_exponent(atom_type.cg)}"
        )
    lines.append(f"temperature {_exponent(spec.temperature)}")
    lines.append(f"vibtemp {_exponent(MOL_VIBTEMP)}")
    lines.append(f"angtemp {_exponent(MOL_ANGTEMP)}")
    lines.append(f"dihtemp {_exponent(MOL_DIHTEMP)}")
    lines.append(f"ecoredcore {_fixed(0.0, 10, 5)} {_fixed(0.0, 10, 5)}")
    lines.append(f"density {_exponent(MOL_DENSITY)}")
    return "".join(f"{line}\n" for line in lines)


def mixato_answers(spec: SystemSpec) -> list[str]:
    """Standard input of ``mixato``: n types; (template, count) per type; density; box name."""
    if spec.kind == "pdb":
        return []
    answers = [str(len(spec.molecules))]
    for molecule in spec.molecules:
        answers += [molecule.name, str(molecule.count)]
    answers += [format_value(float(spec.density)), spec.box_name]
    return answers


# --------------------------------------------------------------------------- .NWTS.dat
def _abundances(dataset: Dataset, label: str) -> str:
    return " ".join(
        f"{isotope.mass_number} {isotope.fraction:g}"
        for isotope in dataset.composition(label)
    )


def render_nwts(dataset: Dataset, spec: SystemSpec) -> str:
    """The ``.NWTS.dat`` text of one data set.

    Global keywords: ``fnameato <box>.ato``, ``output 1``, ``nsamples 1``,
    ``xfrac 0.0``, ``normtot 0``, ``normdif 0``, ``ncomponent <n types>``; one
    ``component`` block per atom type in specification order with ``atom``,
    ``atsymbol`` (element), ``iexchange``, ``abundances`` (``0 1`` natural, ``2 1``
    deuterium, ``1 0.5 2 0.5`` a 50:50 mixture), ``isubs 0``, ``abund2 0 1``.
    """
    values = {
        "fnameato": f"{spec.box_name}.ato",
        "output": 1,
        "nsamples": 1,
        "xfrac": 0.0,
        "normtot": 0,
        "normdif": 0,
        "ncomponent": len(spec.atom_types),
    }
    blocks = [
        {
            "atom": atom_type.label,
            "atsymbol": atom_type.element,
            "iexchange": 1 if atom_type.label in dataset.exchangeable else 0,
            "abundances": _abundances(dataset, atom_type.label),
            "isubs": 0,
            "abund2": "0 1",
        }
        for atom_type in spec.atom_types
    ]
    # write_keyword_file holds the layout EPSR expects; render through a temporary file
    with tempfile.TemporaryDirectory() as temporary:
        path = Path(temporary) / dataset.nwts_file
        epsr_tools.write_keyword_file(
            path, "nwts", values, blocks, title=f"{dataset.label}.NWTS"
        )
        return path.read_text(encoding="latin-1")


# --------------------------------------------------------------------------- .EPSR.inp
def _template_entries() -> (
    tuple[str, list[tuple[str, str, str]], list[tuple[str, str, str]]]
):
    """``(title, global lines, data-block lines)`` of :data:`EPSR_INP_TEMPLATE`, each
    line as ``(keyword, default value, description)``."""
    title = ""
    global_entries: list[tuple[str, str, str]] = []
    block_entries: list[tuple[str, str, str]] = []
    current = global_entries
    for line in EPSR_INP_TEMPLATE.splitlines():
        words = line.split()
        if not words or words == ["q"]:
            continue
        if line.rstrip().endswith(TITLE_MARKER):
            title = words[0]
            continue
        if words[0] == "data" and len(words) == 2 and words[1].isdigit():
            current = block_entries
            continue
        keyword = words[0]
        rest = line.lstrip()[len(keyword) :].lstrip()
        value, _, description = rest.partition(DESCRIPTION_GAP)
        current.append((keyword, value.strip(), description.strip()))
    return title, global_entries, block_entries


def auto_rshmin(spec: SystemSpec) -> float:
    """``rshmin`` when a data set gives none: 0.7 Å for molecular and pdb systems
    (the test-case value, below the shortest bond), 1.2 Å for atomic systems
    (below the shortest expected distance, as for silica); if the shortest bond in
    the molecules is below 0.7 Å, 0.2 Å less than that bond."""
    if spec.kind == "atomic":
        return ATOMIC_RSHMIN
    bonds = [
        molecule.bond_length(first, second)
        for molecule in spec.molecules
        for first, second in molecule.bonds
    ]
    shortest = min(bonds) if bonds else MOLECULAR_RSHMIN
    if shortest < MOLECULAR_RSHMIN:
        return round(shortest - 0.2, 3)
    return MOLECULAR_RSHMIN


def dataset_keywords(dataset: Dataset, spec: SystemSpec) -> dict[str, Any]:
    """The values of one ``data n`` block (``datafile``, ``wtsfile`` and
    :data:`~mellor.spec.DATASET_KEYWORDS`), ``rshmin`` filled in when absent."""
    values: dict[str, Any] = {
        "datafile": dataset.file_name,
        "wtsfile": dataset.weights_file,
    }
    for keyword in DATASET_KEYWORDS:
        values[keyword] = getattr(dataset, keyword)
    if values["rshmin"] is None:
        values["rshmin"] = auto_rshmin(spec)
    return values


def render_inp(spec: SystemSpec, *, quick: bool = False) -> str:
    """The scratch-start ``.EPSR.inp`` text (``ireset 1``, ``iinit 1``, ``nsumt -1``, ``niter 1``).

    ``quick`` is accepted for symmetry with the schedule: the scratch-start file is
    the same for the quick and the full schedule (the iteration counts of the later
    runs are overrides of the script).
    """
    del quick
    _, global_entries, block_entries = _template_entries()
    substitutions: dict[str, Any] = {
        "fnameato": f"{spec.box_name}.ato",
        "fnamepcof": f"{spec.box_name}.pcof",
        "rho": format_value(float(spec.density)),
        "ntimes": spec.schedule.ntimes,
        "ireset": 1,
        "iinit": 1,
        "nsumt": -1,
        "niter": 1,
        "ndata": len(spec.datasets),
    }
    known = {keyword for keyword, _, _ in global_entries}
    unknown = sorted(set(spec.epsr_keywords) - known)
    if unknown:
        raise SpecError([f"epsr_keywords {unknown} are not in the .EPSR.inp template"])
    substitutions.update(spec.epsr_keywords)

    lines = [f"{spec.name}.EPSR{DESCRIPTION_GAP}{TITLE_MARKER}"]
    for keyword, default, description in global_entries:
        value = substitutions.get(keyword, default)
        lines.append(_keyword_line(keyword, format_value(value), description))
    for number, dataset in enumerate(spec.datasets, start=1):
        values = dataset_keywords(dataset, spec)
        lines += ["", f"data {number:3d}", ""]
        for keyword, default, description in block_entries:
            value = values.get(keyword, default)
            lines.append(_keyword_line(keyword, format_value(value), description))
    lines += ["", "q"]
    return "".join(f"{line}\n" for line in lines)


# --------------------------------------------------------------------------- pdb systems
def readpdb_answers(spec: SystemSpec) -> list[str]:
    """Standard input of ``readpdb``: the van der Waals factor and ``y``/``n`` for dihedrals."""
    if spec.pdb is None:
        return []
    return [
        format_value(float(spec.pdb.vdw_factor)),
        "y" if spec.pdb.dihedrals else "n",
    ]


def changeato_commands(spec: SystemSpec) -> list[str]:
    """Standard input of ``changeato``: ``temp``, ``tol``, step sizes, per atom type
    ``atom i`` / ``epsilon`` / ``sigma`` / ``charge`` (and ``atmass`` when it is not the
    element default), extra commands, ``save``, ``y``, ``e``."""
    pdb = spec.pdb
    if pdb is None:
        return []
    commands = [
        f"temp {spec.temperature:g}",
        f"tol {format_value(float(pdb.tolerance))}",
    ]
    ordered = [key for key in CHANGEATO_STEP_KEYS if key in pdb.step_sizes]
    ordered += [key for key in pdb.step_sizes if key not in CHANGEATO_STEP_KEYS]
    commands += [f"{key} {format_value(float(pdb.step_sizes[key]))}" for key in ordered]
    for number, atom_type in enumerate(spec.atom_types, start=1):
        commands += [
            f"atom {number}",
            f"epsilon {format_value(float(atom_type.epsilon))}",
            f"sigma {format_value(float(atom_type.sigma))}",
            f"charge {format_value(float(atom_type.charge))}",
        ]
        element_mass = EPSR_ELEMENT_MASSES.get(atom_type.element)
        if element_mass is None or abs(atom_type.mass - element_mass) > 1e-6:
            commands.append(f"atmass {format_value(float(atom_type.mass))}")
    commands += list(pdb.extra_commands)
    commands += ["save", "y", "e"]
    return commands


def pdb_atom_labels(path: str | Path) -> list[str]:
    """The atom-type labels of a PDB file in the order ``readpdb`` numbers them:
    the atom names (columns 13-16) of the ``ATOM`` records, by first appearance."""
    labels: list[str] = []
    for line in Path(path).read_text(encoding="latin-1").splitlines():
        if line.startswith("ATOM"):
            label = line[12:16].strip()
            if label and label not in labels:
                labels.append(label)
    return labels


def render_vdw_radii() -> str:
    """The ``vanderWaalsRadii.txt`` text ``readpdb`` reads: ``number symbol radius`` per
    element (:data:`VAN_DER_WAALS_RADII`)."""
    return "".join(
        f"{number}   {symbol}   {radius:g}\n"
        for number, (symbol, radius) in enumerate(VAN_DER_WAALS_RADII.items(), start=1)
    )


def _check_pdb_atom_types(path: Path, spec: SystemSpec) -> None:
    """``changeato`` addresses atom types by number, so the specification's atom
    types must be the PDB's, in PDB order."""
    labels = pdb_atom_labels(path)
    if not labels:
        raise SpecError([f"PDB file {path} has no ATOM records"])
    if labels != spec.labels:
        raise SpecError(
            [
                f"atom types {spec.labels} must be the PDB atom names in order of first "
                f"appearance, {labels} (readpdb numbers them so; changeato uses the numbers)"
            ]
        )
