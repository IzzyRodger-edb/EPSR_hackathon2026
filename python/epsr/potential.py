"""Reading and perturbing the EPSR potentials of a live :class:`~epsr.Simulation`.

EPSR's interatomic potential is the sum of a *reference potential*
(Lennard-Jones + Coulomb, from the ``.ato`` file) and an *empirical potential*
(EP) refined from the diffraction data. The EP of each atom-type pair is a sum of
basis functions (Poisson or Gaussian, spacing ``pdstep``) with coefficients
``pcoefp``; the r-space table is rebuilt from the coefficients by the engine
(``generate_ep``). Perturbations therefore act on the coefficients.

Changes are made in memory. They drive the next iteration and are saved to the
``.pcof`` file by the engine at the end of the run.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from .simulation import PairKey, Simulation

__all__ = ["Potential", "PotParams", "TypeParams"]

F32 = npt.NDArray[np.float32]


@dataclass(frozen=True)
class PotParams:
    """Potential parameters (read-only snapshot; set them with ``overrides``)."""

    ereq: float  # absolute_energy: requested EP amplitude (keyword `ereq`)
    potfac: float  # EP refinement factor (0 disables refinement)
    potfacamp: float
    refpotfac: float  # factor on the reference potential
    feedback: float  # confidence factor
    ncoeffp: int  # number of coefficients per pair, including the flag
    pdmin: float
    pdstep: float  # spacing of the basis functions (Å)
    pdmax: float  # range of the EP (Å)
    psigma2: float  # width of the basis functions
    gaussian: bool  # True: Gaussian basis, False: Poisson

    @classmethod
    def from_fortran(cls, d: Mapping[str, Any]) -> PotParams:
        return cls(
            ereq=float(d["absolute_energy"]),
            potfac=float(d["potfac"]),
            potfacamp=float(d["potfacamp"]),
            refpotfac=float(d["refpotfac"]),
            feedback=float(d["feedback"]),
            ncoeffp=int(d["ncoeffp"]),
            pdmin=float(d["pdmin"]),
            pdstep=float(d["pdstep"]),
            pdmax=float(d["pdmax"]),
            psigma2=float(d["psigma2"]),
            gaussian=bool(d["gaussian"]),
        )


@dataclass(frozen=True)
class TypeParams:
    """Reference-potential parameters of one atom type (from the ``.ato`` file)."""

    label: str
    epsilon: float  # Lennard-Jones epsilon (kJ/mol)
    sigma: float  # Lennard-Jones sigma (Å)
    charge: float  # (e)
    mass: float  # (amu)
    fraction: float  # atom fraction of this type


class Potential:
    """Potentials of a :class:`~epsr.Simulation` (``sim.potential``).

    Pair arguments accept an index or two atom-type labels, in either order:
    ``pot.empirical("O", "H")``, ``pot.empirical(("O", "H"))``, ``pot.empirical(3)``.
    Without a pair, arrays are ``(n_pairs, nr)``.
    """

    def __init__(self, sim: Simulation) -> None:
        self._sim = sim

    def __repr__(self) -> str:
        p = self.params
        return f"<epsr.Potential {self._sim.n_pairs} pairs, ereq={p.ereq}, ncoeffp={p.ncoeffp}>"

    def _pick(self, arr: F32, a: str | PairKey | None, b: str | None) -> F32:
        return arr if a is None else arr[self._sim.pair_index(a, b)]

    def pair_index(self, a: str | PairKey, b: str | None = None) -> int:
        return self._sim.pair_index(a, b)

    @property
    def r(self) -> F32:
        """r grid (Å) of the potential tables."""
        return self._sim.r

    # ---------------------------------------------------------------- read
    def empirical(self, a: str | PairKey | None = None, b: str | None = None) -> F32:
        """Empirical potential (kJ/mol, ``xavs``)."""
        return self._pick(self._sim._sim.empirical_potential(), a, b)

    def reference(self, a: str | PairKey | None = None, b: str | None = None) -> F32:
        """Reference potential, Lennard-Jones + Coulomb (kJ/mol, ``ecoul``)."""
        return self._pick(self._sim._sim.reference_potential(), a, b)

    def total(self, a: str | PairKey | None = None, b: str | None = None) -> F32:
        """Total potential used in the Monte Carlo, reference + empirical (``xavtot``)."""
        return self._pick(self._sim._sim.total_potential(), a, b)

    def ep_coefficients(self) -> F32:
        """``(n_pairs, ncoeffp - 1)`` EP coefficients (the refinement-flag column removed)."""
        return self._sim._sim.ep_coefficients()[:, 1:].copy()

    def ep_flags(self) -> F32:
        """``(n_pairs,)`` refinement flags (``pcoefp(1)``) of each pair."""
        return self._sim._sim.ep_coefficients()[:, 0].copy()

    def shared_pairs(self) -> list[list[int]]:
        """Groups of pair indices that share one set of EP coefficients in the ``.pcof`` file."""
        groups: dict[int, list[int]] = {}
        for ip, col in enumerate(self._sim._sim.ep_columns()):
            groups.setdefault(col, []).append(ip)
        return [g for g in groups.values() if len(g) > 1]

    @property
    def params(self) -> PotParams:
        return PotParams.from_fortran(self._sim._sim.pot_params())

    @property
    def type_params(self) -> list[TypeParams]:
        return [TypeParams(**d) for d in self._sim._sim.type_params()]

    # ---------------------------------------------------------------- write
    def _set_full(self, full: npt.ArrayLike) -> None:
        arr = np.ascontiguousarray(full, dtype=np.float32)
        for group in self.shared_pairs():
            rows = arr[group]
            if not np.array_equal(rows, np.broadcast_to(rows[0], rows.shape)):
                labels = [self._sim.pair_labels[i] for i in group]
                raise ValueError(
                    f"pairs {labels} share one set of EP coefficients; give them identical rows "
                    "(perturb_ep updates all of them)"
                )
        self._sim._sim.set_ep_coefficients(arr)

    def set_ep_coefficients(self, coef: npt.ArrayLike) -> None:
        """Replace the EP coefficients (``(n_pairs, ncoeffp - 1)``) and rebuild the EP."""
        coef = np.asarray(coef, dtype=np.float32)
        full = self._sim._sim.ep_coefficients()
        if coef.shape != (full.shape[0], full.shape[1] - 1):
            raise ValueError(f"expected shape {(full.shape[0], full.shape[1] - 1)}, got {coef.shape}")
        full[:, 1:] = coef
        self._set_full(full)

    def set_ep_flags(self, flags: npt.ArrayLike) -> None:
        """Set the per-pair refinement flags (``pcoefp(1)``)."""
        full = self._sim._sim.ep_coefficients()
        full[:, 0] = np.broadcast_to(np.asarray(flags, dtype=np.float32), (full.shape[0],))
        self._set_full(full)

    def perturb_ep(self, pair: str | PairKey, delta: npt.ArrayLike, b: str | None = None) -> None:
        """Add ``delta`` (scalar or ``(ncoeffp - 1,)``) to one pair's EP coefficients.

        Pairs sharing that pair's coefficients (see :meth:`shared_pairs`) change too.
        """
        ip = self._sim.pair_index(pair, b)
        coef = self.ep_coefficients()
        group = next((g for g in self.shared_pairs() if ip in g), [ip])
        coef[group] += np.asarray(delta, dtype=np.float32)
        self.set_ep_coefficients(coef)

    def scale_ep(self, factor: float) -> None:
        """Multiply every EP coefficient by ``factor`` (``multiply_pot_coeffs``)."""
        self._sim._sim.scale_ep(float(factor))

    def zero_ep(self) -> None:
        """Set the empirical potential to zero (like ``ireset 2``)."""
        self.scale_ep(0.0)
