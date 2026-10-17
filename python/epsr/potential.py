"""Reading and perturbing the EPSR potentials of a live :class:`~epsr.Simulation`.

EPSR's interatomic potential is the sum of a *reference potential*
(Lennard-Jones + Coulomb, from the ``.ato`` file) and an *empirical potential*
(EP) refined from the diffraction data. The EP of each atom-type pair is a sum of
basis functions (Poisson or Gaussian, spacing ``pdstep``) with coefficients
``pcoefp``; the r-space table is rebuilt from the coefficients by the engine
(``generate_ep``). Perturbations therefore act on the coefficients.

Changes are made in memory. They drive the next iteration and are saved to the
``.pcof`` file by the engine at the end of the run.

**The ereq amplitude limit.** Every time the EP is rebuilt, the engine measures
its absolute energy and rescales the coefficients so that it does not exceed the
requested amplitude ``ereq`` (``generate_ep``; with ``ereq = 0`` the EP is set to
zero). Scaling an EP down, zeroing it or perturbing it within the limit takes
effect as given; scaling it beyond ``ereq`` is undone by the engine. Raise the
limit with ``overrides={"ereq": ...}`` when the Simulation is created.
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

Float32Array = npt.NDArray[np.float32]

#: Column of the engine's coefficient table that holds the per-pair refinement flag.
FLAG_COLUMN = 0


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
    def from_fortran(cls, fortran_values: Mapping[str, Any]) -> PotParams:
        """Build from the engine's dict (Fortran variable names)."""
        return cls(
            ereq=float(fortran_values["absolute_energy"]),
            potfac=float(fortran_values["potfac"]),
            potfacamp=float(fortran_values["potfacamp"]),
            refpotfac=float(fortran_values["refpotfac"]),
            feedback=float(fortran_values["feedback"]),
            ncoeffp=int(fortran_values["ncoeffp"]),
            pdmin=float(fortran_values["pdmin"]),
            pdstep=float(fortran_values["pdstep"]),
            pdmax=float(fortran_values["pdmax"]),
            psigma2=float(fortran_values["psigma2"]),
            gaussian=bool(fortran_values["gaussian"]),
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
    ``pot.empirical("O1", "H1")``, ``pot.empirical(("O1", "H1"))``, ``pot.empirical(3)``.
    Without a pair, arrays are ``(n_pairs, nr)``.
    """

    def __init__(self, simulation: Simulation) -> None:
        self._simulation = simulation
        # The Rust engine object, which holds the potential tables and coefficients.
        self._engine = simulation._engine

    def __repr__(self) -> str:
        params = self.params
        return (
            f"<epsr.Potential {self._simulation.n_pairs} pairs, "
            f"ereq={params.ereq}, ncoeffp={params.ncoeffp}>"
        )

    def _select_pair(
        self, all_pairs: Float32Array, pair: str | PairKey | None, other_type: str | None
    ) -> Float32Array:
        """All rows, or the row of one pair."""
        if pair is None:
            return all_pairs
        return all_pairs[self._simulation.pair_index(pair, other_type)]

    def pair_index(self, pair_or_type: str | PairKey, other_type: str | None = None) -> int:
        """Index of a pair, as :meth:`Simulation.pair_index`."""
        return self._simulation.pair_index(pair_or_type, other_type)

    @property
    def r(self) -> Float32Array:
        """r grid (Å) of the potential tables."""
        return self._simulation.r

    # ---------------------------------------------------------------- read
    def empirical(
        self, pair: str | PairKey | None = None, other_type: str | None = None
    ) -> Float32Array:
        """Empirical potential (kJ/mol, ``xavs``)."""
        return self._select_pair(self._engine.empirical_potential(), pair, other_type)

    def reference(
        self, pair: str | PairKey | None = None, other_type: str | None = None
    ) -> Float32Array:
        """Reference potential, Lennard-Jones + Coulomb (kJ/mol, ``ecoul``)."""
        return self._select_pair(self._engine.reference_potential(), pair, other_type)

    def total(
        self, pair: str | PairKey | None = None, other_type: str | None = None
    ) -> Float32Array:
        """Total potential used in the Monte Carlo, reference + empirical (``xavtot``)."""
        return self._select_pair(self._engine.total_potential(), pair, other_type)

    def ep_coefficients(self) -> Float32Array:
        """``(n_pairs, ncoeffp - 1)`` EP coefficients (the refinement-flag column removed)."""
        coefficients_with_flags = self._engine.ep_coefficients()
        return coefficients_with_flags[:, FLAG_COLUMN + 1 :].copy()

    def ep_flags(self) -> Float32Array:
        """``(n_pairs,)`` refinement flags (``pcoefp(1)``) of each pair."""
        coefficients_with_flags = self._engine.ep_coefficients()
        return coefficients_with_flags[:, FLAG_COLUMN].copy()

    def shared_pairs(self) -> list[list[int]]:
        """Groups of pair indices that share one set of EP coefficients in the ``.pcof`` file."""
        pairs_by_pcof_column: dict[int, list[int]] = {}
        for pair_index, pcof_column in enumerate(self._engine.ep_columns()):
            pairs_by_pcof_column.setdefault(pcof_column, []).append(pair_index)
        return [pairs for pairs in pairs_by_pcof_column.values() if len(pairs) > 1]

    @property
    def params(self) -> PotParams:
        return PotParams.from_fortran(self._engine.pot_params())

    @property
    def type_params(self) -> list[TypeParams]:
        return [TypeParams(**values) for values in self._engine.type_params()]

    # ---------------------------------------------------------------- write
    def _write_coefficients_with_flags(self, coefficients_with_flags: npt.ArrayLike) -> None:
        """Send a full ``(n_pairs, ncoeffp)`` table to the engine, which rebuilds the EP."""
        table = np.ascontiguousarray(coefficients_with_flags, dtype=np.float32)
        for shared_group in self.shared_pairs():
            rows = table[shared_group]
            if not np.array_equal(rows, np.broadcast_to(rows[0], rows.shape)):
                labels = [self._simulation.pair_labels[index] for index in shared_group]
                raise ValueError(
                    f"pairs {labels} share one set of EP coefficients; give them identical rows "
                    "(perturb_ep updates all of them)"
                )
        self._engine.set_ep_coefficients(table)

    def set_ep_coefficients(self, coefficients: npt.ArrayLike) -> None:
        """Replace the EP coefficients (``(n_pairs, ncoeffp - 1)``) and rebuild the EP.

        If the new EP exceeds ``ereq`` the engine scales all coefficients down to
        the limit; read them back with :meth:`ep_coefficients`.
        """
        coefficients = np.asarray(coefficients, dtype=np.float32)
        table = self._engine.ep_coefficients()
        n_pairs, n_columns = table.shape
        expected_shape = (n_pairs, n_columns - 1)
        if coefficients.shape != expected_shape:
            raise ValueError(f"expected shape {expected_shape}, got {coefficients.shape}")
        table[:, FLAG_COLUMN + 1 :] = coefficients
        self._write_coefficients_with_flags(table)

    def set_ep_flags(self, flags: npt.ArrayLike) -> None:
        """Set the per-pair refinement flags (``pcoefp(1)``): one value, or one per pair."""
        table = self._engine.ep_coefficients()
        n_pairs = table.shape[0]
        table[:, FLAG_COLUMN] = np.broadcast_to(np.asarray(flags, dtype=np.float32), (n_pairs,))
        self._write_coefficients_with_flags(table)

    def perturb_ep(self, pair: PairKey, delta: npt.ArrayLike) -> None:
        """Add ``delta`` (scalar or ``(ncoeffp - 1,)``) to one pair's EP coefficients.

        ``pair`` is an index or a label tuple such as ``("O1", "H1")``. Pairs
        sharing that pair's coefficients (see :meth:`shared_pairs`) change too.
        """
        pair_index = self._simulation.pair_index(pair)
        coefficients = self.ep_coefficients()
        affected_pairs = next(
            (group for group in self.shared_pairs() if pair_index in group), [pair_index]
        )
        coefficients[affected_pairs] += np.asarray(delta, dtype=np.float32)
        self.set_ep_coefficients(coefficients)

    def scale_ep(self, factor: float) -> None:
        """Multiply every EP coefficient by ``factor`` (``multiply_pot_coeffs``) and rebuild the EP.

        The engine then re-applies the ``ereq`` amplitude limit (see the module
        documentation), so the EP never grows beyond ``ereq``.
        """
        self._engine.scale_ep(float(factor))

    def zero_ep(self) -> None:
        """Set the empirical potential to zero (like ``ireset 2``)."""
        self.scale_ep(0.0)
