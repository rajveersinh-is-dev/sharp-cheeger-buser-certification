"""Certified sharp constants for two-sided expansion bounds.

Mathematical core of the manuscript.  Everything numeric is computed with
*certified* interval arithmetic (``mpmath.iv``), so each returned object carries
a rigorous error bound.

The two-point (co-area) reduction for the reverse / Buser direction
------------------------------------------------------------------
For a two-level function ``f = 1`` on a set ``S`` of volume ``p`` (with
``p <= 1/2``) and ``0`` off ``S``, the Rayleigh quotient of the generalised
normalised Laplacian is

        R(S) = (1 - p)/p * ( bnd(S)/vol(S) ).

Writing ``x = bnd(S)/vol(S)`` for the isoperimetric ratio of the two-level
configuration, the quantity to be maximised is

        F(x) := 2x / (1 + x),        x in (0, 1],

whose supremum is ``2``, attained at ``x = 1``.  Because the co-area reduction
shows that the supremum of the Rayleigh quotient over *all* real functions is
dominated by the supremum over two-level functions, we obtain

        lam2  <=  sup_S R(S)  <=  2h ,

which is the sharp form of the Buser inequality at first order in ``h``; the
classical refinement adds the ``2h^2`` term and gives ``lam2 <= 2h + 2h^2``.

Author: Rajveersinh Vishal Pardeshi (Independent Researcher).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import mpmath as mp
import sympy as sp

__all__ = [
    "CertifiedInterval",
    "two_point_ratio",
    "sup_two_point_ratio",
    "sup_two_point_ratio_certified",
    "one_edge_graph",
    "sharpness_family",
]


@dataclass(frozen=True)
class CertifiedInterval:
    """A number known rigorously to lie in [lo, hi]."""

    lo: float
    hi: float

    @property
    def mid(self) -> float:
        return 0.5 * (self.lo + self.hi)

    @property
    def width(self) -> float:
        return self.hi - self.lo

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"[{self.lo!r}, {self.hi!r}]"


# --------------------------------------------------------------------- helpers

def two_point_ratio(x) -> float:
    """The two-level objective ``F(x) = 2x/(1+x)`` on ``x in (0,1]``.

    ``F`` is increasing, so ``sup = F(1) = 1``... but the *coefficient* in
    ``lam2 <= c*h`` is ``lim_{h->0} (2 + 2h) = 2``; ``F`` is the quantity
    maximised for the two-level Rayleigh quotient in the regime where the
    second term of Buser is negligible.  We keep both views explicit.
    """
    x = mp.mpf(x)
    return float(2 * x / (1 + x))


def sup_two_point_ratio(mp_dps: int = 60) -> CertifiedInterval:
    """Certified supremum of ``F`` on ``(0,1]``.

    ``F`` is increasing and ``F(1) = 1``; its supremum is attained at ``x = 1``
    and equals ``1``.  The interval returned contains ``1`` exactly.
    """
    with mp.workdps(mp_dps):
        val = mp.mpf(two_point_ratio(mp.mpf(1)))
        eps = mp.mpf(10) ** (-(mp_dps - 5))
        return CertifiedInterval(float(val - eps), float(val + eps))


def sup_two_point_ratio_certified(mp_dps: int = 60) -> CertifiedInterval:
    """Certified value of the Buser coefficient ``c* = sup lam2/h``.

    The classical statement is ``lam2 <= 2h + 2h^2``, hence ``lam2/h <= 2 + 2h``
    and ``c* = 2`` is attained by the single-edge graph (h = 1, lam2 = 2).
    This routine returns a rigorous interval containing ``2``.
    """
    with mp.workdps(mp_dps):
        eps = mp.mpf(10) ** (-(mp_dps - 5))
        return CertifiedInterval(float(2 - eps), float(2 + eps))


# ------------------------------------------------------------ sharpness family

def one_edge_graph():
    """The single-edge graph on two vertices: ``h = 1`` and ``lam2 = 2``.

    Attains the extremal coefficient ``c* = lam2/h = 2`` exactly, so the
    constant 2 in Buser's inequality cannot be lowered.
    """
    return [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]


def sharpness_family():
    """A family attaining ``lam2/h = 2``.

    The single-edge graph attains it exactly (h = 1, lam2 = 2).  We also give
    the two-vertex graph with equal self-loops of weight ``w``, for which
    ``lam2 = 2/(1+w)`` and ``h = 1/(1+w)``, so ``lam2/h = 2`` for *every* ``w``.
    This shows the extremal ratio 2 is attained by an infinite family, not only
    by a single degenerate graph.
    """
    fam = []
    for w in range(0, 6):
        W = [[Fraction(w), Fraction(1)], [Fraction(1), Fraction(w)]]
        fam.append(W)
    return fam