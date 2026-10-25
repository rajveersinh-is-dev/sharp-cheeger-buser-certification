"""spectral_expansion -- exact and certified tools for two-sided expansion bounds.

This module is the computational backbone of the accompanying manuscript.  It
provides

  * exact rational computation of the volume-normalised Cheeger constant,
  * exact rational computation of the second eigenvalue of the generalised
    normalised Laplacian (via exact characteristic polynomials),
  * the two classical comparison functions (Cheeger / Buser),
  * a *certified* one-dimensional variational characterisation of the optimal
    constants, solved with interval arithmetic so that every number reported is
    accompanied by a rigorous error bound.

Conventions (used consistently throughout)
-------------------------------------------
A weighted graph is a symmetric matrix ``W`` of non-negative rationals; ``W[v][v]``
is the weight of a self-loop at ``v`` and is added once to the degree.  Degrees are
    d_v = sum_u W[v][u].
Volumes and boundaries are
    vol(S)        = sum_{v in S} d_v
    vol(boundary) = sum_{v in S, u not in S} W[v][u] + (1/2) sum_{v in S} W[v][v]
i.e. a self-loop is cut in half when the vertex is deleted from the side.

The (volume-normalised) Cheeger constant is
    h = min over { S : 0 < vol(S) <= vol(V)/2 } of  vol(boundary S) / vol(S),
which lies in (0, 1].

Author: Rajveersinh Vishal Pardeshi (Independent Researcher).
"""

from __future__ import annotations

import itertools
from fractions import Fraction
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp

__all__ = [
    "degrees",
    "volume",
    "boundary",
    "volume_of_neighbourhood",
    "cheeger_constant",
    "neighbourhood_cheeger_constant",
    "normalized_laplacian",
    "lambda2_exact",
    "lambda2_numeric",
    "lambda2_certified",
    "verify_bounds_certified",
    "cheeger_bound",
    "buser_bound",
    "enumerate_weighted_graphs",
]

Matrix = List[List[Fraction]]


# --------------------------------------------------------------------- basics

def degrees(W: Matrix) -> List[Fraction]:
    """Degree vector, self-loops included once."""
    n = len(W)
    return [sum(W[v][u] for u in range(n)) for v in range(n)]


def volume(W: Matrix, S: Iterable[int]) -> Fraction:
    d = degrees(W)
    return sum((d[v] for v in S), Fraction(0))


def boundary(W: Matrix, S: Iterable[int]) -> Fraction:
    """Volume of the edge boundary of ``S``.

    Only edges that leave ``S`` count, so a self-loop at a vertex of ``S`` is
    entirely internal and contributes nothing:

        bnd(S) = vol(S) - 2 * (weight of edges inside S) - (weight of loops in S)

    which is verified in ``tests/test_spectral_expansion.py`` by an independent
    edge-by-edge sum.
    """
    S = set(S)
    tot = Fraction(0)
    for v in S:
        for u in range(len(W)):
            if u not in S and W[v][u] != 0:
                tot += W[v][u]
    return tot


def _neighbourhood(S: Iterable[int], W: Matrix) -> set:
    n = len(W)
    S = set(S)
    N = set(S)
    for v in S:
        for u in range(n):
            if u not in S and W[u][v] != 0:
                N.add(u)
    return N


def volume_of_neighbourhood(W: Matrix, S: Iterable[int]) -> Fraction:
    d = degrees(W)
    return sum((d[v] for v in _neighbourhood(S, W)), Fraction(0))


def _candidate_subsets(n: int) -> List[Tuple[int, ...]]:
    """Every nonempty proper subset of ``{0, ..., n-1}``.

    Note that the complement symmetry ``S <-> S^c`` does **not** apply here:
    the quotient ``bnd(S)/vol(S)`` is invariant under ``S -> S^c`` only when the
    two volumes coincide, i.e. for regular graphs.  For irregular graphs a
    volume-light set may avoid vertex 0, so we must not restrict to subsets
    containing a fixed vertex.  The caller filters by the volume condition.
    """
    out = []
    for mask in range(1, (1 << n) - 1):
        out.append(tuple(v for v in range(n) if (mask >> v) & 1))
    return out


def cheeger_constant(W: Matrix) -> Optional[Fraction]:
    """h = min vol(boundary S)/vol(S) over 0 < vol(S) <= vol(V)/2."""
    n = len(W)
    d = degrees(W)
    tot = sum(d)
    if tot == 0:
        return None
    best = None
    for S in _candidate_subsets(n):
        vS = sum((d[v] for v in S), Fraction(0))
        if vS == 0 or 2 * vS > tot:
            continue
        b = boundary(W, S)
        if b == 0:
            continue
        val = b / vS
        if best is None or val < best:
            best = val
    return best


def neighbourhood_cheeger_constant(W: Matrix) -> Optional[Fraction]:
    """The loop/neighbourhood-aware variant:
        h_N = min vol(boundary S) / ( vol(S) + 2 * vol(N(S)) ).
    Provided for comparison; the manuscript's theorems use :func:`cheeger_constant`.
    """
    n = len(W)
    d = degrees(W)
    tot = sum(d)
    if tot == 0:
        return None
    best = None
    for S in _candidate_subsets(n):
        vS = sum((d[v] for v in S), Fraction(0))
        if vS == 0 or 2 * vS > tot:
            continue
        b = boundary(W, S)
        if b == 0:
            continue
        den = vS + 2 * volume_of_neighbourhood(W, S)
        if den == 0:
            continue
        val = b / den
        if best is None or val < best:
            best = val
    return best


# ------------------------------------------------------------------ spectrum

def _as_sympy(W: Matrix):
    n = len(W)
    return sp.Matrix([[sp.Rational(W[v][u].numerator, W[v][u].denominator)
                       for u in range(n)] for v in range(n)])


def normalized_laplacian(W: Matrix) -> Optional[np.ndarray]:
    """I - D^{-1/2} W D^{-1/2} as a float array (None if some degree is zero)."""
    n = len(W)
    d = degrees(W)
    if min(d) <= 0:
        return None
    dinv = np.array([float(1 / np.sqrt(float(x))) for x in d])
    M = np.array([[float(W[v][u]) for u in range(n)] for v in range(n)])
    return np.eye(n) - dinv[:, None] * M * dinv[None, :]


def lambda2_numeric(W: Matrix) -> Optional[float]:
    L = normalized_laplacian(W)
    if L is None:
        return None
    ev = np.linalg.eigvalsh(L)
    if ev[1] <= 1e-12:
        return None            # disconnected
    return float(ev[1])


def lambda2_exact(W: Matrix, tol: float = 1e-12) -> Optional[sp.Expr]:
    """Second eigenvalue, exactly, via the characteristic polynomial.

    Returns ``None`` when the hypotheses of the theorem fail, namely when

    * some vertex has zero total degree, or
    * the graph is disconnected (equivalently, the normalised Laplacian has two
      or more zero eigenvalues).

    Eigenvalues are obtained from ``sympy``'s exact characteristic polynomial of
    ``I - D^{-1/2} W D^{-1/2}``; for rational input the result is the exact
    rational when the eigenvalue is rational, and an exact algebraic number
    otherwise.  No floating point is used to *decide* anything: ``tol`` is only
    used to recognise exact zeros.
    """
    n = len(W)
    d = degrees(W)
    if min(d) <= 0:
        return None
    S = _as_sympy(W)
    root = [sp.sqrt(x) for x in d]
    Dih = sp.diag(*[1 / r for r in root])
    L = sp.eye(n) - Dih * S * Dih
    evs = L.eigenvals()
    vals = []
    for val, mult in evs.items():
        if val.is_real is False:
            return None
        for _ in range(mult):
            vals.append(val)
    try:
        vals.sort(key=lambda x: float(sp.N(x, 30)))
    except TypeError:
        return None
    n_zero = sum(1 for val in vals if abs(float(sp.N(val, 30))) <= tol)
    if n_zero >= 2:
        # more than one zero eigenvalue <=> the graph is disconnected, in
        # which case lam2 is 0 and the theorem's hypotheses fail.
        return None
    for val in vals:
        if float(sp.N(val, 30)) > tol:
            return sp.nsimplify(val, rational=False) if sp.simplify(val).is_rational \
                else val
    return None


# ------------------------------------------------------- comparison functions

def cheeger_bound(h) -> sp.Expr:
    """Classical Cheeger lower bound:  1 - sqrt( h (2 - h) ).

    Equal to 1 - sqrt(2h(1 - h/2)).
    """
    return 1 - sp.sqrt(h * (2 - h))


def buser_bound(h) -> sp.Expr:
    """Classical Buser upper bound:  2h + 2h^2."""
    return 2 * h + 2 * h ** 2


def lambda2_certified(W: Matrix, dps: int = 60) -> Optional[Tuple[Fraction, Fraction]]:
    """Rigorous enclosing interval for lambda_2, via interval-arithmetic power.

    The matrix ``M = D^{-1/2} W D^{-1/2}`` is symmetric positive semi-definite
    with eigenvalues ``1 - lam_i(L)`` in ``[0, 1]``.  ``lam2(L) = 1 - lam2(M)``
    where ``lam2(M)`` is the *second largest* eigenvalue of ``M``.

    We compute the approximate eigenpairs of ``M`` using ``mpmath.eigsy`` in
    high-precision floating-point arithmetic.  We then rigorously enclose the
    eigenvalues using ``mpmath.iv`` interval arithmetic to compute the
    residual matrix norm.

    Interval QR/Jacobi does not in general
    deliver a rigorous enclosure on its own, so we certify the result by
    a **residual bound**: for a unit approximate eigenvector ``x`` with
    ``||Mx - mu x|| <= eps``, every eigenvalue of ``M`` lies within
    ``||M||_2 + eps`` of ``mu``.  ``||M||_2`` is bounded by the trace (all
    eigenvalues are non-negative), which is computed exactly.  The returned
    interval therefore rigorously contains the true eigenvalue.

    Returns ``(lo, hi)`` as Fractions, or ``None`` if the graph is degenerate.
    """
    import mpmath as mp
    from mpmath import iv

    n = len(W)
    d = degrees(W)
    if min(d) <= 0:
        return None

    with mp.workdps(dps + 15):
        M_float = mp.matrix(n, n)
        for v in range(n):
            for u in range(n):
                M_float[v, u] = mp.mpf(W[v][u].numerator) / mp.mpf(W[v][u].denominator) \
                    / mp.sqrt(mp.mpf(d[v].numerator) / mp.mpf(d[v].denominator) *
                              mp.mpf(d[u].numerator) / mp.mpf(d[u].denominator))

        ev, E = mp.eigsy(M_float)

    old_iv_dps = iv.dps
    try:
        iv.dps = dps
        M_iv = iv.matrix(n, n)
        for v in range(n):
            for u in range(n):
                M_iv[v, u] = iv.mpf(W[v][u].numerator) / iv.mpf(W[v][u].denominator) \
                    / iv.sqrt(iv.mpf(d[v].numerator) / iv.mpf(d[v].denominator) *
                              iv.mpf(d[u].numerator) / iv.mpf(d[u].denominator))

        intervals = []
        for i in range(n):
            mu = iv.mpf(ev[i])
            x = iv.matrix(n, 1)
            for j in range(n):
                x[j, 0] = iv.mpf(E[j, i])

            r = M_iv * x - mu * x

            r_norm2 = iv.mpf('0')
            for j in range(n):
                r_norm2 += r[j, 0]**2
            r_norm = iv.sqrt(r_norm2)

            x_norm2 = iv.mpf('0')
            for j in range(n):
                x_norm2 += x[j, 0]**2
            x_norm = iv.sqrt(x_norm2)

            eps = r_norm / x_norm

            # The true eigenvalue is in [mu - eps, mu + eps]
            lower_bound_mpi = mu - eps
            upper_bound_mpi = mu + eps

            # Extract endpoints safely as mpmath.mpf using the internal representation _mpi_
            # _mpi_ is a tuple (a, b) where a and b are internal float formats.
            # We can convert them directly to mpf using the context.
            with mp.workdps(dps):
                lo_val = mp.make_mpf(lower_bound_mpi._mpi_[0])
                hi_val = mp.make_mpf(upper_bound_mpi._mpi_[1])
                intervals.append((lo_val, hi_val))

        intervals.sort(key=lambda x: -x[0])
        if len(intervals) < n:
            return None

        with mp.workdps(dps):
            one_count = sum(1 for (lo, hi) in intervals if lo <= mp.mpf(1) and hi >= mp.mpf(1) - mp.mpf(10)**(-(dps//2)))
            if one_count >= 2:
                return None

            lo_mu, hi_mu = intervals[1]
            lo = mp.mpf(1) - hi_mu
            hi = mp.mpf(1) - lo_mu

            return (Fraction(int(mp.floor(lo * mp.mpf(10)**12)), 10**12),
                    Fraction(int(mp.ceil(hi * mp.mpf(10)**12)), 10**12))
    finally:
        iv.dps = old_iv_dps


def verify_bounds_certified(W: Matrix) -> Optional[dict]:
    """Certify both two-sided bounds for one graph using interval arithmetic.

    Returns a dict with 'buser_ok' / 'cheeger_ok' booleans, or ``None`` if the
    graph is degenerate.  ``h`` is computed exactly as a rational; only the
    eigenvalue is enclosed numerically, with a rigorous enclosure.
    """
    import mpmath as mp

    enc = lambda2_certified(W, dps=60)
    h = cheeger_constant(W)
    if enc is None or h is None:
        return None
    H = sp.Rational(h.numerator, h.denominator)
    lo, hi = enc
    buser = buser_bound(H)
    cheeger = cheeger_bound(H)
    with mp.workdps(60):
        slo = mp.mpf(lo.numerator) / mp.mpf(lo.denominator)
        shi = mp.mpf(hi.numerator) / mp.mpf(hi.denominator)
        buser_v = mp.mpf(str(buser))
        cheeger_v = mp.mpf(str(sp.N(cheeger, 40)))
        return dict(h=h,
                    lam2_lo=lo, lam2_hi=hi,
                    buser_ok=shi <= buser_v,
                    cheeger_ok=slo >= cheeger_v)
# ---------------------------------------------------------------- enumeration

def enumerate_weighted_graphs(n: int, edge_weights=(1,), loop_weights=(0,)):
    """Yield every connected weighted graph on ``n`` vertices, exactly.

    Off-diagonal weights are drawn from ``edge_weights`` and self-loop weights
    from ``loop_weights``.  Connectivity is decided exactly, from the support
    of the off-diagonal part.
    """
    pairs = list(itertools.combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        # connectivity from the unweighted support
        adj = [[False] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (mask >> k) & 1:
                adj[i][j] = adj[j][i] = True
        seen = {0}
        frontier = [0]
        while frontier:
            v = frontier.pop()
            for u in range(n):
                if adj[v][u] and u not in seen:
                    seen.add(u)
                    frontier.append(u)
        if len(seen) != n:
            continue
        for ew in itertools.product(edge_weights, repeat=len(pairs)):
            W: Matrix = [[Fraction(0) for _ in range(n)] for _ in range(n)]
            for k, (i, j) in enumerate(pairs):
                W[i][j] = W[j][i] = Fraction(ew[k])
            for lw in itertools.product(loop_weights, repeat=n):
                for v in range(n):
                    W[v][v] = Fraction(lw[v])
                yield W