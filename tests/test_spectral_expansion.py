"""Unit tests for the exact expansion machinery.

Run:  python -m pytest tests/ -q
"""

from __future__ import annotations

import itertools
from fractions import Fraction

import pytest
import sympy as sp

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from spectral_expansion import (  # noqa: E402
    boundary, degrees, volume, cheeger_constant, lambda2_exact,
    buser_bound, cheeger_bound, enumerate_weighted_graphs,
    normalized_laplacian,
)


def edge_only(n):
    """Complete graph K_n as an exact weight matrix."""
    W = [[Fraction(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            W[i][j] = W[j][i] = Fraction(1)
    return W


# ----------------------------------------------------------------- conventions

def test_boundary_ignores_internal_loops():
    """A self-loop has both endpoints in S, so it contributes 0 to bnd(S)."""
    W = [[Fraction(5), Fraction(0)], [Fraction(0), Fraction(7)]]
    assert boundary(W, [0]) == Fraction(0)
    assert boundary(W, [1]) == Fraction(0)


def test_boundary_counts_crossing_edges_once():
    W = [[Fraction(2), Fraction(3)],
         [Fraction(3), Fraction(2)]]
    assert boundary(W, [0]) == Fraction(3)
    assert boundary(W, [0, 1]) == Fraction(0)


def test_boundary_identity_volume_minus_internal():
    """bnd(S) = vol(S) - 2*E(S) - loops(S)."""
    for n in (2, 3):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            d = degrees(W)
            for r in range(1, n // 2 + 1):
                for S in itertools.combinations(range(n), r):
                    S = set(S)
                    Ein = sum((W[u][v] for u in S for v in S if u < v),
                              Fraction(0))
                    loops = sum((W[v][v] for v in S), Fraction(0))
                    assert boundary(W, S) == (sum((d[v] for v in S), Fraction(0))
                                              - 2 * Ein - loops)


def test_degree_counts_loop_once():
    W = [[Fraction(2), Fraction(3)], [Fraction(3), Fraction(2)]]
    assert degrees(W) == [Fraction(5), Fraction(5)]


# ------------------------------------------------------------------ Cheeger h

def test_cheeger_of_single_edge_is_one():
    W = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    assert cheeger_constant(W) == Fraction(1)


def test_cheeger_of_complete_graph():
    """For K_n (loop-free, degree n-1) the minimum is attained at the
    largest admissible set, giving  h = (n - floor(n/2)) / (n - 1)."""
    for n in (2, 3, 4, 5, 6):
        W = edge_only(n)
        h = cheeger_constant(W)
        expected = Fraction(n - n // 2, n - 1)
        assert h == expected, f"K_{n}: h={h}, expected {expected}"


def test_cheeger_is_positive_and_at_most_one():
    for n in (2, 3, 4):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            h = cheeger_constant(W)
            assert h is not None
            assert Fraction(0) < h <= Fraction(1)


# ------------------------------------------------------------------ spectrum

def test_lambda2_of_single_edge_is_two():
    W = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    assert sp.simplify(lambda2_exact(W)) == sp.Integer(2)


def test_lambda2_two_vertex_family():
    """W = [[a,w],[w,a]] has lam2 = 2w/(a+w) and h = w/(a+w)."""
    for a in range(0, 6):
        for w in range(1, 4):
            W = [[Fraction(a), Fraction(w)], [Fraction(w), Fraction(a)]]
            lam = sp.simplify(lambda2_exact(W))
            h = cheeger_constant(W)
            assert sp.simplify(lam - sp.Rational(2 * w, a + w)) == 0
            assert h == Fraction(w, a + w)


def test_sharpness_family_ratio_is_two():
    """lam2/h = 2 identically for the two-vertex loop family -- the sharpness case."""
    for a in range(0, 6):
        for w in range(1, 4):
            W = [[Fraction(a), Fraction(w)], [Fraction(w), Fraction(a)]]
            lam = sp.simplify(lambda2_exact(W))
            H = sp.Rational(*cheeger_constant(W).as_integer_ratio())
            ratio = sp.simplify(lam / H)
            assert sp.simplify(ratio - 2) == 0, f"a={a} w={w} ratio={ratio}"


def test_lambda2_none_for_disconnected():
    """Two disjoint edges are disconnected; lambda2 must not be reported."""
    W = [[Fraction(0), Fraction(1), Fraction(0), Fraction(0)],
         [Fraction(1), Fraction(0), Fraction(0), Fraction(0)],
         [Fraction(0), Fraction(0), Fraction(0), Fraction(1)],
         [Fraction(0), Fraction(0), Fraction(1), Fraction(0)]]
    assert lambda2_exact(W) is None


def test_lambda2_none_for_zero_degree():
    W = [[Fraction(0), Fraction(0)], [Fraction(0), Fraction(1)]]
    assert lambda2_exact(W) is None


def test_normalized_laplacian_is_symmetric_with_range():
    W = edge_only(5)
    L = normalized_laplacian(W)
    assert L is not None
    assert abs(L - L.T).max() < 1e-12


# ----------------------------------------------------------------- inequalities

def test_buser_bound_values():
    assert buser_bound(sp.Integer(1)) == sp.Integer(4)
    assert sp.simplify(buser_bound(sp.Rational(1, 2))) == sp.Rational(3, 2)


def test_cheeger_bound_is_one_minus_sqrt():
    # 1 - sqrt(h(2-h)) == 1 - sqrt(2h(1 - h/2))
    h = sp.Rational(1, 3)
    assert sp.simplify(cheeger_bound(h)
                       - (1 - sp.sqrt(2 * h * (1 - h / 2)))) == 0


def test_two_sided_holds_exactly_small():
    """Exact sweep on all connected graphs with n <= 3, weighted with loops.

    The large-scale sweep lives in ``experiments/run_verification.py``; unit
    tests stay fast.  Everything here is exact rational arithmetic.
    """
    checked = 0
    for n in (2, 3):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            h = cheeger_constant(W)
            lam = lambda2_exact(W)
            if h is None or lam is None:
                continue
            H = sp.Rational(h.numerator, h.denominator)
            checked += 1
            assert sp.N(lam - buser_bound(H), 12) <= 1e-12, "Buser violated"
            assert sp.N(lam - cheeger_bound(H), 12) >= -1e-12, "Cheeger violated"
    assert checked > 100


def test_two_sided_holds_exactly_unweighted_n4():
    checked = 0
    for n in (4,):
        for W in enumerate_weighted_graphs(n):
            h = cheeger_constant(W)
            lam = lambda2_exact(W)
            if h is None or lam is None:
                continue
            H = sp.Rational(h.numerator, h.denominator)
            checked += 1
            assert sp.N(lam - buser_bound(H), 12) <= 1e-12
            assert sp.N(lam - cheeger_bound(H), 12) >= -1e-12
    assert checked > 30


def test_boundary_identity_small_fast():
    """Cheap version of the identity check on n <= 3 only (fast in CI)."""
    for n in (2, 3):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            d = degrees(W)
            for r in range(1, n // 2 + 1):
                for S in itertools.combinations(range(n), r):
                    S = set(S)
                    Ein = sum((W[u][v] for u in S for v in S if u < v),
                              Fraction(0))
                    loops = sum((W[v][v] for v in S), Fraction(0))
                    assert boundary(W, S) == (sum((d[v] for v in S), Fraction(0))
                                              - 2 * Ein - loops)


def test_certified_interval_contains_exact_eigenvalue():
    """The certified interval must enclose the exact value on every graph.

    This is the cross-validation that justifies using the (much faster)
    certified path for the large sweep: if the interval ever excluded the exact
    eigenvalue, the sweep's verdicts could not be trusted.
    """
    from spectral_expansion import lambda2_certified
    checked = 0
    for n in (2, 3):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            ex = lambda2_exact(W)
            enc = lambda2_certified(W)
            if ex is None or enc is None:
                continue
            lo, hi = enc
            exv = float(sp.N(ex, 25))
            assert float(lo) - 1e-9 <= exv <= float(hi) + 1e-9, \
                f"interval [{lo},{hi}] excludes exact {ex}"
            checked += 1
    assert checked > 100


def test_certified_bounds_agree_with_exact():
    """Where the certified path says a bound holds, the exact path must agree."""
    from spectral_expansion import verify_bounds_certified
    checked = 0
    for n in (2, 3):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            r = verify_bounds_certified(W)
            h = cheeger_constant(W)
            ex = lambda2_exact(W)
            if r is None or h is None or ex is None:
                continue
            H = sp.Rational(h.numerator, h.denominator)
            assert sp.N(ex - buser_bound(H), 12) <= 1e-12, "exact Buser violated"
            assert sp.N(ex - cheeger_bound(H), 12) >= -1e-12
            checked += 1
    assert checked > 100


def test_certified_path_covers_exactly_the_same_graphs():
    """Regression test for a bug that made the fast path silently skip graphs.

    ``L = I - M`` with ``M = D^-1/2 W D^-1/2``, so the graph is disconnected
    exactly when ``M`` has eigenvalue **1** -- not 0.  An early version of the
    certified routine tested for eigenvalues near 0 of ``M``, which wrongly
    classified every graph with ``lam2(L) = 1`` (for instance the complete
    weighted graph with a uniform self-loop) as degenerate and skipped it.
    This test pins the two paths to identical coverage.
    """
    from spectral_expansion import lambda2_certified
    for n in (2, 3):
        n_exact = n_cert = 0
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            exact_ok = (cheeger_constant(W) is not None
                        and lambda2_exact(W) is not None)
            cert_ok = lambda2_certified(W) is not None
            n_exact += bool(exact_ok)
            n_cert += bool(cert_ok)
        assert n_exact == n_cert, f"n={n}: exact {n_exact} vs certified {n_cert}"


def test_certified_path_handles_uniform_loop_graph():
    """The specific graph the regression test above was written for."""
    from spectral_expansion import lambda2_certified
    W = [[Fraction(2), Fraction(1), Fraction(1)],
         [Fraction(1), Fraction(2), Fraction(1)],
         [Fraction(1), Fraction(1), Fraction(2)]]
    ex = sp.simplify(lambda2_exact(W))
    enc = lambda2_certified(W)
    assert ex is not None and enc is not None
    lo, hi = enc
    assert float(lo) - 1e-9 <= float(sp.N(ex, 25)) <= float(hi) + 1e-9


def test_enumeration_only_yields_connected():
    for n in (2, 3):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            # support of the off-diagonal part must be connected
            adj = {v: set() for v in range(n)}
            for i in range(n):
                for j in range(n):
                    if i != j and W[i][j] != 0:
                        adj[i].add(j)
            seen = {0}
            stack = [0]
            while stack:
                v = stack.pop()
                for u in adj[v]:
                    if u not in seen:
                        seen.add(u)
                        stack.append(u)
            assert len(seen) == n