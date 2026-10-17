"""Re-verification after the boundary() fix, with an independent cross-check.

The first version of ``boundary`` incorrectly ADDED a half self-loop to the
boundary of a set.  That bug was invisible for loop-free graphs (the bulk of
the early experiments) but corrupted every loop-graph verification.  This
script

  1. cross-checks ``boundary`` against an independent edge-by-edge definition,
  2. re-runs the exhaustive exact verification of BOTH classical bounds,
  3. reports the corrected conclusions.
"""

from __future__ import annotations

import itertools
import sys
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "src")
from spectral_expansion import (  # noqa: E402
    boundary, degrees, volume, cheeger_constant, lambda2_exact,
    buser_bound, cheeger_bound, enumerate_weighted_graphs,
)


def independent_boundary(W, S):
    """bnd(S) = sum of weights of edges with exactly one endpoint in S.

    Self-loops have both endpoints in S (or neither), so they contribute 0.
    """
    S = set(S)
    tot = Fraction(0)
    n = len(W)
    for u in range(n):
        for v in range(u + 1, n):
            if (u in S) != (v in S):
                tot += W[u][v]
    return tot


def cross_check():
    print("=== 1. cross-check boundary() against an independent definition ===")
    bad = 0
    checked = 0
    for n in (2, 3, 4):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            for r in range(1, n // 2 + 1):
                for S in itertools.combinations(range(n), r):
                    checked += 1
                    if boundary(W, S) != independent_boundary(W, S):
                        bad += 1
    print(f"  checked {checked} (graph, set) pairs; mismatches = {bad}")
    assert bad == 0, "boundary() disagrees with the independent definition"
    print("  OK\n")


def identity_check():
    print("=== 2. identity  bnd(S) = vol(S) - 2E(S) - loops(S) ===")
    bad = 0
    for n in (2, 3, 4):
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            d = degrees(W)
            for r in range(1, n // 2 + 1):
                for S in itertools.combinations(range(n), r):
                    S = set(S)
                    Ein = sum((W[u][v] for u in S for v in S if u < v),
                              Fraction(0))
                    loops = sum((W[v][v] for v in S), Fraction(0))
                    rhs = sum((d[v] for v in S), Fraction(0)) - 2 * Ein - loops
                    if boundary(W, S) != rhs:
                        bad += 1
    print(f"  mismatches = {bad}")
    assert bad == 0
    print("  OK\n")


def verify_bounds(n_max_unw, n_max_wt):
    print("=== 3. exhaustive EXACT verification of both classical bounds ===")
    for tag, nmax, ew, lw in (
            ("unweighted", n_max_unw, (1,), (0,)),
            ("weighted+loops", n_max_wt, (1, 2), (0, 1, 2)),
    ):
        ntot = 0
        buser_bad = 0
        cheeger_bad = 0
        maxratio = None
        examples = []
        for n in range(2, nmax + 1):
            cnt = 0
            for W in enumerate_weighted_graphs(n, edge_weights=ew, loop_weights=lw):
                h = cheeger_constant(W)
                if h is None:
                    continue
                lam = lambda2_exact(W)
                if lam is None:
                    continue
                cnt += 1
                H = sp.Rational(h.numerator, h.denominator)
                bs = sp.simplify(lam - buser_bound(H))
                cs = sp.simplify(lam - cheeger_bound(H))
                if sp.N(bs, 12) > 1e-12:
                    buser_bad += 1
                    if len(examples) < 3:
                        examples.append(("BUSER", W, h, sp.simplify(lam), bs))
                if sp.N(cs, 12) < -1e-12:
                    cheeger_bad += 1
                    if len(examples) < 6:
                        examples.append(("CHEEGER", W, h, sp.simplify(lam), cs))
                r = sp.N(lam / H, 15)
                if maxratio is None or r > maxratio[0]:
                    maxratio = (r, n, h, sp.simplify(lam))
            ntot += cnt
            print(f"    {tag} n={n}: {cnt} graphs", flush=True)
        print(f"\n  [{tag}] total exact graphs = {ntot}")
        print(f"  [{tag}] Buser   lam2 <= 2h+2h^2        violations = {buser_bad}")
        print(f"  [{tag}] Cheeger lam2 >= 1-sqrt(h(2-h)) violations = {cheeger_bad}")
        if maxratio is not None:
            print(f"  [{tag}] max lam2/h = {maxratio[0]}  "
                  f"(n={maxratio[1]}, h={maxratio[2]}, lam2={maxratio[3]})")
        for e in examples:
            print(f"     {e[0]}: W={[[str(x) for x in r] for r in e[1]]} "
                  f"h={e[2]} lam2={e[3]} slack={e[4]}")
        print()


if __name__ == "__main__":
    cross_check()
    identity_check()
    n_unw = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    n_wt = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    verify_bounds(n_unw, n_wt)
    print("DONE")