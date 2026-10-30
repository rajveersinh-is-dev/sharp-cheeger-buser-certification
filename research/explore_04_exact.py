"""Exact (rational) exploration of Cheeger--Buser statements for graphs WITH self-loops.

Everything here is exact: weights are sympy Rationals, eigenvalues come from exact
characteristic polynomials, and Cheeger constants are ratios of integers.  The
verdict on each inequality is therefore not subject to floating-point error.

Conventions
-----------
A weighted graph is a symmetric matrix W >= 0; W[v][v] = self-loop weight at v.
    D_v        = sum_u W[v][u]        (a loop adds w_vv to the degree)
    vol(S)     = sum_{v in S} D_v
    boundary   vol(dS) = sum_{v in S, u not in S} W[v][u] + (1/2) sum_{v in S} W[v][v]
                         (a loop is cut in half)
    general normalized Laplacian  L = I - D^{-1/2} W D^{-1/2}

Two Cheeger constants are compared:
  h_plain(S) = vol(dS) / vol(S)
  h_loop(S)  = vol(dS) / ( vol(S) + 2 * vol(N(S)) )
        N(S) = S union { u not in S : exists v in S with W[u][v] > 0 }
(the second is the "loop-aware" constant of Lee--Oveis Gharan--Trevisan, which
divides by the volume of S together with twice the volume of its neighbourhood.)

Statements tested:
  (C) Cheeger : lam2 >= 1 - sqrt( 2 h (1 - h/2) )
  (B) Buser   : lam2 <=  2 h + 2 h^2

Run:  python research/explore_04_exact.py [n_max]
"""

from __future__ import annotations

import itertools
import sys
from fractions import Fraction

import sympy as sp


def volume(W, S):
    return sum(sum(W[v][u] for u in range(len(W))) for v in S)


def boundary(W, S):
    n = len(W)
    tot = Fraction(0)
    for v in S:
        for u in range(n):
            if u in S:
                if u == v:
                    tot += Fraction(W[v][v], 2)   # loop cut in half
            else:
                tot += Fraction(W[v][u])
    return tot


def neighbourhood_volume(W, S):
    n = len(W)
    N = set(S)
    for v in S:
        for u in range(n):
            if u not in S and W[u][v] != 0:
                N.add(u)
    return sum(sum(W[v][u] for u in range(n)) for v in N)


def cheeger_constants(W):
    """Return (h_plain, h_loop) minimised over 0 < vol(S) <= vol(V)/2."""
    n = len(W)
    tot = volume(W, range(n))
    best_plain = None
    best_loop = None
    for r in range(1, n // 2 + 1):
        for S in itertools.combinations(range(n), r):
            vS = volume(W, S)
            if vS == 0 or vS > tot / 2:
                continue
            bd = boundary(W, S)
            if bd == 0:
                continue
            hp = Fraction(bd, 1) / Fraction(vS, 1)
            den = vS + 2 * neighbourhood_volume(W, S)
            hl = Fraction(bd, 1) / Fraction(den, 1)
            if best_plain is None or hp < best_plain:
                best_plain = hp
            if best_loop is None or hl < best_loop:
                best_loop = hl
    return best_plain, best_loop


def gen_random_walk_laplacian(W):
    n = len(W)
    M = sp.zeros(n)
    for v in range(n):
        for u in range(n):
            M[v, u] = sp.Rational(W[v][u]).limit_denominator(10**6) \
                if isinstance(W[v][u], Fraction) else sp.Rational(W[v][u])
    d = [sum(M[v, u] for u in range(n)) for v in range(n)]
    if any(x <= 0 for x in d):
        return None          # a vertex with zero total degree: undefined
    for v in range(n):
        for u in range(n):
            M[v, u] = M[v, u] / d[v]
    return sp.eye(n) - M


def lam2_exact(W):
    L = gen_random_walk_laplacian(W)
    if L is None:
        return None
    x = sp.Symbol('x')
    vals = L.charpoly(x).all_roots()
    if any(v.is_real is False for v in vals):
        raise ValueError("non-real eigenvalue encountered")
    vals = sorted(vals, key=lambda x: float(sp.N(x, 25)))
    # discard the trivial eigenvalue 0
    pos = [v for v in vals if float(sp.N(v, 25)) > 1e-12]
    return pos[0] if pos else None


def check(W, label):
    n = len(W)
    hp, hl = cheeger_constants(W)
    if hp is None:
        return None
    lam2 = lam2_exact(W)
    if lam2 is None:
        return None
    out = dict(label=label, n=n, lam2=lam2, h_plain=hp, h_loop=hl)
    for tag, h in (("plain", hp), ("loop", hl)):
        hf = sp.Rational(h.numerator, h.denominator)
        cheeger_slack = sp.simplify(lam2 - (1 - sp.sqrt(2 * hf * (1 - hf / 2))))
        buser_slack = sp.simplify(lam2 - (2 * hf + 2 * hf ** 2))
        out[f"cheeger_slack_{tag}"] = cheeger_slack
        out[f"buser_slack_{tag}"] = buser_slack
    return out


def fmt(x):
    return sp.N(x, 12)


def main():
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    recs = []
    print("Enumerating weighted graphs WITH self-loops, exact rational arithmetic.")
    print(f"n_max = {n_max}\n")

    for n in range(2, n_max + 1):
        pairs = list(itertools.combinations(range(n), 2))
        # loop weights in {0, 1, 2} ; off-diagonal weights in {1, 2}
        combos = 1
        for _ in range(n):
            combos *= 3
        for edge_mask in range(1 << len(pairs)):
            for ew in itertools.product((1, 2), repeat=len(pairs)):
                pass  # too many; handle small n only
        # simpler: only n <= 3 exhaustive over weights, larger n -> simple families
        if n <= 3:
            for emask in range(1 << len(pairs)):
                for offw in itertools.product((1, 2), repeat=len(pairs)):
                    for loops in itertools.product((0, 1, 2), repeat=n):
                        W = [[0] * n for _ in range(n)]
                        for k, (i, j) in enumerate(pairs):
                            if (emask >> k) & 1:
                                W[i][j] = W[j][i] = offw[k]
                        for v in range(n):
                            W[v][v] = loops[v]
                        rec = check(W, f"n={n} e={emask} w={offw} l={loops}")
                        if rec:
                            recs.append(rec)
            print(f"  n={n}: {sum(1 for r in recs if r['n']==n)} exact records")
        else:
            break

    print(f"\ntotal exact records: {len(recs)}")

    print("\n=== sign summary of slacks (negative slack = VIOLATION) ===")
    for key in ("cheeger_slack_plain", "cheeger_slack_loop",
                "buser_slack_plain", "buser_slack_loop"):
        neg = [r for r in recs if float(sp.N(r[key], 15)) < -1e-12]
        zero = [r for r in recs if abs(float(sp.N(r[key], 15))) <= 1e-12]
        pos = [r for r in recs if float(sp.N(r[key], 15)) > 1e-12]
        print(f"  {key:>22}: violations {len(neg):>4}   equal {len(zero):>4}   strict {len(pos):>4}")

    print("\n=== examples of Buser violations (loop-aware h) ===")
    bad = [r for r in recs if float(sp.N(r["buser_slack_loop"], 15)) < -1e-12]
    for r in bad[:6]:
        print(f"  {r['label']}: lam2={fmt(r['lam2'])} h_loop={r['h_loop']} "
              f"slack={fmt(r['buser_slack_loop'])}")
    if not bad:
        print("  none")

    print("\n=== examples of Cheeger violations (loop-aware h) ===")
    badc = [r for r in recs if float(sp.N(r["cheeger_slack_loop"], 15)) < -1e-12]
    for r in badc[:6]:
        print(f"  {r['label']}: lam2={fmt(r['lam2'])} h_loop={r['h_loop']} "
              f"slack={fmt(r['cheeger_slack_loop'])}")
    if not badc:
        print("  none")


if __name__ == "__main__":
    main()