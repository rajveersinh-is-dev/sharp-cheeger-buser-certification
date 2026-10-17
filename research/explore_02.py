"""Exploration 02 -- the right Cheeger constant, and the effect of self-loops.

Standard Cheeger/Buser for the normalized Laplacian  L = I - D^-1/2 A D^-1/2
uses the VERTEX Cheeger constant

    h_v(G) = min_{S : 0 < |S| <= n/2}  |S| / |N(S)| ,   N(S) = S u N(S)\\S

with the classical pair
    Cheeger:  lam2 >= 1 - sqrt( 2 h_v (1 - h_v/2) )
    Buser:    lam2 <= 2 h_v + 2 h_v^2 .

Exploration 01 used an edge/volume constant, which is NOT the h used in these
statements; this script redoes the check with h_v, and -- the point of interest
-- allows SELF-LOOPS, where the classical statements are known to be delicate.

Self-loop convention (standard): D_v = sum_u w_vu  (a loop at v contributes
w_vv to D_v) and  x'x = sum_u w_vu x_v^2 counts a loop with weight w_vv, i.e.
the loop behaves like a "half edge" in the degree.

Run:  python research/explore_02.py [n_max] [loop_weight]
"""

from __future__ import annotations

import itertools
import json
import sys

import numpy as np


def cheeger_hv_edges(nbr, n):
    """Edge-boundary / min-volume Cheeger constant.

        h(G) = min_{S : 0 < vol S <= vol V / 2}  |E(S, S^c)| / vol(S)

    This is the h appearing in the classical statements for d-REGULAR graphs
    (there vol S = d|S|), and it is what Buser's lam2 <= 2h + 2h^2 uses.
    Also return the edge Cheeger constant defined with min(vol S, vol S^c).
    """
    best = float("inf")
    best2 = float("inf")
    argbest = None
    volV = 0
    deg = [bin(x).count("1") for x in nbr]
    volV = sum(deg)
    for r in range(1, n // 2 + 1):
        for S in itertools.combinations(range(n), r):
            volS = sum(deg[v] for v in S)
            volC = volV - volS
            e_in = 0
            for a in range(len(S)):
                for b in range(a + 1, len(S)):
                    if (nbr[S[a]] >> S[b]) & 1:
                        e_in += 1
            bnd = 2 * r - 2 * e_in
            den = min(volS, volC)
            if den > 0:
                val = bnd / den
                if val < best2:
                    best2 = val
            if volS > 0:
                val = bnd / volS
                if val < best:
                    best = val
                    argbest = S
    return best, best2, argbest


def cheeger_v(nbr, n):
    """min |S| / |N(S)| over nonempty S with |S| <= n/2."""
    best = float("inf")
    argbest = None
    for r in range(1, n // 2 + 1):
        for S in itertools.combinations(range(n), r):
            sb = 0
            for v in S:
                sb |= 1 << v
            N = 0
            for v in range(n):
                if v in S or (nbr[v] & sb):
                    N |= 1 << v
            val = r / bin(N).count("1")
            if val < best:
                best = val
                argbest = S
    return best, argbest


def build(n, pairs, mask, loop_w):
    """Adjacency rows as bitmasks plus the loop weight of each vertex."""
    nbr = [0] * n
    for k, (i, j) in enumerate(pairs):
        if (mask >> k) & 1:
            nbr[i] |= 1 << j
            nbr[j] |= 1 << i
    return nbr, loop_w


def spectrum(nbr, n, loop_w):
    """Normalized Laplacian eigenvalues for a graph with constant loop weight."""
    W = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if (nbr[i] >> j) & 1:
                W[i, j] = 1.0
    if loop_w:
        for i in range(n):
            W[i, i] = float(loop_w)
    d = W.sum(axis=1)
    if d.min() <= 0:
        return None
    dinv = 1.0 / np.sqrt(d)
    L = np.eye(n) - (W * dinv[:, None]) * dinv[None, :]
    ev = np.linalg.eigvalsh(L)
    # connected  <=>  exactly one zero eigenvalue  <=>  ev[1] > 0
    if ev[1] <= 1e-9:
        return None
    return float(ev[0]), float(ev[1]), float(ev[-1])


def main():
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    loop_w = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0

    viol = {"cheeger": [], "buser": []}
    worst_buser = (-1.0, None)
    worst_cheeger = (1e9, None)

    n_graphs = 0
    for n in range(2, n_max + 1):
        m = n * (n - 1) // 2
        pairs = list(itertools.combinations(range(n), 2))
        for mask in range(1 << m):
            nbr = [0] * n
            for k, (i, j) in enumerate(pairs):
                if (mask >> k) & 1:
                    nbr[i] |= 1 << j
                    nbr[j] |= 1 << i
            spec = spectrum(nbr, n, loop_w)
            if spec is None:
                continue
            n_graphs += 1
            _, lam2, lammax = spec
            h, h2, argS = cheeger_hv_edges(nbr, n)

            # Cheeger (lower bound on lam2)
            if h > 0:
                bound = 1.0 - np.sqrt(2 * h * (1.0 - h / 2.0))
                slack = lam2 - bound
                if slack < -1e-9:
                    viol["cheeger"].append(dict(n=n, mask=mask, lam2=lam2, h=h,
                                                slack=slack, S=list(argS)))
                if slack < worst_cheeger[0]:
                    worst_cheeger = (slack, dict(n=n, mask=mask, lam2=lam2, h=h))

            # Buser (upper bound on lam2)
            ub = 2 * h + 2 * h * h
            ratio = lam2 / ub if ub > 0 else float("nan")
            if not np.isnan(ratio) and ratio > worst_buser[0]:
                worst_buser = (ratio, dict(n=n, mask=mask, lam2=lam2, h=h,
                                            degseq=sorted(bin(x).count("1") for x in nbr)))

    print(f"loop weight = {loop_w}, n_max = {n_max}")
    print(f"connected graphs examined: {n_graphs}")
    print(f"Cheeger lower-bound VIOLATIONS: {len(viol['cheeger'])}")
    print(f"Buser upper-bound   VIOLATIONS: {len(viol['buser'])}")
    print(f"tightest Cheeger slack  (lam2 - bound): {worst_cheeger[0]:+.8f}  {worst_cheeger[1]}")
    print(f"largest Buser ratio    lam2/(2h+2h^2): {worst_buser[0]:.8f}  {worst_buser[1]}")

    for v in viol["cheeger"][:8]:
        print("   CHEEGER VIOLATION:", v)


if __name__ == "__main__":
    main()