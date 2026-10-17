"""Exploration 03 -- where is Buser's constant 2 sharp?

Cheeger/Buser for the normalized Laplacian (Chung, Spectral Graph Theory, for
d-REGULAR graphs, h = min |E(S,S^c)| / vol(S)):

    Cheeger (lower):  lam2 >= 1 - sqrt( 2h (1 - h/2) )
    Buser   (upper):  lam2 <=  2h + 2h^2

This script:
  (a) enumerates d-regular graphs (n <= 7) and samples them for n = 8, 9, 10,
  (b) also samples random WEIGHTED d-regular graphs (all degree-1 weights 1/2),
  (c) computes the empirical extremal curve sup{lam2 : h} and compares with
      the classical 2h + 2h^2.

Purpose: find out empirically whether the constant 2 is attained, approached,
or exceeded.  Run:  python research/explore_03.py
"""

from __future__ import annotations

import itertools
import json
import random

import numpy as np

SEED = 20261017


def subsets_upto_half(n):
    out = []
    for r in range(1, n // 2 + 1):
        for S in itertools.combinations(range(1, n), r - 1):
            out.append((0,) + S)
    return out


def cheeger_regular(A, subs):
    d = A.sum(axis=1)
    best = float("inf")
    for S in subs:
        idx = list(S)
        inside = 0.0
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                inside += A[idx[a], idx[b]]
        cut = d[idx].sum() - 2 * inside
        volS = d[idx].sum()
        if volS > 0:
            val = cut / volS
            if val < best:
                best = val
    return best


def lam2_normalized(A):
    d = A.sum(axis=1)
    if d.min() <= 0:
        return None
    dinv = 1.0 / np.sqrt(d)
    L = np.eye(A.shape[0]) - (A * dinv[:, None]) * dinv[None, :]
    ev = np.linalg.eigvalsh(L)
    if ev[1] <= 1e-12:
        return None
    return float(ev[1])


def sample_regular(n, d, rng, tries=4000):
    """Yield distinct d-regular simple graphs via the configuration model."""
    m = n * d // 2
    seen = set()
    for _ in range(tries):
        stubs = [v for v in range(n) for _ in range(d)]
        rng.shuffle(stubs)
        edges, ok = [], True
        used = set()
        for k in range(m):
            u, v = stubs[2 * k], stubs[2 * k + 1]
            if u == v:
                ok = False
                break
            key = (u, v) if u < v else (v, u)
            if key in used:
                ok = False
                break
            used.add(key)
            edges.append(key)
        if not ok:
            continue
        fs = frozenset(edges)
        if fs in seen:
            continue
        seen.add(fs)
        A = np.zeros((n, n))
        for u, v in edges:
            A[u, v] = 1.0
            A[v, u] = 1.0
        yield A


def perturb_weights(A, rng):
    """Multiply each off-diagonal weight by 1 or 2 (stays d-regular)."""
    B = A.copy()
    n = A.shape[0]
    for u in range(n):
        for v in range(u + 1, n):
            if A[u, v] > 0:
                w = float(rng.choice([1.0, 2.0]))
                B[u, v] = B[v, u] = w
    return B


def main():
    rng = random.Random(SEED)
    records = []

    for n in range(3, 11):
        subs = subsets_upto_half(n)
        for d in range(1, n):
            if (n * d) % 2:
                continue
            if d == n - 1:
                continue  # complete graph: h=1, uninformative
            cnt = 0
            limit = 600 if n <= 7 else (400 if n <= 8 else 250)
            for A in sample_regular(n, d, rng, tries=60 * limit):
                lam2 = lam2_normalized(A)
                if lam2 is None:
                    continue
                h = cheeger_regular(A, subs)
                records.append(dict(kind="unweighted", n=n, d=d, h=h, lam2=lam2))
                # weighted variants
                for _ in range(2):
                    B = perturb_weights(A, rng)
                    l2 = lam2_normalized(B)
                    if l2 is not None:
                        records.append(dict(kind="weighted", n=n, d=d,
                                            h=cheeger_regular(B, subs), lam2=l2))
                cnt += 1
                if cnt >= limit:
                    break
        print(f"  n={n}: total records so far {len(records)}", flush=True)

    print(f"\ntotal records: {len(records)}")
    for r in records:
        h = r["h"]
        r["ratio"] = r["lam2"] / (2 * h + 2 * h * h) if h > 0 else float("nan")

    print("\n=== max ratio lam2/(2h+2h^2) by kind ===")
    for kind in sorted({r["kind"] for r in records}):
        sel = [r for r in records if r["kind"] == kind]
        w = max(sel, key=lambda r: r["ratio"])
        print(f"  {kind:>10}: max ratio {w['ratio']:.6f}  n={w['n']} d={w['d']} "
              f"h={w['h']:.6f} lam2={w['lam2']:.6f}")

    print("\n=== top 15 overall by ratio ===")
    for w in sorted(records, key=lambda r: -r["ratio"])[:15]:
        print(f"  {w['kind']:>10} n={w['n']:>2} d={w['d']:>2} h={w['h']:.6f} "
              f"lam2={w['lam2']:.6f} ratio={w['ratio']:.6f}")

    print("\n=== sup lam2 in bins of h vs Buser bound ===")
    pts = sorted(((r["h"], r["lam2"]) for r in records))
    lo, hi = pts[0][0], pts[-1][0]
    print(f"{'h bin':>22} {'#':>6} {'sup lam2':>10} {'2h+2h^2':>10} {'ratio':>8}")
    for i in range(20):
        a = lo + (hi - lo) * i / 20
        b = lo + (hi - lo) * (i + 1) / 20
        sel = [p[1] for p in pts if (a < p[0] <= b) or (i == 0 and p[0] == lo)]
        if not sel:
            continue
        hc = (a + b) / 2
        ub = 2 * hc + 2 * hc * hc
        print(f"  [{a:.4f}, {b:.4f}] {len(sel):>6} {max(sel):>10.5f} {ub:>10.5f} "
              f"{max(sel)/ub:>8.4f}")

    with open("results/explore_03_records.json", "w") as fh:
        json.dump(records, fh)
    print("\nwrote results/explore_03_records.json")


if __name__ == "__main__":
    main()