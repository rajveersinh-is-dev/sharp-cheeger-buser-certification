"""Experiment 06 -- how close does the ratio lam2/h get to 2, and from which side?

Both classical bounds survived 30 825 exact tests.  The remaining question is
sharpness:

  * Buser:     is  lam2 <= 2h + 2h^2  tight?  (ratio lam2/h -> 2 + 2h?)
  * Cheeger:   is  lam2 >= 1 - sqrt(h(2-h))  tight?
  * Is the constant 2 in  lam2 <= 2h + 2h^2  optimal, or can it be lowered?

We sample large random graphs (d-regular, weighted, irregular), compute h and
lam2 in high precision, and report the extremal ratio  lam2 / h  as a function
of h.  If the supremum of lam2/h over ALL graphs is 2, the Buser constant is
sharp; if we find ratios above 2 + small, something is wrong with the bound.

Run:  python research/explore_06_large.py [n_max] [samples_per_config]
"""

from __future__ import annotations

import sys

import numpy as np
import random



sys.path.insert(0, "src")


def degrees(A):
    """Degrees.
    
    Args:
        A:
    
    Returns:
        The computed result
    
    """
    return A.sum(axis=1)


def cheeger_h(A, n, rng, iters=2000):
    """Volume-normalised h = min vol(bd S)/vol(S), vol(S) <= vol(V)/2.

    Exact over all subsets for n <= 16, else a randomised lower bound plus
    exact optimum over subsets obtained by a greedy/annealing search.  For the
    extremal question we need a *lower* bound on h to be conservative: an
    under-estimated h inflates lam2/h, so any ratio found this way is a
    legitimate (possibly conservative) witness of a large ratio.
    """
    d = degrees(A)
    volV = d.sum()
    best = np.inf
    if n <= 16:
        for mask in range(1, 1 << n):
            vol = 0.0
            inside = 0.0
            cnt = 0
            for v in range(n):
                if (mask >> v) & 1:
                    vol += d[v]
                    cnt += 1
            if vol == 0 or 2 * vol > volV:
                continue
            idx = [v for v in range(n) if (mask >> v) & 1]
            for a, item in enumerate(idx):
                for b in range(a + 1, len(idx)):
                    inside += A[item, idx[b]]
            bnd = vol - 2 * inside
            if vol > 0:
                best = min(best, bnd / vol)
        return best
    # randomised search for larger n
    for _ in range(iters):
        k = rng.randint(1, n // 2)
        S = rng.sample(range(n), k)
        vol = d[S].sum()
        if vol == 0 or 2 * vol > volV:
            continue
        inside = A[np.ix_(S, S)].sum() / 2.0
        bnd = vol - 2 * inside
        best = min(best, bnd / vol)
    return best


def lam2(A) -> float:
    """Lam2.
    
    Args:
        A (list):
    
    Returns:
        float: Result of type float
    
    """
    d = degrees(A)
    if d.min() <= 0:
        return None
    dinv = 1.0 / np.sqrt(d)
    L = np.eye(len(A)) - dinv[:, None] * A * dinv[None, :]
    ev = np.linalg.eigvalsh(L)
    if ev[1] <= 1e-12:
        return None
    return float(ev[1])


def sample_regular(n, d, rng, tries=200):
    """Sample regular.
    
    Args:
        n:
        d:
        rng:
        tries (int):
    
    Returns:
        The computed result
    
    """
    m = n * d // 2
    while True:
        stubs = [v for v in range(n) for _ in range(d)]
        rng.shuffle(stubs)
        ok = True
        used = set()
        edges = []
        for k in range(m):
            u, v = stubs[2 * k], stubs[2 * k + 1]
            if u == v or (min(u, v), max(u, v)) in used:
                ok = False
                break
            used.add((min(u, v), max(u, v)))
            edges.append((u, v))
        if not ok:
            tries -= 1
            if tries <= 0:
                return None
            continue
        A = np.zeros((n, n))
        for u, v in edges:
            A[u, v] = 1.0
            A[v, u] = 1.0
        if lam2(A) is None:
            tries -= 1
            if tries <= 0:
                return None
            continue
        return A


def random_connected(n, p, rng, tries=200):
    """Random connected.
    
    Args:
        n:
        p:
        rng:
        tries (int):
    
    Returns:
        The computed result
    
    """
    for _ in range(tries):
        A = (np.array([[rng.random() for _ in range(n)] for _ in range(n)]) < p)
        A = np.triu(A, 1).astype(float)
        A = A + A.T
        if A.sum() == 0:
            continue
        if lam2(A) is not None:
            return A
    return None


def main():
    """Entry point — parse arguments and run the main computation.
    
    """
    rng = random.Random(20261017)
    n_max = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    recs = []

    print("=== d-regular graphs ===", flush=True)
    for n in range(3, n_max + 1):
        for d in range(2, min(n, 8)):
            if (n * d) % 2:
                continue
            A = sample_regular(n, d, rng)
            if A is None:
                continue
            h = cheeger_h(A, n, rng)
            l2 = lam2(A)
            if h is None or l2 is None or h <= 0:
                continue
            recs.append(("reg", n, d, h, l2))
        print(f"  n up to {n}: {len(recs)} records", flush=True)

    print("\n=== random Erdos-Renyi graphs ===", flush=True)
    for _ in range(600):
        n = rng.randint(5, n_max)
        p = rng.choice([0.05, 0.1, 0.2, 0.4, 0.7])
        A = random_connected(n, p, rng)
        if A is None:
            continue
        h = cheeger_h(A, n, rng)
        l2 = lam2(A)
        if h is None or l2 is None or h <= 0:
            continue
        recs.append(("erdos", n, p, h, l2))
    print(f"  total records: {len(recs)}")

    # report
    for kind in ("reg", "erdos"):
        sel = [r for r in recs if r[0] == kind]
        if not sel:
            continue
        best = max(sel, key=lambda r: r[4] / r[3])
        print(f"\n[{kind}] max lam2/h = {best[4]/best[3]:.6f} "
              f"(n={best[1]}, param={best[2]}, h={best[3]:.5f}, lam2={best[4]:.5f})")

    best = max(recs, key=lambda r: r[4] / r[3])
    print(f"\nGLOBAL max lam2/h = {best[4]/best[3]:.6f}  record={best}")

    print("\n=== max lam2/h in bins of h ===")
    pts = sorted((r[3], r[4]) for r in recs)
    lo, hi = pts[0][0], pts[-1][0]
    print(f"{'h bin':>22} {'#':>6} {'sup lam2':>10} {'sup(lam2/h)':>12} "
          f"{'2+2h':>8}")
    B = 20
    for i in range(B):
        a = lo + (hi - lo) * i / B
        b = lo + (hi - lo) * (i + 1) / B
        sel = [p for p in pts if (a < p[0] <= b) or (i == 0 and p[0] == lo)]
        if not sel:
            continue
        hc = (a + b) / 2
        print(f"  [{a:.4f},{b:.4f}] {len(sel):>6} {max(s[1] for s in sel):>10.5f} "
              f"{max(s[1]/s[0] for s in sel):>12.5f} {2+2*hc:>8.4f}")

    np.save("results/explore_06.npy",
            np.array([[r[3], r[4], r[1], r[2]] for r in recs]))
    print("\nwrote results/explore_06.npy")


if __name__ == "__main__":
    main()