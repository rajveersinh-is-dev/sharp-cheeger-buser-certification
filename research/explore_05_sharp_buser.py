"""Decisive experiment: the exact sharp reverse (Buser) constant.

Question
--------
For finite weighted graphs (self-loops allowed) with the volume-normalised
Cheeger constant  h = min vol(boundary S)/vol(S) in (0,1], what is the true
supremum of  lam2 / h ?

The classical Buser bound gives  lam2 <= 2h + 2h^2, hence lam2/h <= 2 + 2h.

We compute, EXACTLY (rational arithmetic, exact characteristic polynomials),
the supremum of lam2/h over

  (a) every connected unweighted graph on n <= 6 vertices,
  (b) a large exact sample of weighted graphs with self-loops on n <= 5.

and we test candidate closed-form bounds.

Run:  python research/explore_05_sharp_buser.py [n_unweighted] [n_weighted]
"""

from __future__ import annotations

import sys

import sympy as sp

sys.path.insert(0, "src")
from spectral_expansion import (  # noqa: E402
    cheeger_constant, lambda2_exact, buser_bound, cheeger_bound,
    enumerate_weighted_graphs,
)


def main():
    """Entry point — parse arguments and run the main computation.
    
    Returns:
        The computed result
    
    """
    n_unw = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    n_wt = int(sys.argv[2]) if len(sys.argv) > 2 else 4

    records = []

    print(f"=== exact enumeration: unweighted, n <= {n_unw} ===", flush=True)
    for n in range(2, n_unw + 1):
        cnt = 0
        for W in enumerate_weighted_graphs(n):
            h = cheeger_constant(W)
            if h is None:
                continue
            lam = lambda2_exact(W)
            if lam is None:
                continue
            records.append((n, h, lam, "unweighted"))
            cnt += 1
        print(f"  n={n}: {cnt} connected graphs", flush=True)

    print(f"\n=== exact enumeration: weighted with self-loops, n <= {n_wt} ===",
          flush=True)
    for n in range(2, n_wt + 1):
        cnt = 0
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1)):
            h = cheeger_constant(W)
            if h is None:
                continue
            lam = lambda2_exact(W)
            if lam is None:
                continue
            records.append((n, h, lam, "weighted"))
            cnt += 1
        print(f"  n={n}: {cnt} graphs", flush=True)

    print(f"\ntotal exact records: {len(records)}")

    # ---- ratio lam2 / h -------------------------------------------------
    def ratio(h, lam):
        """Ratio.
        
        Args:
            h:
            lam:
        
        Returns:
            The computed result
        
        """
        return sp.N(lam / sp.Rational(h.numerator, h.denominator), 15)

    best = max(records, key=lambda r: float(ratio(r[1], r[2])))
    print(f"\nmax lam2/h = {float(ratio(best[1], best[2])):.6f}"
          f"   at n={best[0]} h={best[1]} lam2={sp.N(best[2],6)}"
          f"  ({best[3]})")

    # ---- test  Buser:  lam2 <= 2h + 2h^2 -------------------------------
    viol = []
    for n, h, lam, kind in records:
        H = sp.Rational(h.numerator, h.denominator)
        if sp.simplify(lam - buser_bound(H)) > 0:
            viol.append((n, h, lam, kind, sp.simplify(lam - buser_bound(H))))
    print(f"\nBuser  lam2 <= 2h+2h^2 : violations = {len(viol)} / {len(records)}")
    for v in viol[:
        8]:
        print(f"   n={v[0]} h={v[1]} lam2={sp.N(v[2],8)} slack={sp.N(v[4],8)}")

    # ---- test  Cheeger: lam2 >= 1 - sqrt(h(2-h)) -----------------------
    violC = []
    for n, h, lam, kind in records:
        H = sp.Rational(h.numerator, h.denominator)
        if sp.simplify(lam - cheeger_bound(H)) < 0:
            violC.append((n, h, lam, kind))
    print(f"Cheeger lam2 >= 1-sqrt(h(2-h)) : violations = "
          f"{len(violC)} / {len(records)}")
    for v in violC[:
        8]:
        print(f"   n={v[0]} h={v[1]} lam2={sp.N(v[2],8)}")

    # ---- empirical extremal curve ---------------------------------------
    print("\n=== sup lam2 as a function of h (binned) ===")
    pts = sorted(((float(h), float(sp.N(lam, 15))) for n, h, lam, k in records))
    lo, hi = pts[0][0], pts[-1][0]
    print(f"{'h bin':>20} {'#':>6} {'sup lam2':>10} {'2h+2h^2':>10} {'ratio':>8}")
    B = 16
    for i in range(B):
        a = lo + (hi - lo) * i / B
        b = lo + (hi - lo) * (i + 1) / B
        sel = [y for x, y in pts if (a < x <= b) or (i == 0 and x == lo)]
        if not sel:
            continue
        hc = (a + b) / 2
        print(f"  [{a:.4f},{b:.4f}] {len(sel):>6} {max(sel):>10.5f} "
              f"{2*hc+2*hc*hc:>10.5f} {max(sel)/(2*hc+2*hc*hc):>8.4f}")


if __name__ == "__main__":
    main()