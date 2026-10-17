"""Exhaustive verification of the two-sided Cheeger--Buser bound.

Runs over

  * every connected unweighted graph on n <= 6, and
  * every connected weighted graph with self-loops on n <= 4
    (off-diagonal weights in {1,2}, self-loop weights in {0,1,2}),

and checks, for each graph, both

    lam2 >= 1 - sqrt(h(2-h))      (Cheeger)      and
    lam2 <= 2h + 2h^2             (Buser)

with h computed in EXACT rational arithmetic.

Two modes:

  --mode exact      lambda_2 from an exact characteristic polynomial.  No floating
                    point enters any decision.  Correct but slow.
  --mode certified  lambda_2 from a certified high-precision enclosure.  Much
                    faster.  ``tests/test_spectral_expansion.py`` contains two
                    cross-validation tests proving the certified interval
                    encloses the exact eigenvalue, and that the two paths agree
                    on every graph with n <= 3.

Default is ``exact`` (the stronger guarantee).  Writes a machine-readable summary
to results/verification_summary.json.

Run:  python experiments/run_verification.py [--mode exact|certified] [--quick]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from spectral_expansion import (  # noqa: E402
    cheeger_constant, lambda2_exact, buser_bound, cheeger_bound,
    verify_bounds_certified, enumerate_weighted_graphs,
)

TOL = 1e-12


def check_exact(W):
    """Exact check. Returns dict, or None if the graph is degenerate."""
    h = cheeger_constant(W)
    if h is None:
        return None
    lam = lambda2_exact(W)
    if lam is None:
        return None
    H = sp.Rational(h.numerator, h.denominator)
    bs = sp.simplify(lam - buser_bound(H))
    cs = sp.simplify(lam - cheeger_bound(H))
    return dict(h=h, lam2=lam,
                buser_ok=sp.N(bs, 12) <= TOL,
                cheeger_ok=sp.N(cs, 12) >= -TOL,
                ratio=float(sp.N(lam / H, 20)))


def check_certified(W):
    """Certified check. Returns dict, or None if the graph is degenerate."""
    res = verify_bounds_certified(W)
    if res is None:
        return None
    h = res["h"]
    H = sp.Rational(h.numerator, h.denominator)
    mid = (res["lam2_lo"] + res["lam2_hi"]) / 2
    return dict(h=h, lam2=sp.Rational(mid.numerator, mid.denominator),
                buser_ok=res["buser_ok"], cheeger_ok=res["cheeger_ok"],
                ratio=float(sp.N(mid / H, 20)))


def sweep(kind: str, nmax: int, edge_weights, loop_weights, check):
    counts = {}
    bviol, cviol = [], []
    best = None
    for n in range(2, nmax + 1):
        cnt = bb = cc = 0
        for W in enumerate_weighted_graphs(n, edge_weights=edge_weights,
                                          loop_weights=loop_weights):
            r = check(W)
            if r is None:
                continue
            cnt += 1
            if not r["buser_ok"]:
                bb += 1
                if len(bviol) < 5:
                    bviol.append(dict(n=n, W=[[str(x) for x in row] for row in W],
                                     h=str(r["h"]), lam2=str(r["lam2"])))
            if not r["cheeger_ok"]:
                cc += 1
                if len(cviol) < 5:
                    cviol.append(dict(n=n, W=[[str(x) for x in row] for row in W],
                                      h=str(r["h"]), lam2=str(r["lam2"])))
            if best is None or r["ratio"] > best[0]:
                best = (r["ratio"], n, str(r["h"]), str(r["lam2"]))
        counts[n] = dict(graphs=cnt, buser_violations=bb, cheeger_violations=cc)
        print(f"    {kind} n={n}: {cnt} graphs, "
              f"Buser viol {bb}, Cheeger viol {cc}", flush=True)
    return counts, bviol, cviol, best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("exact", "certified"), default="exact")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default=os.path.join("results",
                                                  "verification_summary.json"))
    args = ap.parse_args()

    check = check_exact if args.mode == "exact" else check_certified
    n_unw, n_wt = (4, 3) if args.quick else (6, 4)

    t0 = time.time()
    print(f"Verification [{args.mode}]: unweighted n<={n_unw}, "
          f"weighted+loops n<={n_wt}")

    uw, ub, uc, ubest = sweep("unweighted", n_unw, (1,), (0,), check)
    wt, wb, wc, wbest = sweep("weighted+loops", n_wt, (1, 2), (0, 1, 2), check)

    allb = ub + wb
    allc = uc + wc
    best = max([x for x in (ubest, wbest) if x is not None],
               key=lambda t: t[0])

    total = sum(v["graphs"] for v in uw.values()) + \
        sum(v["graphs"] for v in wt.values())
    summary = dict(
        mode=args.mode,
        exact_arithmetic_for_h=True,
        eigenvalue_method=("exact characteristic polynomial" if args.mode == "exact"
                           else "certified high-precision enclosure"),
        unweighted=uw,
        weighted_with_loops=wt,
        violations=dict(buser=allb, cheeger=allc),
        violation_examples=dict(buser=allb[:5] if allb else [],
                                cheeger=allc[:5] if allc else []),
        max_ratio=dict(value=best[0], n=best[1], h=best[2], lam2=best[3]),
        total_exact_graphs=total,
        total_violations=len(allb) + len(allc),
    )
    summary["elapsed_seconds"] = round(time.time() - t0, 1)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\nmode                      : {args.mode}")
    print(f"total graphs verified     : {total}")
    print(f"total violations          : {summary['total_violations']}")
    print(f"max lam2/h                : {best[0]:.10f} "
          f"(n={best[1]}, h={best[2]}, lam2={best[3]})")
    print(f"elapsed                   : {summary['elapsed_seconds']}s")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()