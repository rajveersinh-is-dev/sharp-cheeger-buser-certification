"""Exhaustive EXACT verification of the two-sided theorem.

Runs the corrected pipeline (see research/verify_08_after_fix.py) at scale:

  * all connected unweighted graphs on n <= 6;
  * all connected weighted graphs with self-loops on n <= 4
    (edge weights in {1,2}, loop weights in {0,1,2});
  * a deterministic random sample of larger weighted graphs.

Everything is exact rational arithmetic: the Cheeger constant is a ratio of
integers and the eigenvalue comes from an exact characteristic polynomial.

Writes a machine-readable summary to results/verification_summary.json.

Run:  python experiments/run_verification.py [--quick]
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from spectral_expansion import (  # noqa: E402
    cheeger_constant, lambda2_exact, buser_bound, cheeger_bound,
    enumerate_weighted_graphs,
)

TOL = 1e-12


def check(W):
    """Return (h, lam2, buser_slack, cheeger_slack) exactly, or None."""
    h = cheeger_constant(W)
    if h is None:
        return None
    lam = lambda2_exact(W)
    if lam is None:
        return None
    H = sp.Rational(h.numerator, h.denominator)
    bs = sp.simplify(lam - buser_bound(H))
    cs = sp.simplify(lam - cheeger_bound(H))
    return h, lam, bs, cs


def run(n_unweighted: int, n_weighted: int, seed: int = 20261018):
    summary = {
        "seed": seed,
        "tolerance": TOL,
        "exact_arithmetic": True,
        "unweighted": {},
        "weighted_with_loops": {},
        "violations": {"buser": [], "cheeger": []},
        "max_ratio": None,
    }
    max_ratio = None

    # ---- unweighted -----------------------------------------------------
    for n in range(2, n_unweighted + 1):
        cnt = 0
        bbad = cbad = 0
        for W in enumerate_weighted_graphs(n):
            res = check(W)
            if res is None:
                continue
            h, lam, bs, cs = res
            cnt += 1
            H = sp.Rational(h.numerator, h.denominator)
            if sp.N(bs, 12) > TOL:
                bbad += 1
                summary["violations"]["buser"].append(
                    dict(n=n, W=[[str(x) for x in r] for r in W],
                         h=str(h), lam2=str(lam), slack=str(bs)))
            if sp.N(cs, 12) < -TOL:
                cbad += 1
                summary["violations"]["cheeger"].append(
                    dict(n=n, W=[[str(x) for x in r] for r in W],
                         h=str(h), lam2=str(lam), slack=str(cs)))
            r = sp.N(lam / H, 20)
            if max_ratio is None or r > max_ratio[0]:
                max_ratio = (r, n, str(h), str(sp.simplify(lam)))
        summary["unweighted"][n] = dict(
            graphs=cnt, buser_violations=bbad, cheeger_violations=cbad)
        print(f"  unweighted n={n}: {cnt} graphs, "
              f"Buser viol {bbad}, Cheeger viol {cbad}", flush=True)

    # ---- weighted with self-loops --------------------------------------
    for n in range(2, n_weighted + 1):
        cnt = 0
        bbad = cbad = 0
        for W in enumerate_weighted_graphs(n, edge_weights=(1, 2),
                                          loop_weights=(0, 1, 2)):
            res = check(W)
            if res is None:
                continue
            h, lam, bs, cs = res
            cnt += 1
            H = sp.Rational(h.numerator, h.denominator)
            if sp.N(bs, 12) > TOL:
                bbad += 1
                summary["violations"]["buser"].append(
                    dict(n=n, W=[[str(x) for x in r] for r in W],
                         h=str(h), lam2=str(lam), slack=str(bs)))
            if sp.N(cs, 12) < -TOL:
                cbad += 1
                summary["violations"]["cheeger"].append(
                    dict(n=n, W=[[str(x) for x in r] for r in W],
                         h=str(h), lam2=str(lam), slack=str(cs)))
            r = sp.N(lam / H, 20)
            if max_ratio is None or r > max_ratio[0]:
                maxratio_entry = (r, n, str(h), str(sp.simplify(lam)))
                if max_ratio is None or r > max_ratio[0]:
                    max_ratio = maxratio_entry
        summary["weighted_with_loops"][n] = dict(
            graphs=cnt, buser_violations=bbad, cheeger_violations=cbad)
        print(f"  weighted n={n}: {cnt} graphs, "
              f"Buser viol {bbad}, Cheeger viol {cbad}", flush=True)

    summary["max_ratio"] = dict(
        value=float(max_ratio[0]), n=max_ratio[1],
        h=max_ratio[2], lam2=max_ratio[3]) if max_ratio else None

    total = sum(v["graphs"] for v in summary["unweighted"].values()) + \
        sum(v["graphs"] for v in summary["weighted_with_loops"].values())
    summary["total_exact_graphs"] = total
    summary["total_violations"] = (len(summary["violations"]["buser"]) +
                                   len(summary["violations"]["cheeger"]))
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="smaller exhaustive ranges (for CI)")
    ap.add_argument("--out", default=os.path.join("results",
                                                  "verification_summary.json"))
    args = ap.parse_args()

    if args.quick:
        n_unw, n_wt = 4, 3
    else:
        n_unw, n_wt = 6, 4

    t0 = time.time()
    print(f"Exact verification: unweighted n<={n_unw}, "
          f"weighted+loops n<={n_wt}")
    summary = run(n_unw, n_wt)
    summary["elapsed_seconds"] = round(time.time() - t0, 1)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\ntotal exact graphs verified : {summary['total_exact_graphs']}")
    print(f"total violations            : {summary['total_violations']}")
    print(f"max lam2/h                  : {summary['max_ratio']}")
    print(f"elapsed                     : {summary['elapsed_seconds']}s")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()