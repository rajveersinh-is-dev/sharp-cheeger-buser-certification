'''Exhaustive verification of the two‑sided Cheeger–Buser bound.

The script enumerates all connected graphs up to a given size and checks the
two inequalities

    λ₂ ≥ 1 - sqrt(h·(2‑h))   (Cheeger)   and
    λ₂ ≤ 2·h + 2·h²          (Buser)

where ``h`` is the Cheeger constant computed in exact rational arithmetic.
Two modes are supported:

* ``exact`` – λ₂ is obtained from an exact characteristic polynomial (slow but
  fully rigorous).
* ``certified`` – λ₂ is obtained from a certified high‑precision enclosure
  (much faster).

The script writes a machine‑readable summary to ``results/verification_summary.json``.
'''  # noqa: D400, D401

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple, Optional

import sympy as sp

# Ensure the repository's ``src`` directory is on ``sys.path``.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from spectral_expansion import (  # noqa: E402
    cheeger_constant,
    lambda2_exact,
    buser_bound,
    cheeger_bound,
    verify_bounds_certified,
    enumerate_weighted_graphs,
)

# ---------------------------------------------------------------------------
# Constants & type aliases
# ---------------------------------------------------------------------------
TOL: float = 1e-12
Matrix = List[List[int]]
CheckFn = Callable[[Matrix], Optional[Dict[str, Any]]]

__all__ = [
    "check_exact",
    "check_certified",
    "sweep",
    "run_verification",
    "parse_args",
    "main",
]

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _safe_ratio(lam: sp.Rational, h: Fraction) -> float:
    """Return ``float(lam / h)`` handling the ``h == 0`` edge case.

    The Cheeger constant can be zero for a disconnected graph, but such graphs
    are filtered out earlier.  Guarding against division by zero makes the code
    robust if the underlying library ever returns ``0``.
    """
    if h == 0:
        return float("inf")
    return float(sp.N(lam / sp.Rational(h.numerator, h.denominator), 20))


# ---------------------------------------------------------------------------
# Core verification functions
# ---------------------------------------------------------------------------

def check_exact(W: Matrix) -> Optional[Dict[str, Any]]:
    """Perform an *exact* verification of the Cheeger–Buser bounds.

    Parameters
    ----------
    W:
        Weighted adjacency matrix of the graph (integer entries).

    Returns
    -------
    dict or None
        ``None`` is returned for degenerate graphs (e.g. disconnected) where the
        Cheeger constant or the second eigenvalue cannot be defined.
    """
    h: Optional[Fraction] = cheeger_constant(W)
    if h is None:
        return None
    lam: Optional[sp.Rational] = lambda2_exact(W)
    if lam is None:
        return None

    # Convert ``h`` to a SymPy Rational for exact arithmetic.
    H = sp.Rational(h.numerator, h.denominator)
    buser_diff = sp.simplify(lam - buser_bound(H))
    cheeger_diff = sp.simplify(lam - cheeger_bound(H))

    return {
        "h": h,
        "lam2": lam,
        "buser_ok": sp.N(buser_diff, 12) <= TOL,
        "cheeger_ok": sp.N(cheeger_diff, 12) >= -TOL,
        "ratio": _safe_ratio(lam, h),
    }


def check_certified(W: Matrix) -> Optional[Dict[str, Any]]:
    """Perform a *certified* verification using interval arithmetic.

    The function ``verify_bounds_certified`` returns a dictionary containing a
    lower and upper bound for ``λ₂`` together with boolean flags indicating whether
    the Buser and Cheeger inequalities hold.  This wrapper converts the interval
    midpoint to a rational number and assembles the same result structure as
    :func:`check_exact`.
    """
    res: Optional[Dict[str, Any]] = verify_bounds_certified(W)
    if res is None:
        return None
    h: Fraction = res["h"]
    H = sp.Rational(h.numerator, h.denominator)
    # Midpoint of the certified interval.
    mid = (res["lam2_lo"] + res["lam2_hi"]) / 2
    lam_mid = sp.Rational(mid.numerator, mid.denominator)
    return {
        "h": h,
        "lam2": lam_mid,
        "buser_ok": res["buser_ok"],
        "cheeger_ok": res["cheeger_ok"],
        "ratio": _safe_ratio(lam_mid, h),
    }


# ---------------------------------------------------------------------------
# Sweep utilities
# ---------------------------------------------------------------------------

def sweep(
    kind: str,
    nmax: int,
    edge_weights: Tuple[int, ...],
    loop_weights: Tuple[int, ...],
    check: CheckFn,
) -> Tuple[Dict[int, Dict[str, int]], List[Dict[str, Any]], List[Dict[str, Any]], Optional[Tuple[float, int, str, str]]]:
    """Iterate over all graphs of a given ``kind`` and collect statistics.

    Parameters
    ----------
    kind:
        Human‑readable identifier used only for progress output.
    nmax:
        Maximum number of vertices to enumerate.
    edge_weights:
        Allowed off‑diagonal weights.
    loop_weights:
        Allowed diagonal (self‑loop) weights.
    check:
        Callable performing the actual bound verification.

    Returns
    -------
    counts:
        Mapping ``n → {graphs, buser_violations, cheeger_violations}``.
    bviol, cviol:
        Lists of up to five violation examples for each inequality.
    best:
        Tuple ``(ratio, n, h, lam2)`` describing the graph with the largest
        ``λ₂ / h`` ratio, or ``None`` if no graph was processed.
    """
    counts: Dict[int, Dict[str, int]] = {}
    bviol: List[Dict[str, Any]] = []
    cviol: List[Dict[str, Any]] = []
    best: Optional[Tuple[float, int, str, str]] = None

    for n in range(2, nmax + 1):
        cnt = bcount = ccount = 0
        for W in enumerate_weighted_graphs(
            n, edge_weights=edge_weights, loop_weights=loop_weights
        ):
            result = check(W)
            if result is None:
                continue
            cnt += 1
            if not result["buser_ok"]:
                bcount += 1
                if len(bviol) < 5:
                    bviol.append(
                        {
                            "n": n,
                            "W": [[str(x) for x in row] for row in W],
                            "h": str(result["h"]),
                            "lam2": str(result["lam2"]),
                        }
                    )
            if not result["cheeger_ok"]:
                ccount += 1
                if len(cviol) < 5:
                    cviol.append(
                        {
                            "n": n,
                            "W": [[str(x) for x in row] for row in W],
                            "h": str(result["h"]),
                            "lam2": str(result["lam2"]),
                        }
                    )
            if best is None or result["ratio"] > best[0]:
                best = (
                    result["ratio"],
                    n,
                    str(result["h"]),
                    str(result["lam2"]),
                )
        counts[n] = {"graphs": cnt, "buser_violations": bcount, "cheeger_violations": ccount}
        print(f"    {kind} n={n}: {cnt} graphs, Buser viol {bcount}, Cheeger viol {ccount}", flush=True)
    return counts, bviol, cviol, best


# ---------------------------------------------------------------------------
# High‑level orchestration
# ---------------------------------------------------------------------------

def run_verification(mode: str, quick: bool, out_path: Path) -> Dict[str, Any]:
    """Execute the full verification suite and write the JSON summary.

    Parameters
    ----------
    mode:
        Either ``"exact"`` or ``"certified"``.
    quick:
        If ``True`` a reduced graph set is used (``n_unw=4``, ``n_wt=3``).
    out_path:
        Destination file for the JSON summary.
    """
    check: CheckFn = check_exact if mode == "exact" else check_certified
    n_unw, n_wt = (4, 3) if quick else (6, 4)

    start = time.time()
    print(
        f"Verification [{mode}]: unweighted n<={n_unw}, weighted+loops n<={n_wt}"
    )

    uw, ub, uc, ubest = sweep("unweighted", n_unw, (1,), (0,), check)
    wt, wb, wc, wbest = sweep("weighted+loops", n_wt, (1, 2), (0, 1, 2), check)

    all_b = ub + wb
    all_c = uc + wc
    best_candidates = [c for c in (ubest, wbest) if c is not None]
    best = max(best_candidates, key=lambda t: t[0]) if best_candidates else None

    total_graphs = sum(v["graphs"] for v in uw.values()) + sum(v["graphs"] for v in wt.values())
    summary: Dict[str, Any] = {
        "mode": mode,
        "exact_arithmetic_for_h": True,
        "eigenvalue_method": (
            "exact characteristic polynomial" if mode == "exact" else "certified high-precision enclosure"
        ),
        "unweighted": uw,
        "weighted_with_loops": wt,
        "violations": {"buser": all_b, "cheeger": all_c},
        "violation_examples": {
            "buser": all_b[:5] if all_b else [],
            "cheeger": all_c[:5] if all_c else [],
        },
        "max_ratio": {
            "value": best[0] if best else None,
            "n": best[1] if best else None,
            "h": best[2] if best else None,
            "lam2": best[3] if best else None,
        },
        "total_exact_graphs": total_graphs,
        "total_violations": len(all_b) + len(all_c),
        "elapsed_seconds": round(time.time() - start, 1),
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    # Human‑readable final report.
    print("\nmode                      :", mode)
    print("total graphs verified     :", total_graphs)
    print("total violations          :", summary["total_violations"])
    if best:
        print(
            f"max lam2/h                : {best[0]:.10f} (n={best[1]}, h={best[2]}, lam2={best[3]})"
        )
    print("elapsed                   :", summary["elapsed_seconds"], "s")
    print("wrote", out_path)
    return summary


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command‑line arguments.

    The function is isolated to make unit‑testing easier.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("exact", "certified"),
        default="exact",
        help="Verification mode (default: exact)",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run a reduced test suite for faster execution",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results", "verification_summary.json"),
        help="Output JSON file (default: results/verification_summary.json)",
    )
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    run_verification(args.mode, args.quick, args.out)


if __name__ == "__main__":
    main()
