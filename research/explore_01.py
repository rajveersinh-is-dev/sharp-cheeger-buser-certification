"""Exploration 01 -- what does the data say?

Enumerate every connected simple graph on n <= 6 vertices and measure, exactly:

  * lam2      : 2nd eigenvalue of the normalised Laplacian  I - D^-1/2 A D^-1/2
  * lam2_comb : 2nd eigenvalue of the combinatorial Laplacian  D - A
  * h_edge    : edge Cheeger constant, min_S |dS| / min(vol S, vol S^c)
  * h_vert    : vertex boundary Cheeger constant, min_S |N(S) \\ S| / min(|S|,|S^c|)

and report which configurations extremise  lam2 / (2h + 2h^2)  (Buser) and
lam2 / h (the bare ratio).  Purpose: discover the real extremal structure
before committing to a theorem, rather than guessing.

Pure enumeration; isoperimetric quantities are ratios of integers (exact in
floating point here); eigenvalues via numpy symmetric eigensolver.

Run:  python research/explore_01.py
"""

from __future__ import annotations

import itertools
import json
import sys

import numpy as np




def subsets_containing(n: int, pivot: int = 0):
    """All subsets of {0..n-1} containing `pivot`.

    The Cheeger quotient |dS| / min(vol S, vol S^c) is invariant under
    S <-> S^c, so these 2^(n-1) subsets suffice to realise every value.
    """
    others = [v for v in range(n) if v != pivot]
    out = []
    for r in range(len(others) + 1):
        for combo in itertools.combinations(others, r):
            out.append((pivot,) + combo)
    return out


def analyse_graph(nbr, subs) -> tuple:
    """Return (lam2, lam2_comb, h_edge, h_vert) or None if disconnected."""
    n = len(nbr)
    deg = [bin(x).count("1") for x in nbr]
    if min(deg) == 0:
        return None

    A = np.zeros((n, n))
    for i in range(n):
        r = nbr[i]
        while r:
            b = r & -r
            j = b.bit_length() - 1
            A[i, j] = 1.0
            r ^= b
    D = np.diag([float(x) for x in deg])
    L = D - A
    dinv = 1.0 / np.sqrt(np.array(deg, dtype=float))
    NL = (L * dinv[:, None]) * dinv[None, :]

    ev = np.linalg.eigvalsh(NL)
    if ev[0] > 1e-9:
        return None
    lam2 = float(ev[1])
    lam2_comb = float(np.linalg.eigvalsh(L)[1])

    volV = sum(deg)
    h_edge = float("inf")
    h_vert = float("inf")

    for S in subs:
        sb = 0
        for v in S:
            sb |= 1 << v
        volS = sum(deg[v] for v in S)
        volC = volV - volS

        e_in = 0
        for a, item in enumerate(S):
            for b in range(a + 1, len(S)):
                if (nbr[item] >> S[b]) & 1:
                    e_in += 1
        bnd = volS - 2 * e_in
        den = min(volS, volC)
        if den > 0:
            r = bnd / den
            if r < h_edge:
                h_edge = r

        vb = 0
        for v in range(n):
            if not ((sb >> v) & 1) and (nbr[v] & sb):
                vb += 1
        k = len(S)
        denv = min(k, n - k)
        if denv > 0:
            r = vb / denv
            if r < h_vert:
                h_vert = r

    return lam2, lam2_comb, h_edge, h_vert


def adj_string(nbr, n):
    """Adj string.
    
    Args:
        nbr:
        n:
    
    Returns:
        The computed result
    
    """
    rows = []
    for i in range(n):
        rows.append("".join("1" if (nbr[i] >> j) & 1 else "0" for j in range(n)))
    return "/".join(rows)


def main():
    """Entry point — parse arguments and run the main computation.
    
    """
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    recs = []
    for n in range(2, max_n + 1):
        subs = subsets_containing(n)
        m = n * (n - 1) // 2
        pairs = list(itertools.combinations(range(n), 2))
        count = 0
        for mask in range(1 << m):
            nbr = [0] * n
            for k, (i, j) in enumerate(pairs):
                if (mask >> k) & 1:
                    nbr[i] |= 1 << j
                    nbr[j] |= 1 << i
            res = analyse_graph(nbr, subs)
            if res is None:
                continue
            lam2, lam2_comb, h_edge, h_vert = res
            deg = sorted(bin(x).count("1") for x in nbr)
            recs.append(
                dict(n=n, mask=mask, lam2=lam2, lam2_comb=lam2_comb,
                     h_edge=h_edge, h_vert=h_vert, degseq=deg,
                     adj=adj_string(nbr, n))
            )
            count += 1
        print(f"n={n}: {count} connected labelled graphs", flush=True)

    print(f"\ntotal connected graphs enumerated: {len(recs)}\n")

    # attach the three ratios BEFORE any sorting
    for r in recs:
        h = r["h_edge"]
        r["r_buser"] = r["lam2"] / (2 * h + 2 * h * h) if h > 0 else float("nan")
        r["r_ratio"] = r["lam2"] / h if h > 0 else float("nan")
        r["r_vratio"] = r["lam2"] / r["h_vert"] if r["h_vert"] > 0 else float("nan")

    def show(title, key, k=10, rev=True, filt=None):
        """Show.
        
        Args:
            title:
            key:
            k (int):
            rev (bool):
            filt:
        
        """
        rs = recs if filt is None else [r for r in recs if filt(r)]
        rs = sorted(rs, key=lambda r: r[key], reverse=rev)[:k]
        print(f"--- top {k} by {key} ({title}) ---")
        for r in rs:
            print(f"  n={r['n']} deg={r['degseq']} lam2={r['lam2']:.6f} "
                  f"h_edge={r['h_edge']:.6f} h_vert={r['h_vert']:.4f} "
                  f"adj={r['adj']}")
        print()

    show("Buser ratio lam2/(2h+2h^2)", "r_buser", 12)
    show("bare ratio lam2/h_edge", "r_ratio", 12)
    show("lam2 with h_vert", "r_vratio", 12)

    # global maxima of the Buser ratio
    print("=== maxima over the whole enumeration ===")
    print(f"  max lam2/(2h+2h^2) = {max(r['r_buser'] for r in recs):.6f}")
    print(f"  max lam2/h_edge    = {max(r['r_ratio'] for r in recs):.6f}")
    print(f"  max lam2/h_vert    = {max(r['r_vratio'] for r in recs):.6f}")

    with open("results/explore_01_records.json", "w") as fh:
        json.dump(recs, fh)
    print("\nwrote results/explore_01_records.json")


if __name__ == "__main__":
    main()