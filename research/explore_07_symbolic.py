"""Experiment 07 (fast) -- exact symbolic analysis of the two-vertex loop family.

Question
--------
With the *plain* volume isoperimetric constant

    h(S) = bnd(S)/vol(S),
    bnd(S) = sum_{v in S, u not in S} W[v][u] + (1/2) sum_{v in S} W[v][v],

is  lam2 <= 2h + 2h^2  (Buser)  valid for general weighted graphs with self-loops?
is  lam2 >= 1 - sqrt(h(2-h))  (Cheeger)  valid?

We analyse the two-parameter family

        W = [[a, w],
             [w, a]]        (a = loop weight, w = edge weight)

symbolically and exactly.  Everything below is sympy rational algebra.

Run:  python research/explore_07_symbolic.py
"""

from __future__ import annotations


import sympy as sp



def analyse(a_name="a", w_name="w", a_pos=True, w_pos=True) -> dict:
    """Analyse.
    
    Args:
        a_name (str):
        w_name (str):
        a_pos (bool):
        w_pos (bool):
    
    Returns:
        dict: Result of type dict
    
    """
    a = sp.Symbol(a_name, positive=True) if a_pos else sp.Symbol(a_name)
    w = sp.Symbol(w_name, positive=True) if w_pos else sp.Symbol(w_name)

    M = sp.Matrix([[a, w], [w, a]])
    D = sp.diag(a + w, a + w)
    Dih = D ** sp.Rational(-1, 2)
    L = sp.simplify(sp.eye(2) - Dih * M * Dih)

    evs = L.eigenvals()
    vals = []
    for val, mult in evs.items():
        for _ in range(mult):
            vals.append(sp.simplify(val))
    probe = {a: sp.Integer(1), w: sp.Integer(1)}
    vals.sort(key=lambda x: float(sp.N(x.subs(probe), 25)))
    lam1, lam2 = vals[0], vals[1]

    # degrees are (a+w); the only admissible non-trivial cut is S = {0},
    # vol(S) = a + w = vol(V)/2, bnd(S) = w  (the half-loop at 0 is internal).
    h = sp.simplify(w / (a + w))

    buser = 2 * h + 2 * h ** 2
    cheeger = 1 - sp.sqrt(h * (2 - h))

    print(f"--- family  W = [[a, w], [w, a]]  "
          f"(a = self-loop weight, w = edge weight) ---")
    print("eigenvalues of L :", lam1, ",", lam2)
    print("h (plain)        :", h)
    print()
    bs = sp.simplify(lam2 - buser)
    print("Buser   lam2 - (2h+2h^2) =", sp.factor(sp.together(bs)))
    num, den = sp.fraction(sp.together(bs))
    print("  numerator         =", sp.factor(num))
    print("  => Buser FAILS when", sp.factor(num), "> 0")

    cs = sp.simplify(lam2 - cheeger)
    print()
    print("Cheeger lam2 - (1 - sqrt(h(2-h))) =", sp.simplify(cs))
    print("  sign check: is it >= 0 for all a,w > 0?")
    print("  numerator rationalised =",
          sp.simplify(sp.radsimp(cs)))

    return dict(lam1=lam1, lam2=lam2, h=h, buser_slack=bs, cheeger_slack=cs)


def star_family():
    """Star / 'lollipop' family: irregular graphs without self-loops.

    The path P_n with a heavy loop at one end is the simplest source of
    irregularity.  Here we use the two-vertex graph with an extra pendant
    attached, but kept symbolic in the pendant weight.
    """
    print()
    print("--- three-vertex path with self-loops (t = loop weight at the ends) ---")
    t = sp.Symbol("t", positive=True)
    W = sp.Matrix([[t, 1, 0],
                   [1, 0, 1],
                   [0, 1, t]])
    D = sp.diag(t + 1, 2, t + 1)
    Dih = D ** sp.Rational(-1, 2)
    L = sp.simplify(sp.eye(3) - Dih * W * Dih)
    evs = L.eigenvals()
    vals = []
    for val, mult in evs.items():
        for _ in range(mult):
            vals.append(sp.simplify(val))
    vals.sort(key=lambda x: float(sp.N(x.subs({t: sp.Integer(1)}), 25)))
    print("eigenvalues :", [sp.N(v.subs({t: sp.Integer(1)}), 10) for v in vals])
    print("charpoly    :", sp.factor(sp.together(L.charpoly(t).as_expr())))

    # Cheeger constant: candidates are the singletons {0},{2} (vol t+1, bnd 1)
    # and {1} (vol 2, bnd 2), and pairs {0,1},{1,2} (vol t+3 > vol(V)/2).
    print("h({0}) = 1/(t+1) ,  h({2}) = 1/(t+1) ,  h({1}) = 2/2 = 1")
    print("h(S) = min( 1/(t+1), 1 ) = 1/(t+1)")


if __name__ == "__main__":
    analyse()
    star_family()