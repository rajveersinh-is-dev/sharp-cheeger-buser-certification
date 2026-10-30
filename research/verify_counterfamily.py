"""Exact verification of the central phenomenon: Buser's inequality fails with self-loops.

Claim under test
----------------
For a finite connected graph with self-loops and its generalized normalized
Laplacian

    L = I - D^{-1/2} A D^{-1/2},     D_v = sum_u w_vu  (a loop at v adds w_vv),

let  h(G) = min_{S : 0 < vol(S) <= vol(V)/2}  vol(boundary S) / vol(S).

Classical Buser (d-regular, no loops):      lam2 <= 2h + 2 h^2.
Hypothesis: the same bound FAILS once self-loops are permitted.

Exact test family
-----------------
Two vertices, one edge of weight 1, a self-loop of weight w at each vertex.

    W = [[w, 1], [1, w]],   D = (1+w, 1+w),
    L = I - (1/(1+w)) * [[w,1],[1,w]].

Eigenvalues of (1+w)^{-1} [[w,1],[1,w]] are 1 (for (1,1)) and (w-1)/(1+w).
So lam2 = 1 - (w-1)/(1+w) = 2/(1+w) ... let us just compute symbolically below.

h(G): the only nontrivial cut is S={0}, vol(S) = D_0 = 1+w,
boundary weight = 1, so h = 1/(1+w).

We compute lam2 exactly with sympy as a function of w and compare with 2h+2h^2.

Everything here is exact (sympy rationals). No floating point is used for the
verdict.
"""

from __future__ import annotations

import sympy as sp

w = sp.symbols("w", positive=True)

W = sp.Matrix([[w, 1], [1, w]])
D = sp.diag(*[sum(W[i, j] for j in range(2)) for i in range(2)])
Dinv = D ** -1
L = sp.eye(2) - Dinv * W

print("W =", W)
print("degree vector =", [sp.simplify(D[i, i]) for i in range(2)])
print()
print("L =")
sp.pprint(L)
print()
print("charpoly of L :", sp.factor(L.charpoly(w).as_expr()))
ev = [sp.simplify(r) for r in L.eigenvals().keys()]
print("eigenvalues   :", ev)

lam2 = sp.nsimplify(sp.simplify([e for e in ev if e > 0][0]))
print("lam2 (exact)  :", lam2)

h = sp.Rational(1) / (1 + w)
print("h (exact)     :", h)

buser = sp.simplify(2 * h + 2 * h ** 2)
print("2h + 2h^2     :", buser)

slack = sp.simplify(lam2 - buser)
print("slack = lam2 - (2h+2h^2) :", sp.factor(slack))

print("\n--- where does the bound fail? ---")
# slack > 0  <=>  violation
num, den = sp.fraction(sp.together(slack))
print("slack > 0  <=>  ", sp.factor(num), "> 0")
roots = sp.solve(sp.Eq(num, 0), w)
print("slack = 0 at w =", roots)
crit = [r for r in roots if r.is_real]
if crit:
    print("positive critical value :", sp.simplify(crit[0]))
    print("numerically             :", sp.N(crit[0], 20))
    print("\nVerdict: bound holds for w <", sp.simplify(crit[0]),
          "and FAILS for w >", sp.simplify(crit[0]))
else:
    print("\nVerdict: bound never fails (slack is always <= 0 for w > 0)")