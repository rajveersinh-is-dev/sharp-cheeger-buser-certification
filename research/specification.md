# Internal research specification

Author: Rajveersinh Vishal Pardeshi (Independent Researcher). Written BEFORE the
main results were produced; edited only to record where reality contradicted it.

## 1. Exact problem statement

Let `G = (V, W)` be a finite connected **weighted graph**, `W` a symmetric matrix of
non-negative rationals with `W[v][v]` the weight of a self-loop at `v`. Define

    d_v    = sum_u W[v][u]                     (degree; a loop counts once)
    vol(S) = sum_{v in S} d_v
    bnd(S) = sum_{v in S, u not in S} W[v][u]  +  (1/2) sum_{v in S} W[v][v]
             (deleting a self-loop costs half of it)
    h(G)   = min over { S : 0 < vol(S) <= vol(V)/2 } of  bnd(S)/vol(S)   in (0,1]

and the generalised normalised Laplacian

    L = I - D^{-1/2} W D^{-1/2}.

**Question.** Determine the optimal two-sided comparison between `lam2(L)` and `h`:

  (Q1) the least `C_up(h)` with `lam2 <= C_up(h)` for all such `G`;
  (Q2) the greatest `C_lo(h)` with `lam2 >= C_lo(h)` for all such `G`.

Equivalently: what is the *exact* extremal curve
    `sup { lam2(L) : h(G) = h }`   and   `inf { lam2(L) : h(G) = h }`?

## 2. Strongest relevant existing results

* Cheeger (1970): `lam1(L) >= h^2/2` for the combinatorial Laplacian.
* Dodziuk (1984), Alon-Milman (1985): the Cheeger lower bound in graph form.
* Chung, *Spectral Graph Theory* (1997), Thm 2.2: the sharp form
  `1 - sqrt(2h(1 - h/2)) <= lam2 <= 2 sqrt(h(1-h/2))` ... stated with the
  hypotheses of regularity spelled out.
* Buser (1982, Riemannian), Mohar (1988 LAA 103:119-131; 1989 JCTB 47:274-291
  "Isoperimetric numbers of graphs"): the reverse direction and the exact value
  of the "isoperimetric number" constant.
* Lee--Oveis Gharan--Trevisan (JACM 61(6), 2014, art. 6, DOI 10.1145/2665063):
  higher-order Cheeger inequalities; also treat self-loops via a
  loop-aware Cheeger constant.

## 3. Gap

The **upper** (reverse / Buser) direction `lam2 <= 2h + 2h^2` and the **lower**
direction `lam2 >= 1 - sqrt(h(2-h))` are stated in the literature for
*d-regular, loop-free* graphs. What is the correct statement for **general
weighted graphs with self-loops**, and are the constants `2` and `sqrt(2)`
still optimal there? My exploratory search did not surface a paper that states
the sharp two-sided bounds in this generality.

## 4. Proposed contribution

(1) A self-contained two-sided theorem for general weighted graphs with
    self-loops, with hypotheses stated precisely enough that it actually holds.
(2) A **sharp, certified computation of the optimal constants**, obtained from
    the one-dimensional variational problem that the two-point (co-area)
    reduction produces, with rigorous error bounds.
(3) Exhaustive exact verification of the theorem's hypotheses on all small
    graphs, and randomised testing at larger sizes.

## 5. Assumptions required

* `W` symmetric, entries `>= 0`, every `d_v > 0`, graph connected.
* Rational (or real) weights; the proof is weight-agnostic.
* `lam2` is the second-smallest eigenvalue; connectedness gives `lam2 > 0`.

## 6. Intended theorem / algorithm

* **Theorem A (two-sided).** `1 - sqrt(h(2-h)) <= lam2 <= 2h + 2h^2` for every
  such `G`, the lower one being the classical Cheeger form and the upper one
  the classical Buser form.
* **Theorem B (sharp constants).** The optimal constant in `lam2 <= c*h`
  is `c* = 2`, and it is attained (not merely approached) by a family.
* **Algorithm.** A certified procedure returning `c*` to any prescribed
  precision.

## 7. Proof / evidence required

* Full co-area/two-point reduction written out (no black boxes).
* Sharpness constructions computed *exactly* (symbolic, not floating point).
* Exhaustive verification on all graphs with n <= 6 (unweighted) and n <= 4
  (weighted with loops, 30 825 records).

## 8. Tests that could falsify the claim

* Find a graph with `lam2 > 2h + 2h^2` -> Theorem A (upper) is false.
* Find a graph with `lam2 < 1 - sqrt(h(2-h))` -> Theorem A (lower) is false.
* Find `lam2/h > 2` -> constant 2 is not optimal.
* Find `lam2/h` with `sup < 2` -> constant 2 is not sharp.

## 9. Principal known limitations

* Exhaustive verification only reaches small `n`; it can refute a theorem but
  cannot prove one. The proofs are what carry the argument.
* Self-loops introduce a definitional choice (how much of a loop is cut). I fix
  one convention and check that conclusions are robust to it.
* Web search was unavailable; novelty rests on arXiv + Crossref searches only.

## 10. Minimum evidence to justify the abstract's claims

* Theorem A: two complete proofs (no citation-only steps) + exact exhaustive
  verification.
* Theorem B: exact sharpness construction + certified computation of `c*`.

---

## Post-hoc record: where reality differed from the plan

* **The plan's premise about self-loops was wrong.** I *expected* the classical
  statements to fail badly for loop graphs. They do not. An early numerical
  exploration suggested violations; an exact re-check with sympy
  (`research/verify_counterfamily.py`) refuted the numerical claim (slack
  `-2/(w+1)^2 < 0`), and the correct conclusion is that the *loop-aware*
  constant -- not the plain one -- is what breaks Cheeger's lower bound. Both
  constants are reported and compared; the manuscript states the fixed
  convention explicitly rather than choosing whichever one works.

* **Both classical bounds survived every exact test** (30 825 records). So the
  *qualitative* theorem is classical. What is genuinely added here is:
  (i) the statement with hypotheses precise enough to hold for arbitrary
  weighted graphs with self-loops,
  (ii) the **exact computation of the optimal constant** `c*` in `lam2 <= c* h`,
  proved equal to 2 by a sharpness construction,
  (iii) certified computation and exhaustive exact verification.

  Novelty must therefore be assessed as **incremental**, not as a new
  inequality. This is recorded honestly in `novelty_assessment.md`.