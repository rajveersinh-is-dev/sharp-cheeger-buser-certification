# Candidate generation, screening and selection

Instruction: generate >= 10 serious candidates, screen them, rank with the 100-point
rubric, and explain why the winner was chosen. Search capability limits are recorded
in `00_search_capability_notes.md`.

Scoring rubric (weights): Credible novelty 25, Significance 20, Feasibility in one
session 15, Rigorous verifiability 15, Scholarly audience 10, Reproducibility 10,
Journal fit 5.  Scores are decision aids, not measurements of scientific quality.

## Stage A -- candidate generation

**C1. Sharp Cheeger--Buser constants for the normalized Laplacian (general weighted graphs).**
Q: what is the exact optimal constant in `lam2 <= C h + C h^2` for d-regular graphs?
State of the art: Cheeger 1970, Buser 1982, Alon--Milman 1985, Dodziuk 1984, Chung's
book; the constant 2 is *known* to be asymptotically tight. Gap: none -- it is solved.
Failure mode: cannot claim novelty.

**C2. Higher-order Cheeger inequalities with optimal constants (hypergraphs).**
Q: sharp constant in the k-th order Cheeger bound. Literature is very active
(Fan--Cheeger--Li; Lee--Oveis Gharan--Trevisan; Lou--Tang--Yip; Louis--Makarychev).
Exact optimality is open. Failure mode: too hard to close in one session; high risk.

**C3. Cheeger--Buser inequalities for graphs with SELF-LOOPS (vertex-weighted).**
Q: do the classical statements survive loops, and what are the correct sharp forms?
State of the art: Lee, Oveis Gharan & Trevisan (JCTB 2014) define a loop-aware
Cheeger constant; Salez (Combinatorica 2015) proves the sharp sqrt(2) rounding
constant. So the *Cheeger* direction is largely settled. The **reverse (Buser)
direction** is much less explored. -> selected; see below.

**C4. Second gap from 1 (dual Cheeger) of the normalized Laplacian.**
Q: relate `2 - lam_n` to bipartiteness ratios. Literature: Bauer--Jost, Trevisan,
Pokharanakar (Math. Z. 2026, graphons). Extremely fresh; would risk colliding and the
proof machinery is heavy. Failure mode: infeasible + priority risk.

**C5. Cheeger inequalities with vertex weights / non-uniform measures.**
Q: generalize to vertex-weighted expanders. Active (Lau--Kwok--Tung; Olesker-Tzang).
Gap hard to isolate.

**C6. Sharp order barriers for partitioned Runge--Kutta methods.**
Q: maximal achievable convergence order for a given method class. Active, but needs
specialised background and is far from OUP math journals in style.

**C7. Spectral bounds on multiway Cheeger constants via cyclomatic number.**
Already done in arXiv:2409.07097. Not novel.

**C8. Houdre--Tetali isoperimetric conjecture.**
Settled/counterexemplified by Lau--Tjowasi (arXiv:2407.11357). Not novel.

**C9. Certified exact computation of Cheeger constants / extremal spectral ratios.**
Q: can one *certify* sharp constants numerically with error bounds? Methodologically
useful and very verifiable, but on its own it is software, not a theorem.

**C10. Hitting-time / mixing-time characterisation via vertex conductance.**
Olesker-Tzang & Zanetti (PTRF 2024). Not novel.

**C11. Cheeger for simplicial complexes / Hodge Laplacians.**
Extremely active (Jost--Zhang, Vigano--Birdal--Schaub--Barahona). Overcrowded.

**C12. Sharp constants for the *sweep*/one-sided Cheeger inequality.**
Q: relate the second eigenvalue to the one-sided sweep constant. Ties to C3 but the
standard two-sided constants are already sharp; marginal.

## Stage B -- what the evidence actually showed

Rather than trusting recall, I ran two exhaustive/sampled computational studies
(`research/explore_01..03.py`). Findings that materially changed the choice:

1. **Cheeger's lower bound `lam2 >= 1 - sqrt(2h(1-h/2))` is FALSE for general
   irregular graphs** (1440 violations among all 27475 connected graphs on n <= 6,
   e.g. n=6, lam2=0.2362, h=0.2857, violation 0.0639). This is a *known*
   irregularity artifact -- the classical statements are made for d-regular graphs --
   and confirms that any honest statement must fix the hypotheses precisely.

2. **Buser's upper bound `lam2 <= 2h + 2h^2` is FALSE once self-loops are allowed.**
   Verified numerically on the 2-vertex family with loop weight w:
   `lam2 = w/(1+w)` while `h = 1/(1+w)`, so the bound requires
   `w^2 - w - 4 <= 0`, i.e. it fails for `w > (1+sqrt(17))/2 ~ 2.56`.
   Confirmed numerically (loop weight 3, n <= 6: max ratio 1.143).

3. **Buser's upper bound is also violated by weighted *regular* graphs**
   (sampled max ratio 1.0191), i.e. irregularity is not the whole story.

4. For genuinely d-regular *unweighted* graphs (n <= 10, ~9k graphs) the classical
   Buser constant 2 is never approached from above (max ratio 0.80): consistent with
   the literature, and a sanity check that my code is right.

The interesting, under-explored object is therefore **the reverse (Buser) direction for
weighted graphs with self-loops**: what is the *correct* upper bound on the spectral
gap in terms of a loop-aware isoperimetric constant, and what is the optimal constant?

## Stage C -- rubric ranking

| id | novelty 25 | signif 20 | feas 15 | verif 15 | audience 10 | repro 10 | fit 5 | total |
|----|-----------|----------|---------|----------|------------|----------|------|-------|
| C1  | 3  | 12 | 14 | 14 | 9 | 9 | 5 | 66 |
| C2  | 20 | 18 | 4  | 8  | 9 | 8 | 4 | 71 |
| C3  | 19 | 15 | 13 | 15 | 9 | 10 | 5 | 86 |
| C4  | 18 | 15 | 4  | 7  | 8 | 8 | 5 | 65 |
| C5  | 11 | 12 | 10 | 11 | 8 | 8 | 4 | 64 |
| C6  | 16 | 14 | 8  | 10 | 8 | 8 | 3 | 67 |
| C7  | 1  | 8  | 12 | 12 | 8 | 8 | 3 | 52 |
| C8  | 2  | 9  | 8  | 10 | 8 | 8 | 4 | 49 |
| C9  | 14 | 9  | 14 | 15 | 7 | 10 | 3 | 72 |
| C10 | 2  | 10 | 10 | 10 | 8 | 8 | 3 | 51 |
| C11 | 5  | 14 | 7  | 9  | 9 | 8 | 4 | 56 |
| C12 | 12 | 11 | 12 | 13 | 8 | 9 | 4 | 69 |

## Selection

**Chosen: C3** -- the reverse/Buser direction for weighted graphs with self-loops,
with a loop-aware isoperimetric constant and an exactly determined optimal constant.

Reasons, in order of weight:

* **Verifiability (15/15).** Every claim is a statement about finite weighted graphs.
  It is checkable by exhaustive enumeration with *exact rational arithmetic* plus
  certified numerics. There is no modelling assumption and no data.
* **Feasibility (13/15).** The core mathematics is the co-area (two-point) reduction,
  which is short and can be written out in full with no external black boxes. The
  sharp constant reduces to a one-dimensional optimisation that can be certified to
  arbitrary precision.
* **Novelty (19/25).** Not "nobody has looked": the loop-aware Cheeger constant and the
  sharp `sqrt(2)` rounding constant in the *Cheeger* direction are due to Lee--Oveis
  Gharan--Trevisan and Salez. My claim must be narrower and will be stated narrowly:
  a correct *upper* (reverse) bound in this setting, and the exactly determined optimal
  constant. The searches in Section 5 of `literature_review.md` did not surface a
  reverse-direction treatment of loop graphs, but absence of a hit is not proof.
* **Significance (15/20).** Self-loops / vertex weights are not cosmetic: they model
  laziness, self-transition in Markov chains, degree heterogeneity after normalisation,
  and the multicut / MaxCut setting. A bound that is *wrong* for such graphs (as the
  classical one is) has practical bite.
* **Audience/fit (9/10, 5/5).** Spectral graph theory / discrete geometry has a clear
  readership and clear OUP homes (JLMS, IMRN, Mathematika, QJM).

Runner-up C2 (higher-order sharp constants) scored higher on novelty but far lower on
feasibility and verifiability: closing it properly needs machinery I could not honestly
complete and verify in one session. C9 (certified computation) is retained as a
*methodological component* of C3 rather than as the headline.

**Honesty note.** C1, C7, C8, C10 were rejected as *not novel* because prior work
already contains the result, not because they were hard.