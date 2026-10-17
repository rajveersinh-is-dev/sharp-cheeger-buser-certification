# Hostile internal peer review

Two independent adversarial passes, written in the role of a sceptical referee
with expertise in spectral graph theory. These are **internal adversarial
checks**. No external or actual peer review has taken place, and nothing here
should be read as referee feedback from anyone else.

Each objection is classified **Critical / Major / Moderate / Minor**, with the
evidence that would resolve it and the action taken.

---

# Review 1 — Technical validity

## R1.1 (Critical) — "A new inequality" claim

*Objection.* The abstract of an earlier draft of this work implied the two-sided
bound might be new. It is not: it is Buser's and Cheeger's classical pair,
treated in Chung's *Spectral Graph Theory*. Publishing it as a theorem would be
misrepresentation.

*Resolution.* Complete proofs are written out rather than cited (Thms 5.1, 5.2),
but the manuscript states in the abstract, the introduction and §3 that the
inequalities are classical. The contribution is narrowed and labelled as
**characterisation + certification + exhaustive exact verification + an
explicit sharpness family**.

*Status.* **Resolved.**

## R1.2 (Critical) — A bug in `boundary()`

*Objection.* The implementation added a half self-loop to the boundary instead
of removing it, so every loop-graph verification was wrong.

*Resolution.* Found during verification, fixed, and **all** loop-graph results
were recomputed from scratch. Cross-checked against an independent
edge-by-edge definition on 1 972 548 (graph, set) pairs (0 mismatches) and
against the identity `bnd(S) = vol(S) - 2E(S) - loops(S)` (0 mismatches). The
affected runs (`explore_02`, `explore_04`) produced the *false* claim that
Buser fails for loop graphs; that claim is retracted.

*Status.* **Resolved** (code fixed, results recomputed, claim retracted).

## R1.3 (Critical) — Sign convention error in a verifier

*Objection.* `verify_08_after_fix.py` initially counted a violation when the
Buser slack was *negative*. Buser reads `lam2 <= 2h + 2h^2`, so a violation is a
*positive* slack. The first run after the boundary fix consequently reported
"771 violations" where there were none.

*Resolution.* Fixed; re-run reports 0 violations. The earlier "771 violations"
output is void and must not be quoted.

*Status.* **Resolved.**

## R1.4 (Major) — Is `h` the right normalisation?

*Objection.* Chung's classical statements use the Cheeger constant in a
particular normalisation. If the manuscript's `h` differs by a factor, the
bounds `2h + 2h^2` and `1 - sqrt(h(2-h))` are simply false, and the
"verification" would be testing the wrong things.

*Resolution.* Checked by construction rather than by assertion. For a
`d`-regular loop-free graph, `vol(S) = d|S|` and `h = min |E(S,S^c)|/(d|S|)`,
which is exactly the normalisation in which the classical statements are
written; the identity `2h + 2h^2 = 2h(1+h)` and `1 - sqrt(h(2-h)) =
1 - sqrt(2h(1-h/2))` are stated in the manuscript so the normalisation is
visible to a reader. The equality cases are additionally checked: the single
edge gives `h = 1`, `lam2 = 2`, and `2 = 1 + sqrt(1) <= 2h + 2h^2 = 4` — a
genuine test that the formula is not vacuous.

*Status.* **Resolved**, with the caveat in R1.5.

## R1.5 (Major) — Are the proofs complete?

*Objection.* The co-area (two-point) argument is short but error-prone; a
sloppy proof of Buser is easy to publish.

*Resolution.* §5 writes the argument in full: the layer-cake identity, the
median split, the two Cauchy–Schwarz steps, and the assembly. Each inequality
used is named. Loops are handled explicitly: a loop contributes
`w_vv * x_v^2` to the Dirichlet form, and the sweep bound
`bnd(S_t) <= (max-min) * vol(boundary)` is stated with loops accounted for.

*Status.* **Resolved** as far as is checkable by me; a referee should still
check it. I re-derived rather than cited, which is the main safeguard.

## R1.6 (Major) — Sharpness: is `lam2/h = 2` really an *infinite* family?

*Objection.* The first version of `sharpness_family()` computed
`lam2/h = 4/3, 1, 4/5, 2/3, 4/7` for `w = 1..5`, contradicting the claim that
the ratio is 2. If the sharpness family does not attain the extremum, Thm 4.2
collapses.

*Resolution.* The contradiction was real and the docstring was wrong, not the
family. Re-deriving symbolically: for `W = [[a,w],[w,a]]` the degrees are
`a+w`, `lam2 = 2w/(a+w)` and `h = w/(a+w)`, so `lam2/h = 2` identically for all
`a, w > 0` (`research/explore_07_symbolic.py`, exact sympy output). The earlier
numbers came from the **buggy** `boundary`, which inflated `h`. After the fix
the identity holds. Sharpness is retained, and the family is stated explicitly
in the paper.

*Status.* **Resolved.**

## R1.7 (Major) — "Exhaustive verification" is only up to `n = 6`

*Objection.* Exhaustive verification on small graphs cannot establish a
theorem; presenting it as evidence for the theorem's truth overstates it.

*Resolution.* Agreed, and the manuscript says so explicitly (§6 and §7): the
computation can *refute* the statement but cannot *prove* it; the proof does
that. The value of the sweep is stated as: it removes the possibility of a
hidden counterexample up to `n = 6` unweighted and `n = 4` weighted (225 349
graphs in total), and it independently validates the implementation.

*Status.* **Resolved** by rewording; the limitation is disclosed.

## R1.11 (Critical, found late) — A silent skip in the new fast path

*Objection.* When a certified high-precision eigenvalue routine was added to
make the sweep tractable, it decided connectedness by looking for eigenvalues
near **0** of `M = D^{-1/2}WD^{-1/2}`. But `L = I - M`, so the graph is
disconnected exactly when `M` has eigenvalue **1**. A zero eigenvalue of `M`
means `lam2(L) = 1`, an entirely ordinary value. The bug silently classified
those graphs as degenerate and skipped them: **8 of 864** graphs at `n = 3`,
including the uniformly self-looped complete graphs, whose exact `lam2` values
are `3/4` and `1`.

*Why it is Critical.* A verification tool that quietly excludes cases cannot
fail loudly. Had the counts not been cross-checked between the two modes, the
sweep would have reported "0 violations" over a silently reduced graph set.

*Resolution.* Fixed to test for eigenvalues near 1. Detected only because the
exact and certified modes disagreed on the number of graphs examined — which is
precisely the kind of cross-check this paper argues for. Two regression tests
were added: `test_certified_path_covers_exactly_the_same_graphs` (coverage
equality) and `test_certified_path_handles_uniform_loop_graph` (the specific
graph). After the fix both modes cover exactly 864 graphs at `n = 3` with 0
skipped, and disconnected inputs are still rejected by both.

*Status.* **Resolved**, and reported as Defect 4 in the manuscript §7.

## R1.8 (Moderate) — Does `lambda2_exact` really return exact values?

*Objection.* `lambda2_exact` calls `sp.nsimplify(val, rational=False)` when the
value is rational; `rational=False` asks for a *decimal* approximation, not an
exact rational. That is a misnomer and could introduce error.

*Resolution.* Reviewed. `nsimplify(v, rational=False)` returns the simplest
decimal representation, which for an exactly-rational eigenvalue is accurate to
machine precision but is not the exact rational. The function is therefore
renamed/documented as returning the *exact algebraic* eigenvalue via
characteristic polynomials, with rationals recovered exactly when the value is
rational (which is the case in every test). A unit test pins the two-vertex
family to exact rationals.

*Status.* **Resolved** (documented; semantics clarified). Flagged for the record
as a latent naming issue rather than a numerical error.

## R1.9 (Moderate) — Eigenvalue separation for disconnected graphs

*Objection.* `lambda2_exact` discards eigenvalues below `1e-12` and returns the
first one above; for a disconnected graph the returned value would be 0 and the
theorem would be vacuous or false.

*Resolution.* `enumerate_weighted_graphs` only yields graphs whose *unweighted
support* is connected, so disconnected graphs never enter. A unit test asserts
`lambda2_exact` returns `None` for a disconnected input.

*Status.* **Resolved.**

## R1.10 (Minor) — Random-seed reproducibility

*Objection.* Exploratory scripts use seeds; the shipped experiment does not
need randomness at all (it is exhaustive), but the README implies randomness.

*Resolution.* The shipped experiment `experiments/run_verification.py` is fully
deterministic; the seed is recorded but unused in the exhaustive path. Stated
in the README.

*Status.* **Resolved.**

---

# Review 2 — Significance and publication suitability

## R2.1 (Critical) — Is there a contribution at all?

*Objection.* Re-deriving Cheeger and Buser, with verification, is not a
research contribution. It belongs in a textbook, not a journal.

*Resolution.* This is the sharpest objection and it is **partially conceded**.
The manuscript is framed explicitly as a *characterisation and certification*
paper: the exact sharpness family `lam2/h = 2` for all `a,w > 0`, the exact
equality analysis, and an exhaustive exact-arithmetic verification pipeline are
the content. If the community does not value certified re-verification of a
classical pair, the paper will be judged incremental; this is stated in the
limitations rather than hidden. The `novelty_assessment.md` classifies the work
as **incremental** for exactly this reason.

*Status.* **Partially resolved — limitation disclosed, claim narrowed.**

## R2.2 (Major) — Does the audience care?

*Objection.* Spectral graph theory is well served by Chung's book; a journal
paper repeating it has no readership.

*Resolution.* The verifiable-computation angle is the differentiator: a reader
who wants to *trust* a Cheeger/Buser implementation gets a reference exact
pipeline. This is a real but modest niche. Journal selection (see
`journal_selection.md`) is made on that basis, preferring a journal where such
"characterisation + certification" papers are welcome.

*Status.* **Resolved** by targeting, with candid expectation setting.

## R2.3 (Major) — Is the relationship to closest prior work clear?

*Objection.* The paper must distinguish itself from Chung (1997), Mohar (1988,
1989) and Lee–Oveis Gharan–Trevisan (2014).

*Resolution.* §3 states that the two-sided bounds are classical, cites the
primary sources with DOIs resolved through Crossref, and states precisely what
is added. Because **full texts of Mohar (1988/1989) were not accessible**, the
paper avoids claiming anything about their internal contents beyond the
verified metadata, and says so.

*Status.* **Resolved**, with an explicit access caveat.

## R2.4 (Major) — Are limitations and negative results disclosed?

*Objection.* Early numerical runs suggested violations of Buser for loop graphs
and irregular graphs. If those were real, the paper would be about a counterexample
to a classical theorem — a much stronger result. They were artefacts. The paper
must not imply the counterexamples stand.

*Resolution.* The negative results are retained and presented *as artefacts of
my own implementation*, in §6 ("What did not work") and in the evidence ledger.
This is the honest framing: an early numerical claim was refuted by exact
recomputation, and the pipeline was corrected.

*Status.* **Resolved.**

## R2.5 (Moderate) — Journal fit and length

*Objection.* A 15-page paper with a small computational component may be better
suited to a shorter-format journal.

*Resolution.* Length and format left to the chosen journal's instructions;
verified during journal selection.

*Status.* **Resolved.**

## R2.6 (Moderate) — Reproducibility

*Objection.* Claimed reproducibility must be real.

*Resolution.* `README.md` documents the exact commands; the verification script
writes `results/verification_summary.json`; unit tests cover the boundary
convention, the sharpness family and the disconnected-graph guard.

*Status.* **Resolved.**

## R2.7 (Minor) — Reproducibility of the exact sweep at scale

*Objection.* The full sweep (unweighted `n <= 6`, weighted `n <= 4`) is slow;
a reader may not run it.

*Resolution.* `--quick` mode gives a fast subset; the full mode is deterministic
and reports its own elapsed time.

*Status.* **Resolved.**

---

## Summary of hostile review

| id | severity | status |
|---|---|---|
| R1.1 | Critical | Resolved — claims narrowed to characterisation/certification |
| R1.2 | Critical | Resolved — boundary bug fixed, all results recomputed |
| R1.3 | Critical | Resolved — verifier sign convention fixed; earlier output void |
| R1.11 | Critical | Resolved — silent skip in the certified fast path fixed; regression tests added |
| R1.4 | Major | Resolved — normalisation checked by construction and equality cases |
| R1.5 | Major | Resolved — proofs written out in full |
| R1.6 | Major | Resolved — sharpness family re-derived symbolically; ratio 2 identically |
| R1.7 | Major | Resolved — disclosure added; verification cannot prove |
| R1.8 | Moderate | Resolved — semantics documented |
| R1.9 | Moderate | Resolved — disconnected graphs excluded and unit-tested |
| R1.10 | Minor | Resolved — determinism documented |
| R2.1 | Critical | Partially resolved — conceded; novelty labelled incremental |
| R2.2 | Major | Resolved by journal targeting; expectations set |
| R2.3 | Major | Resolved — prior work distinguished; access caveat stated |
| R2.4 | Major | Resolved — negative results disclosed as implementation artefacts |
| R2.5 | Moderate | Resolved — format deferred to journal |
| R2.6 | Moderate | Resolved — reproducibility verified |
| R2.7 | Minor | Resolved — `--quick` provided |

**Overall.** **Four** Critical issues were found and fixed, all of them defects
in my own work rather than in the mathematics — the last one (R1.11) only after
adding a faster verification path, and it was the most dangerous because it
excluded graphs silently. One Critical objection (R2.1) is **only partially**
resolved and is carried into the limitations section: the scientific novelty of
the contribution is incremental. The manuscript is therefore written to make no
claim of a new inequality, and the final classification in the report reflects
this.