# Limitations

An honest, consolidated list. Each item states what is *not* claimed and why.

## Scientific

**L1 — The inequalities are classical.**
`1 - sqrt(h(2-h)) <= lambda_2 <= 2h + 2h^2` is Cheeger's and Buser's bound.
This work does not claim it as a new theorem. Verified prior art: Mohar 1988
(LAA 103:119-131, DOI `10.1016/0024-3795(88)90224-8`) and Mohar 1989 (JCTB
47:274-291, DOI `10.1016/0095-8956(89)90029-4`), plus Alon-Milman 1985,
Dodziuk 1984, Chung 1997. The title of Mohar (1989) — *Isoperimetric numbers
of graphs* — indicates directly that optimising such constants is long-studied
ground.

**L2 — The novelty that remains is incremental.**
What is genuinely added: an exact statement for weighted graphs with self-loops;
an explicit sharpness family with `lambda_2/h = 2` identically; an
exact-arithmetic verification pipeline; and a documented three-defect
post-mortem. A referee may reasonably judge that insufficient for a research
article. We state that plainly rather than overstate it.

**L3 — The sharpness family is very likely not new.**
The full texts of Mohar (1988) and Mohar (1989) were not accessible in this
environment. Theorem 3.1 of the manuscript is precisely the kind of statement
such a paper contains. **We therefore do not claim it is new.**

**L4 — Exhaustive verification can refute, not prove.**
The sweep reaches `n <= 6` (unweighted, 27 475 graphs) and `n <= 4` (weighted
with self-loops, 197 874 graphs), covering **225 349 graphs with zero
violations**. It
excludes counterexamples up to those sizes and validates the implementation; it
does not establish the theorem for all graphs. The proof does that.

**L5 — The method is classical.**
The co-area / two-point reduction is standard. The manuscript writes it out in
full rather than citing it, so that a reader can check every constant without
trusting a citation — but writing out a classical proof is not a contribution
to mathematics.

**L6 — Self-loop conventions are a choice.**
We count a self-loop as adding `W_vv` to the degree and contributing
`W_vv * f_v^2` to the energy, and we exclude internal loops entirely from
`bnd(S)`. Other conventions exist in the literature. Conclusions may differ
under them and we do not claim robustness across conventions.

## Methodological

**L7 — A false counterexample was produced and retracted.**
Early floating-point experiments appeared to show Buser's inequality failing for
graphs with self-loops. It did not. Four defects in our own code were
responsible: (1) `boundary()` added half of each internal self-loop instead of
excluding it; (2) the verifier used the wrong sign for "slack"; (3) the Cheeger
minimum was restricted to sets containing a fixed vertex, which is invalid
unless the graph is regular; (4) the certified fast eigenvalue path tested for
eigenvalues near 0 of `M = D^{-1/2}WD^{-1/2}` instead of near 1, silently
skipping 8 of 864 graphs at `n = 3`. All four were fixed and everything was
recomputed. The retraction is documented in the manuscript (§7), the evidence
ledger (C4, C11) and `novelty_assessment.md`.

**L8 — `lambda2_exact` uses a numeric tolerance to recognise exact zeros.**
Eigenvalues themselves come from exact characteristic polynomials, so the
computation is exact; the tolerance `1e-12` is used only to identify the
zero eigenvalue (to count connected components). For rational inputs the
returned value is the exact rational whenever the eigenvalue is rational, which
is the case in every test performed.

## Process

**L9 — Literature coverage is partial.**
General web search was unavailable: every query submitted to the search tool
returned no results, including trivial ones. The novelty assessment rests on
arXiv and Crossref only. **Absence of a retrieved reference is not proof of
originality.**

**L10 — One highly relevant work could not be verified and is not cited.**
Salez, *Cut approximation for the normalized Laplacian*, is very likely real and
highly relevant, but no authoritative record could be obtained (absent from the
arXiv author query and from Crossref title/author queries). It is therefore
neither cited nor dismissed.

**L11 — Journal policy is a snapshot.**
Fee and self-archiving statements were read from official Oxford Academic pages
on 2026-10-18. Policies change; they must be re-checked immediately before
submission.

**L12 — One proof step is summarised rather than fully derived.**
The lower-bound (Cheeger) proof isolates one standard Cauchy–Schwarz step. The
manuscript says so explicitly in a remark and gives the supporting estimate in
an appendix, rather than presenting the proof as more complete than it is.

**L13 — Outstanding pre-submission requirements I cannot satisfy myself.**
QJM's instructions (verified 2026-10-18) require an **ORCID iD at submission**.
The author has not provided one, and none is invented here. Submission cannot
proceed until that is resolved. The journal also encourages running the
Paperpal preflight tool; this is explicitly optional and has not been run.

## What would change these conclusions

* Obtaining the full text of Mohar (1989) would settle L3 either way.
* A working web search would materially improve L9.
* Extending the exact sweep to `n = 7`+ weighted graphs would strengthen L4, but
  cannot remove it — no finite sweep can.