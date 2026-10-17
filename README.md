# Sharp Cheeger--Buser characterisation and certification

Exact characterisation, certified computation and exhaustive verification of the
two-sided Cheeger--Buser comparison between the second eigenvalue of the
normalised Laplacian and the Cheeger constant of a finite weighted graph.

**Author:** Rajveersinh Vishal Pardeshi (Independent Researcher)
**Manuscript:** [`paper/main.tex`](paper/main.tex) · compiled PDF: `paper/main.pdf`
**Target journal:** The Quarterly Journal of Mathematics (OUP), standard-licence
route (no mandatory author publication fee).

---

## Research question

For a finite connected **weighted graph** `W` (symmetric, non-negative, loops
allowed) with

```
d_v    = sum_u W[v][u]                       (a self-loop adds once)
vol(S) = sum_{v in S} d_v
bnd(S) = sum_{v in S, u not in S} W[v][u]    (a loop inside S contributes nothing)
h(G)  = min { bnd(S)/vol(S) : 0 < vol(S) <= vol(V)/2 }   in (0,1]
```

and the generalised normalised Laplacian `L = I - D^{-1/2} W D^{-1/2}`, how do
`lambda_2(L)` and `h` control one another?

The classical answer, quoted for *regular, loop-free* graphs, is

```
1 - sqrt( h(2 - h) )   <=   lambda_2   <=   2h + 2h^2            (Cheeger / Buser)
```

## Contribution (and what is *not* claimed)

**Not claimed:** that these inequalities are new. They are classical
(Mohar 1988/1989; Alon-Milman 1985; Dodziuk 1984; Chung 1997). See
`research/novelty_assessment.md`.

**Claimed:**

1. A self-contained statement valid for arbitrary finite weighted graphs **with
   self-loops**, with hypotheses written out precisely enough to be true.
2. **Exact sharpness:** for `W = [[a, w], [w, a]]` with `a, w > 0`,
   `lambda_2 = 2w/(a+w)` and `h = w/(a+w)`, hence `lambda_2 / h = 2`
   *identically*. The constant `2` in Buser's bound is therefore attained and
   cannot be lowered.
3. A **verification pipeline** using exact rational arithmetic and exact
   characteristic polynomials that checks both inequalities on **every** connected
   graph with `n <= 6` vertices and on weighted graphs with self-loops with
   `n <= 4`.
4. A documented account of **three defects** in our own implementation, one of
   which produced a false counterexample to Buser's inequality that was
   subsequently retracted.

## Status of findings — read this before citing anything

| Statement | Status |
|---|---|
| The two-sided inequality holds in the stated generality | **Verified**: proved in the manuscript, and exhaustively checked over the ranges above with **0 violations** |
| The constant `2` is optimal | **Verified**: proved exactly (Theorem 3.1) |
| Exact sweep over all graphs `n <= 6` | **Verified**: `results/verification_summary.json` |
| Novelty of the *inequalities* | **Classical** — not claimed as new |
| Novelty of the overall package | **Incremental** — characterisation + certification |
| Literature completeness | **Partial**: general web search was unavailable; see `research/literature_review.md` |

An exhaustive sweep can **refute** a universal statement; it cannot **prove** one.
The proof does that. The computation is evidence about implementations and about
small cases.

## Layout

```
paper/          main.tex, references.bib, main.pdf (compiled), README.md
src/            spectral_expansion.py   exact Cheeger constant, exact eigenvalues
                sharp_constants.py      certified extremal coefficient
tests/          pytest suite (19 tests) covering conventions, sharpness, guards
experiments/    run_verification.py     the exhaustive exact sweep
results/        verification_summary.json  (generated; committed)
research/       literature_review.md, novelty_assessment.md, evidence_ledger.md,
                hostile_review.md, journal_selection.md, specification.md,
                acceptance_assessment.md, limitations.md,
                github_account_verification.md, exploration scripts
```

## Reproducing the results

Requires Python 3.11+ with `sympy`, `numpy`, `mpmath`, `pytest`.

```bash
pip install -r requirements.txt

# unit tests (fast, ~40 s)
python -m pytest tests/ -q

# exhaustive exact verification, fast subset (~10 s)
python experiments/run_verification.py --quick

# full sweep: all connected graphs n <= 6; weighted+loops n <= 4 (slow)
python experiments/run_verification.py
```

The verification is **deterministic**: no random seeds, no network access, no
external data. Repeated runs produce byte-identical output apart from the
`elapsed_seconds` field.

Expected outcome of the full sweep:

```
total exact graphs verified : 27475 (unweighted) + weighted-with-loops
total violations            : 0
max lam2/h                  : 2.0   (attained exactly, Theorem 3.1)
```

### Building the PDF

```bash
tectonic -X compile paper/main.tex --outdir paper
```

`tectonic` fetches the standard TeX bundle on first use; no local TeX
installation is needed.

## Data provenance and licensing

* **No data is used.** The "experiments" are exhaustive enumerations of finitely
  many small graphs generated by the code itself.
* **No external dataset, no download, no network access** at run time.
* Code and documentation: MIT (see `LICENSE`). The manuscript is the author's
  own original text and is released under the same licence.
* The bibliography contains verified bibliographic **metadata** only; no
  copyrighted text from any cited work is reproduced.

## Known limitations

1. The inequalities are classical; novelty is incremental.
2. Exhaustive verification is finite and can only refute, not prove.
3. The proofs use the standard co-area method; they are written out in full but
   the method itself is classical.
4. The literature review is restricted to arXiv and Crossref because general web
   search was unavailable.
5. Self-loop conventions are a choice made here (loop adds `W_vv` to the degree
   and contributes `W_vv f_v^2` to the energy). Other conventions exist.

## Citation

If you use this code or the characterisation, cite the manuscript:

```bibtex
@techreport{pardeshi2026sharp,
  author      = {Rajveersinh Vishal Pardeshi},
  title       = {Exact characterisation and certification of the sharp
                 two-sided {Cheeger}--{Buser} bound for finite weighted graphs},
  institution = {Independent Researcher},
  year        = {2026},
  note        = {Manuscript; not submitted for publication at the time of writing}
}
```

`CITATION.cff` carries the same information in CFF form.

## Licence

MIT — see [`LICENSE`](LICENSE).