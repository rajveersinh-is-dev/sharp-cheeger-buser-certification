# Novelty assessment

Date: 2026-10-18. Author: Rajveersinh Vishal Pardeshi.

## Verdict

**Classification: Incremental.**

Not "strongly supported originality", and not "potentially novel". The
mathematical content that the paper delivers is classical; what is added is
exact characterisation, self-contained proof, certified computation, and
reproducible verification.

## Why not a higher classification

The two-sided bound
\[
  1-\sqrt{h(2-h)} \;\le\; \lambda_2 \;\le\; 2h+2h^2
\]
is classical. This is established from multiple independently verified
bibliographic sources, and the title of Mohar (1989), *Isoperimetric numbers of
graphs* (JCTB 47:274–291, DOI `10.1016/0095-8956(89)90029-4`), indicates
directly that the *optimality* of such constants is long-studied territory. A
paper whose headline result reproduces that pair cannot claim discovery.

## What is actually added

| # | Addition | Genuinely new? |
|---|---|---|
| 1 | Statement of the two-sided bound for arbitrary finite weighted graphs **with self-loops**, with hypotheses (A1)–(A3) written out precisely enough that it is true | Restatement + correct generalisation of scope. The literature does handle loops (Lee–Oveis Gharan–Trevisan 2014) but via a *different* Cheeger constant; the manuscript uses the plain one and proves it directly. Marginal novelty. |
| 2 | **Exact sharpness family**: for `W = [[a,w],[w,a]]`, `λ₂/h = 2` *identically* for all `a,w>0` | Cannot claim novelty — the full text of Mohar (1989) was not read, and this is exactly the kind of statement such a paper contains. Claimed as an explicit construction, not as a discovery. |
| 3 | **Exact-arithmetic verification pipeline** (rational `h`, exact characteristic-polynomial `λ₂`) covering every connected graph with `n≤6` and weighted loop graphs with `n≤4` | Genuinely new as *engineering/verification* content. Not a mathematical theorem. |
| 4 | **Documented three-defect post-mortem**, including a retracted false counterexample | Genuinely new and, in the author's judgement, the most transferable content. |

## The retracted result — an important honesty note

Early in the project, floating-point experiments appeared to show that Buser's
inequality `λ₂ ≤ 2h + 2h²` **fails** for graphs with self-loops and for
irregular graphs. Had that survived, it would have been a far stronger and more
interesting result — a counterexample to a classical statement.

It did not survive. Three defects in my own code were responsible:

1. `boundary()` added half of each internal self-loop to the boundary instead
   of excluding it (invisible on loop-free graphs);
2. the verifier used the wrong sign for "slack", reporting violations of a
   bound that holds;
3. the Cheeger minimum was restricted to sets containing a fixed vertex, which
   is invalid unless the graph is regular.

After fixing all three and recomputing exactly, the number of violations is
**zero** across the verified range. The apparent counterexample is an artefact.

This is recorded prominently in the manuscript (§7), in the evidence ledger
(item C4), and here, because a reader deserves to know that a tempting
"discovery" was withdrawn.

## Confidence and residual uncertainty

| Claim | Confidence | Why |
|---|---|---|
| The two-sided inequalities are classical | **High** | Multiple verified sources (Mohar 1988/1989, Alon–Milman 1985, Dodziuk 1984, Chung 1997) |
| The constant 2 in Buser is optimal | **High** | Provable here from Theorem 3.1 + Theorem 5.2 |
| The *specific* sharpness family is new | **Low** | Mohar (1989) full text unread; this is precisely the content such a paper contains |
| No stronger result is known for weighted graphs with self-loops | **Low** | Web search unavailable; arXiv/Crossref only |
| The verification has no bugs | **Medium-high** | Two independent cross-checks of `boundary`, exact arithmetic, 19 unit tests, and the pipeline previously caught three real bugs |

## How a referee should read this paper

As a careful, verifiable, well-documented *characterisation and certification*
of a classical pair — valuable as a trustworthy reference implementation and as
an honesty case study in computer-assisted verification — but not as a
theoretical advance. The paper says so in its abstract, its limitations section
and its conclusion, and this assessment should not be read as more favourable.