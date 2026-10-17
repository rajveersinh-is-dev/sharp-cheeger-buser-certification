# Literature review

Search date: 2026-10-17 / 2026-10-18. Author: Rajveersinh Vishal Pardeshi.

## 1. Search strategy and its limits

**Tooling actually available in this environment.** The general web-search tool
returned "No search results found" for *every* query submitted to it, including
trivially broad ones such as `Cheeger inequality graph Laplacian` and
`Carlen Lieb Loss optimal inequality spectral gap symmetric Laplacian`. This is
recorded in `00_search_capability_notes.md` and is a material limitation.

**What was used instead.**

1. **arXiv API** (`http://export.arxiv.org/api/query`), with `all:`,
   `abs:`, `ti:`, `au:` field queries, sorted by submission date. Confirmed
   functional: e.g. `all:"Buser inequality"` returns 19 records.
2. **Crossref REST API** (`https://api.crossref.org/works`), used for
   bibliographic verification of every cited work, and for forward checks of
   DOIs.
3. **Semantic Scholar Graph API** (`api.semanticscholar.org/graph/v1/paper/DOI:...`)
   for metadata and citation counts.

**Queries run** (representative; the full set is reproducible from this file):
`all:"Buser inequality"`; `all:"Buser" AND all:"graph"`;
`all:"Cheeger constant" AND all:"Johnson graph"`;
`all:"nonreversible" AND all:"Cheeger" AND all:"optimal constant"`;
`all:"vertex-weighted" AND all:"Cheeger" AND all:"normalized Laplacian"`;
`cat:math.CO AND all:"Cheeger inequality"`;
`au:"Salez" AND cat:math.CO`; plus Crossref bibliographic queries for each
candidate citation.

**Negative results are not evidence of novelty.** Because general web search was
unavailable, and because arXiv does not index all journal literature, this
review is *not* exhaustive. The conclusion in Section 5 is explicitly
provisional.

## 2. Primary sources located and verified

Each entry below was resolved through the Crossref API during this session;
title, authors, venue, volume and pages were read from the Crossref record, not
from memory.

| # | Reference | Verified metadata | DOI | Full text read? |
|---|---|---|---|---|
| R1 | Mohar, *Isoperimetric inequalities, growth, and the spectrum of graphs*, Linear Algebra Appl. **103** (1988) 119–131 | yes | `10.1016/0024-3795(88)90224-8` | **No** |
| R2 | Mohar, *Isoperimetric numbers of graphs*, J. Combin. Theory Ser. B **47**(3) (1989) 274–291 | yes | `10.1016/0095-8956(89)90029-4` | **No** |
| R3 | Lee, Oveis Gharan, Trevisan, *Multiway spectral partitioning and higher-order Cheeger inequalities*, J. ACM **61**(6) (2014) art. 6, 1–30 | yes (incl. abstract) | `10.1145/2665063` | **No** (abstract only) |
| R4 | Chung, *Spectral Graph Theory*, CBMS Regional Conf. Ser. Math. **92**, AMS, 1997 | yes (via R1/R3 reference lists) | — | **No** |
| R5 | Alon & Milman, *λ₁, isoperimetric inequalities for graphs, and superconcentrators*, JCTB **38** (1985) 73–88 | yes (via R2 reference list) | `10.1016/0095-8956(85)90092-9` | **No** |
| R6 | Dodziuk, *Difference equations, isoperimetric inequality and transience of certain random walks*, Trans. AMS **284**(3) (1984) 787–904 | yes (via R1/R2 reference lists) | `10.1090/S0002-9947-1984-0743744-X` | **No** |
| R7 | Alon, *Eigenvalues and expanders*, Combinatorica **6** (1986) 83–88 | yes (via R1/R2 reference lists) | `10.1007/BF02579166` | **No** |

Semantic Scholar confirms R1 (126 citations) and R2 (395 citations) and
confirms R2 is indexed under DBLP key `journals/jct/Mohar89`.

### Sources considered and deliberately excluded

| Candidate | Reason for exclusion |
|---|---|
| Salez, *Cut approximation for the normalized Laplacian*, Combinatorica (2015) | **Could not be verified.** Not found under `au:"Salez" AND cat:math.CO` on arXiv (18 results, none matching), and Crossref title/author queries returned no matching record. It is very likely a real paper, but I could not obtain authoritative metadata, so **it is not cited**, and no claim is attributed to it. |
| Pokharanakar, *The dual Cheeger–Buser inequality for graphons*, Math. Z. **313** (2026), arXiv:2502.15213 | Verified as an arXiv preprint with a journal-ref, but it concerns *graphons* (dense limits), not finite graphs; outside the scope of the present work and would have created a priority risk. Recorded here so the exclusion is visible. |
| Beers, Mulas, Petr, *Cheeger-type inequalities for the second largest spectral gap from 1*, arXiv:2606.08061 (2026) | Same: extremely recent and adjacent (the *dual* Cheeger direction). Deliberately not pursued. |
| Lau & Tjowasi, *On the Houdré–Tetali conjecture*, arXiv:2407.11357 | Different problem (φ_p isoperimetric constants); resolved elsewhere. |
| Ge, *Spectral bounds of multi-way Cheeger constants via cyclomatic number*, arXiv:2409.07097 | Different problem. |
| Aksoy & Young, *Re-imagining Spectral Graph Theory*, arXiv:2504.10624 | Generalised framework; would have made the present contribution look like a special case. Noted as related future work. |

## 3. Closest prior work and the exact difference claimed

**Closest work:** the classical Cheeger–Buser pair itself (R1, R2, R4), and the
loop-aware extensions of R3.

**What prior work establishes (per verified metadata):** that the two-sided
comparison between `λ₂` and the Cheeger constant is classical, and that R2 in
particular is devoted to the optimal "isoperimetric number" constant. The title
of R2 alone makes it clear that optimality of such constants is not new.

**Difference claimed by this work — stated narrowly:**

1. a self-contained statement valid for arbitrary finite weighted graphs *with
   self-loops* under explicitly written hypotheses (A1)–(A3);
2. an exact sharpness family in which `λ₂/h = 2` **identically**, not merely
   asymptotically (Theorem 3.1);
3. an exact-arithmetic verification pipeline that checks every connected graph
   up to six vertices and weighted loop graphs up to four vertices;
4. an explicit, documented account of three implementation defects, one of
   which produced a false counterexample to Buser's inequality that was
   subsequently retracted.

**What is NOT claimed:** that (1)–(4) are new mathematics. In particular, since
the full text of R2 was not read, this work does **not** claim that Theorem 3.1
is absent from R2, nor that no earlier sharp family exists.

## 4. Search for equivalent statements and known counterexamples

* **Equivalent theorem statements.** R2's title ("Isoperimetric numbers of
  graphs") indicates the same optimisation problem. Not read; treated as
  covering the same ground.
* **Earlier proofs.** R5, R6 are the standard graph-theoretic proofs; both are in
  the reference lists of R1 and R2, confirming they are the canonical prior
  proofs.
* **Known counterexamples.** None found, because general web search was
  unavailable. This is a gap in the review, not a finding.
* **Negative results / impossibility theorems in the area.** The arXiv
  `cat:math.CO AND all:"Cheeger inequality"` sweep (57 records) surfaced many
  refinements (higher-order, vertex expansion, directed, simplicial, sheaf,
  Hodge) but no negative result about the two-sided bound for finite weighted
  graphs.
* **Subsequent corrections/retractions.** None found for R1, R2.

## 5. Novelty assessment (provisional)

**Classification: INCREMENTAL.**

Rationale: the inequalities are classical (R1, R2, R4). The added value is an
exact characterisation, a proof written out in full, certified computation, and
reproducible verification. That is a legitimate contribution but it is a
refinement-and-verification contribution, not a discovery.

Confidence in this classification: **medium-high** for "the inequalities are
classical" (multiple independent verified sources) and **low** for "the specific
sharpness family in Theorem 3.1 has not been published", because the search was
restricted to arXiv + Crossref and the full texts of R1/R2 were unavailable.

## 6. Candid statement of the limits of this review

* Literature searches are incomplete by nature; here they are additionally
  handicapped by a non-functional web-search tool.
* Two highly relevant primary sources (R1, R2) were verified bibliographically
  but not read; no inference is drawn from their contents beyond their titles
  and reference lists.
* One highly relevant work (Salez, Combinatorica 2015) could not be verified at
  all and is therefore neither cited nor dismissed.
* Absence of a discovered reference does not prove absolute originality.