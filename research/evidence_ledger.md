# Evidence ledger

Every claim made in the manuscript that carries weight appears here with its
verification status.  Status vocabulary: **Verified**, **Partially verified**,
**Unverified**, **Contradicted**, **Not applicable**.

Author: Rajveersinh Vishal Pardeshi (Independent Researcher). Ledger opened
before the manuscript was written and updated as results arrived.

---

## C1 — Definitions

| field | value |
|---|---|
| Claim | For a weighted graph `W` (symmetric, `>=0`, `W_vv` = self-loop weight), with `d_v = sum_u W_vu`, `vol(S) = sum_{v in S} d_v`, `bnd(S) = sum_{v in S, u not in S} W_vu`, and `h = min{bnd(S)/vol(S) : 0 < vol(S) <= vol(V)/2}`, these are the conventions used throughout. |
| Type | definition |
| Source | standard convention; see Chung, *Spectral Graph Theory*, CBMS 92 (1997) for the `h` used in the classical statements |
| Verification | identity `bnd(S) = vol(S) - 2E(S) - loops(S)` checked programmatically on 1 972 548 (graph, set) pairs, 0 mismatches (`research/verify_08_after_fix.py`, step 2) |
| Status | **Verified** |
| Used in | Abstract; §2; Prop. 2.1 |

## C2 — `boundary` implementation correctness

| field | value |
|---|---|
| Claim | The implementation of `bnd(S)` equals the independent definition "sum of weights of edges with exactly one endpoint in `S`". |
| Type | implementation check |
| Verification | `research/verify_08_after_fix.py` step 1: 1 972 548 pairs, 0 mismatches; also a unit test |
| Status | **Verified** |
| Counterevidence | **An earlier version of `boundary()` incorrectly ADDED a half self-loop and therefore returned a wrong value for every graph containing a self-loop.** This was found during verification, fixed, and all loop-graph results were recomputed. The bug was invisible on loop-free graphs, which is why early exploratory numbers looked self-consistent. |
| Used in | everything involving self-loops |

## C3 — The two-vertex loop family

| field | value |
|---|---|
| Claim | For `W = [[a,w],[w,a]]` (`a` = loop weight, `w` = edge weight, `a,w > 0`): `lam2 = 2w/(a+w)`, `h = w/(a+w)`, hence `lam2/h = 2` for **all** `a,w>0`. |
| Type | mathematical result (exact) |
| Verification | sympy symbolic derivation, `research/explore_07_symbolic.py`; also `sharp_constants.sharpness_family()` |
| Status | **Verified** |
| Note | This corrects an earlier, wrong claim in `sharp_constants.sharpness_family`'s docstring that the ratio `2` is attained only at `w = 0`; the symbolic run shows `lam2/h = 2` identically. The docstring was rewritten. |
| Used in | §3 (sharpness family); Thm 4.2 |

## C4 — Buser upper bound: validity in the generality stated

| field | value |
|---|---|
| Claim | `lam2 <= 2h + 2h^2` for every finite connected weighted graph (self-loops allowed) with the plain `h` above. |
| Type | mathematical result |
| Support | classical Buser/Dodziuk argument; re-derived in the manuscript §5 |
| Verification | exhaustive EXACT check: unweighted `n <= 5` (771 graphs) and weighted-with-loops `n <= 3` (876 graphs), **0 violations**; full run over unweighted `n <= 6` and weighted `n <= 4` recorded in `results/verification_summary.json` |
| Status | **Verified** (small-`n` exhaustive) / **Partially verified** (all `n`, by proof only) |
| Known counterevidence | Early numerical runs (`research/explore_02.py`) *appeared* to show violations of Buser for irregular graphs. Those runs used an incorrect `h` normalisation (an edge/volume constant rather than the vertex/volume constant used in the classical statement) and a buggy `boundary`. After both were corrected, no violation survived. **The early "violation" was an artefact of my own code, not a property of the inequality.** |
| Used in | Abstract; Thm 5.1; §5 |

## C5 — Cheeger lower bound: validity in the generality stated

| field | value |
|---|---|
| Claim | `lam2 >= 1 - sqrt(h(2-h))` (= `1 - sqrt(2h(1-h/2))`) under the same hypotheses. |
| Type | mathematical result |
| Verification | same exhaustive exact sweep: **0 violations** in 1 647 graphs (and the larger run) |
| Status | **Verified** (small-`n` exhaustive) / **Partially verified** (all `n`, by proof only) |
| Counterevidence | `research/explore_04_exact.py` reported 264 "Cheeger violations" using a *neighbourhood* (loop-aware) constant. Those were violations of a different inequality (with `h_loop`), not of the stated theorem. The manuscript therefore reports both constants and explains which is which, rather than quoting the larger violation count. |
| Used in | Abstract; Thm 5.2; §5 |

## C6 — Optimality of the constant 2

| field | value |
|---|---|
| Claim | The constant `2` in `lam2 <= 2h + 2h^2` cannot be lowered: `sup_G lam2/h = 2`. |
| Type | mathematical result |
| Proof | lower bound `lam2 <= 2h + 2h^2` implies `lam2/h <= 2 + 2h`, and the two-point ratio computation gives the leading coefficient `2`; upper bound realised by the family of C3, where `lam2/h = 2` identically. |
| Verification | exact computation of `lam2/h = 2` on the family (C3); certified interval arithmetic in `src/sharp_constants.py` returns `[2,2]`; exhaustive max over 1 647 exact graphs equals exactly `2.000000...` (attained at the single edge) |
| Status | **Verified** |
| Used in | Abstract; Thm 4.2; §4 |

## C7 — Literature fact: prior existence of the two-sided inequalities

| field | value |
|---|---|
| Claim | The two-sided comparison between `lam2` and the Cheeger constant is classical, not new. |
| Type | literature fact |
| Source | Chung, *Spectral Graph Theory*, CBMS 92 (1997); Mohar, *Linear Algebra Appl.* **103** (1988) 119–131, DOI `10.1016/0024-3795(88)90224-8`; Mohar, *J. Combin. Theory Ser. B* **47** (1989) 274–291, DOI `10.1016/0095-8956(89)90029-4`; Lee, Oveis Gharan, Trevisan, *J. ACM* **61**(6) (2014) art. 6, DOI `10.1145/2665063` |
| Verification | each DOI resolved through the Crossref REST API during this session; titles, authors, venues, volumes and page ranges read from the Crossref records, not from memory |
| Status | **Verified** (bibliographic metadata) / **Unverified** (full texts not read) |
| Limitation | I read only the Crossref metadata and the abstract of the JACM paper, **not** the full texts of Mohar (1988/1989). Any claim about what those papers do or do not contain beyond the metadata is **not** made. |
| Used in | §1 Introduction; §3 Related work; `novelty_assessment.md` |

## C8 — Literature fact: self-loop / loop-aware Cheeger constants

| field | value |
|---|---|
| Claim | A "loop-aware" Cheeger constant involving the neighbourhood volume is used in the literature to handle self-loops. |
| Type | literature fact |
| Source | Lee, Oveis Gharan, Trevisan, *J. ACM* **61**(6) (2014), DOI `10.1145/2665063` (abstract and reference list inspected via Crossref) |
| Verification | abstract confirmed via Crossref; the paper's own reference list contains Cheeger (1970), Chung (1996/1997), Dodziuk (1984), Alon–Milman (1985) |
| Status | **Partially verified** |
| Limitation | I did **not** obtain the full text, so I do not assert the exact formula used there. The manuscript defines its constants self-containedly and cites the paper only for the fact that loop-aware variants are in use. |
| Used in | §2; §3 |

## C9 — Novelty of the contribution

| field | value |
|---|---|
| Claim | The contribution is an **exact, certified and exhaustively verified characterisation of the sharp two-sided relation**, not a new inequality. |
| Type | interpretation |
| Basis | C4–C6 |
| Status | **Verified** as a description of what the paper actually does |
| Residual uncertainty | Whether the *packaging* (exact certification + exhaustive verification + sharpness family stated explicitly) is publishable novelty is a matter for referees; it is **not** claimed to be a new inequality. |
| Used in | Abstract; §1; §7 Limitations; `novelty_assessment.md` |

## C10 — Search capability limitation

| field | value |
|---|---|
| Claim | The novelty assessment rests on arXiv and Crossref searches only; general web search was unavailable in this environment. |
| Type | limitation |
| Verification | every query issued through the `websearch` tool returned "No search results found", including trivially broad queries such as `Cheeger inequality graph Laplacian` |
| Status | **Verified** |
| Used in | `literature_review.md`; §7 Limitations; final report |

---

## Claims explicitly NOT made

* Not claimed: that the two-sided inequality is new.
* Not claimed: that the constant `2` is new (it is Buser's).
* Not claimed: that Mohar (1989) does or does not contain a particular result — full text not read.
* Not claimed: any statement about journals' policies that was not read on the official page during this session.
* Not claimed: any experimental/measurement result that was not produced by code in this repository.