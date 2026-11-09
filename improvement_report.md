# Improvement Report

## A. Repository and baseline

- Repository and branch: `sharp-cheeger-buser-certification` on the current branch `jules-15988787811602642106-5802ae88` (baseline).
- Starting commit and resulting commit or PR: Evaluated against the latest commit; no new PR generated.
- Scope examined: Mathematical routines in `src/spectral_expansion.py`, `src/sharp_constants.py`, the test suite in `tests/`, and verification scripts in `experiments/` and `research/`. The memory constraints (regarding `charpoly().all_roots()` over `.eigenvals()` for numeric matrices, and random walk Laplacian formulation to avoid square roots) were also audited.

## B. Improvements implemented

No safe code change justified. The repository is already adhering tightly to the domain instructions, memory constraints, and high correctness standards.

- The numeric eigenvalue extraction explicitly and exclusively uses `.charpoly(x).all_roots()` in `lambda2_exact(W)`, avoiding the hanging issue reported for `.eigenvals()`.
- The `.eigenvals()` method is used only for purely symbolic matrices in the `research/` scripts, which perfectly complies with the directive.
- The random-walk Laplacian formulation `I - D^{-1} W` is properly leveraged in `lambda2_exact` and across relevant research scripts.
- The test suite is fast, completely green, and robust, with line coverage of 100% on the core module (`src/spectral_expansion.py`) and constants (`src/sharp_constants.py`).

## C. Verification results

| Check or command | Actual outcome | Relevant output or evidence | Explanation of limitation |
|---|---|---|---|
| `pytest tests/` | PASS | `47 passed in 36.89s` | Evaluated locally; tests fast and reliable. |
| `pytest --cov=src tests/` | PASS | `TOTAL 228 0 100%` coverage | Verified locally. Full statement coverage on `src/`. |
| `python experiments/run_verification.py --quick` | PASS | `total graphs verified: 601`, `total violations: 0` | The fast verification script executed cleanly. |
| Code audit (Memory directives) | PASS | `grep` showed `charpoly(x).all_roots()` in numeric paths. | Found complete compliance. |

## D. Measurable impact

- No before/after measurements are applicable as no code changes were introduced.

## E. Remaining risks

- Exhaustive verifications are inherently limited to a finite domain (e.g. `n <= 6`). The manuscript and limitations file correctly disclose this.
- Outstanding pre-submission requirements (e.g., ORCID iD) identified in `research/limitations.md` must be fulfilled by the human researcher before publication.

## F. Next-run priorities

1. Re-evaluate if new edge cases or symbolic configurations are introduced to the repository.
2. Consider expanding the domain of exact verification if a more efficient computational infrastructure is provided.
3. Keep dependency files (`requirements.txt`) up to date as library releases (`sympy`, `mpmath`) evolve.

## G. Overall verdict

- AUDIT COMPLETED — NO SAFE CODE CHANGE JUSTIFIED
