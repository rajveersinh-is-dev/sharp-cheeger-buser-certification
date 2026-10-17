# Search capability notes (2026-10-17)

## Tools available / unavailable

- `websearch` tool: **returns no results for any query**. Unavailable in this environment.
- `webfetch` against `http://export.arxiv.org/api/query`: **works**. Supports `all:`, `ti:`, `abs:` fields,
  `sortBy=submittedDate&sortOrder=descending`, pagination via `start`. Quoted phrases work.
- `webfetch` against `https://api.crossref.org/works`: **works**. Returns real bibliographic records
  (title, authors, container-title, year, DOI, page). Usable for verifying citations.

**Consequence for the integrity rules:** the novelty assessment in this project is therefore based on
arXiv full-text-field search + Crossref bibliographic search + forward-checks of DOIs, NOT on a general
web search engine. This is a genuine and material limitation on how strongly "not previously published"
can be claimed, and it is disclosed in `literature_review.md` and in the manuscript's limitations.

## Consequence for candidate generation

Candidate generation (Section 4A of the brief) had to be driven by first-principles reasoning about
where provable-and-verifiable gaps plausibly exist, then checked against the searchable indexes, rather
than the reverse (survey the literature and pick off an advertised open problem).

## Verified environment

- `tectonic` (LaTeX engine) present at `C:\Users\davea\.cargo\bin\tectonic.exe`. Smoke test compiled
  a minimal article to PDF successfully (4.3 KiB). A `Fontconfig` warning is emitted but is benign and
  does not prevent PDF output.
- No `pdflatex`, `xelatex`, `lualatex`, or `latexmk` on PATH.
- Python 3.14.6 with numpy 2.4.6, scipy 1.18.0, sympy 1.14.0, matplotlib 3.11.1.
- `git` 2.53.0 and `gh` 2.96.0 available; `gh` authenticated as `daveaddams91-dev` (active).