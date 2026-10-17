# Building the manuscript

## Prerequisites

No local TeX installation is required. The project is built with
[Tectonic](https://tectonic-typesetting.github.io/), which downloads the TeX
bundle itself on first use.

Install Tectonic with, e.g.:

```bash
cargo install tectonic        # or: winget install Tectonic.Tectonic
```

## Build

From the **repository root**:

```bash
tectonic -X compile paper/main.tex --outdir paper
```

This produces `paper/main.pdf` (8 pages) and runs BibTeX automatically.

To keep logs for inspection:

```bash
tectonic -X compile paper/main.tex --outdir paper --keep-logs
```

## Files

| File | Role |
|---|---|
| `main.tex` | the manuscript |
| `references.bib` | bibliography; every entry verified via the Crossref API |
| `main.pdf` | compiled output (tracked in the repository) |
| `main.log`, `main.aux`, `main.blg` | intermediates; **not** tracked (see `.gitignore`) |

## Verified build

The committed `main.pdf` was produced by this command in this repository with

* `tectonic` (XeTeX-based) on Windows,
* 8 pages of output,
* **0 undefined citations and 0 undefined references** in `main.log`,
* one cosmetic `Overfull \hbox` warning in the data-availability section.

Known warning (benign):

```
warning: main.tex:365: Overfull \hbox (64.26088pt too wide) in paragraph at lines 362--365
```

## Citation policy in `references.bib`

Every entry was resolved through the Crossref REST API (title, authors, venue,
volume, issue, pages, DOI read from the authoritative record). Entries whose
full text could not be accessed are used only for claims supported by their
verified metadata. No reference was reconstructed from memory. See
`research/literature_review.md` for the full query log.