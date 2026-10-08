# Draft results section

This folder contains the separate LaTeX results section requested by the
collaborator. It reports the completed historical campaign scans and the
comparison using identical target dates across lags. The section includes
five tables, four individual component figures, descriptive extrema, and
a limited astrophysical comparison with Lu and Xi.

This is a descriptive draft frozen for collaborator review. Initial categorical calibration and the ordering audit are complete. Realistic
observation calibration, flux uncertainty propagation, significance testing of the fixed
support curves, and expanded scientific runs under equation 3.18 remain
pending. The historical observations do not supply a complete supported
curve for every requested target under the retained lag range and guard.
The original full manuscript remains a separate previous study.

## Files

`results_section.tex` is the section to integrate into a manuscript.
`preview.tex` supplies a minimal document wrapper for review.
`references.tex` supplies the six references used by the preview.
The five table files and four images under `figures` are frozen exports
from the completed scans. `descriptive_extrema.csv` records maxima for
all four retained curves. `asset_manifest.json` identifies the numerical
inputs and exported asset hashes. `delivery_manifest.json` records the
delivered PDF and document source hashes.

The compiled preview is at `output/pdf/surd_campaign_results.pdf` relative
to the repository root. It has ten pages. All pages were rendered and
visually inspected, and the final compilation has no layout or unresolved
reference warnings. Numerical export checks verify the full lag grids,
invariant target dates, and normalized component sums.

## Reproduce

From the repository root, with the saved scans available:

```sh
python -m collaborator_rebuild.build_results_assets
mkdir -p tmp/pdfs/surd_section_build output/pdf
tectonic --keep-logs --outdir tmp/pdfs/surd_section_build collaborator_rebuild/results_section/preview.tex
cp tmp/pdfs/surd_section_build/preview.pdf output/pdf/surd_campaign_results.pdf
```

To compile a copied folder by itself, run `tectonic preview.tex` inside it.
The preview defines its asset path as the current document directory.

## Manuscript integration

The host document needs `graphicx`, `amsmath`, `booktabs`, `float`, and
`natbib`. Include `results_section.tex` and keep its tables and figures
together. Set `\SURDResultsPath` to this folder relative to the compiling
document. Its default is `collaborator_rebuild/results_section`.
Add the six `CR` citation keys to the host bibliography; the preview's
`references.tex` records their metadata. Do not insert a second bibliography
into an existing manuscript. The campaign definition refers to Peterson
et al. 2004, ApJ 613, 682, DOI 10.1086/423269.

The preview disables automatic word hyphenation. The same settings can be
applied in a host preamble with `\hyphenpenalty=10000` and
`\exhyphenpenalty=10000`.

The algorithmic interpretation and astrophysical interpretation still require
the scientific review specified by the collaborator. This draft is not a
claim that every requested analysis has been completed.
