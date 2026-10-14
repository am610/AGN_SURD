# SURD analysis of NGC 5548 by season

Per season Synergistic-Unique-Redundant Decomposition (SURD) of the 5100 A
continuum and the total, blue, core and red H beta fluxes of NGC 5548, for the five
historical AGN Watch seasons (1989 to 1993), following the season by season strategy
of Lu et al. 2022 and Xi et al. 2025. The paper is `main.tex` (MNRAS class).

| Location | Content |
| :--- | :--- |
| `main.tex`, `main.pdf`, `refs.bib` | The paper. `main.tex` inputs the results section and the figure atlases below |
| `collaborator_rebuild/` | Analysis code, tests and notes (start with `EPOCH_ANALYSIS.md`) |
| `collaborator_rebuild/results_section/` | The LaTeX results section, its tables, and a standalone `preview.tex` |
| `collaborator_rebuild/epoch_results/` | Scan settings, hashes, every component of every curve, extrema, and the figure atlases `figures/atlas_*.pdf` |
| `collaborator_rebuild/calibration_review/` | Categorical benchmark of the estimator (Table 8 of the paper) |
| `pipeline/` | Retrieval and reconciliation of the AGN Watch inputs |
| `audit/iw_bp_96_tab.txt` | Corrected table of Wanders and Peterson 1996 (line fluxes) |
| `output/` | Replies to the collaborator's instructions (`docx`) and the standalone results section (`pdf`) |

The earlier validation study (pooled window, ICCF/JAVELIN/kNN cross-checks) is the git
tag `validation-study-main-tex`. The full tree before the clean up is the tag
`pre-cleanup-oct14`.

## Reproduce

Use Python 3.11 and the official SURD source (commit
`79dbdea85e6754ec2b5457b3e37204c5d53d1815`, cloned to `SURD/`, which Git ignores).

```sh
python3.11 -m venv .venv && source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python -m pipeline.fetch_inputs                       # AGN Watch inputs, hashes verified
python -m collaborator_rebuild.prepare_historical     # five series per season
python -m collaborator_rebuild.epoch_scan             # 7200 decompositions and figures
python -m collaborator_rebuild.summarize_epoch        # tables for the paper
python -m pytest -q                                   # tests
tectonic main.tex                                     # or compile main.tex on Overleaf
```

Compiling needs the MNRAS class `mnras.cls` and `mnras.bst` (the Overleaf MNRAS template).
