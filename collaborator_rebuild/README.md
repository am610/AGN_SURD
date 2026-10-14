# Analysis code

Start with [EPOCH_ANALYSIS.md](EPOCH_ANALYSIS.md) (design, outputs, findings),
[SPECIFICATION.md](SPECIFICATION.md) (how the September instructions are implemented)
and [DATA_SEARCH.md](DATA_SEARCH.md) (what is and is not publicly available from
Lu et al. 2022 and Xi et al. 2025).

| File | Role |
| :--- | :--- |
| `official_decomposition.py` | Calls the authors' unchanged SURD routine (source hashes verified) and returns all 26 components and normalized leakage |
| `prepare_historical.py` | Builds the five series per season from the corrected published table and the AGN Watch continuum |
| `lag_design.py` | Season intervals and the equation 3.18 predictor construction |
| `epoch_scan.py`, `plot_epoch.py`, `summarize_epoch.py` | The per season scans, figures and tables |
| `calibrate_estimator.py` | Categorical benchmark with known probability laws |
| `build_replies_docx.py` | Inserts the coloured replies into the collaborator's instructions |
| `test_*.py` | Tests |
