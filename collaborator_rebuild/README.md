# Active collaborator study

The user selected historical NGC 5548 observations for the first analysis.
Modern Lu and Xi observations remain an optional second layer, not a
prerequisite for the authorized historical study.

The scope is five observed series by campaign, both specified four variable
subsets, every future target including its own history, individual normalized
SURD components and leakage, common and unequal lags, thesis equation 3.18,
and a LaTeX results section. Use the unchanged authors' SURD decomposition.

Read [progress](PROGRESS.md), [specification](SPECIFICATION.md),
[review](REVIEW.md), and [data search](DATA_SEARCH.md).

The historical runner is exploratory. Total flux provenance and estimation
reliability remain open. Equation 3.18 construction passes tests, but only its
support audit has run. No expanded scientific result or causal detection is
claimed.

Run from the repository root:

```sh
python -m pytest -q collaborator_rebuild/test_lag_design.py
python -m collaborator_rebuild.check_support
python -m collaborator_rebuild.run_historical
python -m collaborator_rebuild.plot_reviewed
```

Downloads and derived exploratory outputs are local and ignored by Git. The
support checker downloads the historical broad line candidate and saves its
provenance. Optional modern preparation requires the two tables described in
DATA_SEARCH.md in sources/ before running:

```sh
python -m collaborator_rebuild.prepare_public_series
```

Existing historical outputs retain their original configuration. Rerunning
the runner regenerates them. Resolve scientific settings before promoting
results to tracked release evidence or writing definitive conclusions.
