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

Corrected published line columns are now available through
`prepare_historical.py`. See `BINNING_REVIEW.md` for the new runs and substantial
bin sensitivity. Use the explicit `--published` runner option for these inputs;
the default runner remains the earlier profile based exploratory path.

The saved scans now have a campaign review in `CAMPAIGN_REVIEW.md`, produced
by `python -m collaborator_rebuild.review_campaigns`. Even the fully supported
1993 campaign remains strongly sensitive to bin choice. The review includes
individual component changes, tied peak sets, occupancy diagnostics, and a
figure marking configurations that fail the numerical support guard.

## Current review boundary

The bounded diagnostic work is complete and the analysis is frozen for
collaborator review. Initial categorical calibration and the predictor ordering
audit are complete. The separate results section now has five tables, four
figures, and a ten page PDF. All pages were visually inspected.

Ordering changes component allocation in 11 of 240 configurations. The maximum
total variation is 0.187814. Total information and leakage remain unchanged to
numerical precision. The 1993 common curves are stable under ordering but remain
strongly sensitive to bin choice. No stable physical delay is established.

Read [review handoff](REVIEW_HANDOFF.md) for completed requirements, remaining
gaps, and the stopping rule. Further scientific runs or changes to the authors'
routine require a new agreed scope. Earlier progress entries below are history.

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
