# Per season SURD analysis (14 October 2026)

This is the analysis that answers the September instructions. It replaces the
native date scans for the main results; those remain in
the earlier review notes, which are in the git tag `pre-cleanup-oct14`.

## Design

| Item | Setting |
| :--- | :--- |
| Epochs | Five historical seasons, JD less 2400000: 47509 to 47809, 47861 to 48179, 48225 to 48534, 48623 to 48898, 48954 to 49255 (Lu et al. 2022, Table 6, Years 1 to 5) |
| Series | F5100, total, blue, core and red H beta (corrected published columns, Wanders and Peterson 1996) |
| Variable sets | {total, blue, core, red} and {F5100, blue, core, red} |
| Targets | Every member of each set in turn; all four members are past predictors |
| Grid | Daily integer grid per epoch; linear interpolation; no extrapolation; no bridging of gaps above 15 days |
| Lags | Base lag 1 to 30 days in one day steps |
| Bins | Two (primary) and three (sensitivity) marginal quantile bins from the full native epoch series |
| Common lags | Every predictor at X(t minus L) |
| Uncommon lags (eq. 3.18) | Each predictor variable has its own lag L plus delta. Blue leading: blue delta 5 days. Red leading: red delta 5 days |
| Estimator | Unchanged authors' SURD routine, source hashes verified |
| Normalization | Atoms divided by joint mutual information (26 atoms, sum 1); leakage divided by target entropy |

7200 decompositions (5 epochs, 2 sets, 4 targets, 3 lag configurations, 30 lags,
2 bins), all with defined normalization, 187200 component rows. Every lag keeps
at least 118 grid tuples. Grid tuples are interpolated and autocorrelated.

The instructions leave the lag values blank ("AYAN"). The 1 to 30 day range
covers the seasonal common lag ranges proposed in `DATA_SEARCH.md`. The five day
offset is a stipulated analysis choice, about the size of published H beta lag
uncertainties; it is not a measured delay.

Equation 3.18 is applied with one lag per predictor variable so that the four
predictors still give the 26 requested components. The literal expansion with a
second copy of every variable would have eight predictors and 502 components.
It is implemented in `lag_design.expanded_predictors` but is not used here.

## Reproduce

From the repository root, after `python -m collaborator_rebuild.prepare_historical`
(writes `collaborator_rebuild/historical_prepared`):

```sh
python -m collaborator_rebuild.epoch_scan        # scans and figures, about 90 s
python -m collaborator_rebuild.summarize_epoch   # tables and extrema
mkdir -p tmp/pdfs/surd_section_build
tectonic --keep-logs --outdir tmp/pdfs/surd_section_build collaborator_rebuild/results_section/preview.tex
cp tmp/pdfs/surd_section_build/preview.pdf output/pdf/surd_epoch_results.pdf
python -m pytest -q collaborator_rebuild
```

## Outputs (`epoch_results/`)

`configuration.json` and `manifest.json` (settings, input and source hashes),
`support_and_leakage.csv` (tuple counts and leakage for every scan),
`individual_components.csv.gz` (all 26 components), `descriptive_extrema.csv`,
`bin_sensitivity.csv`, `mean_leakage.csv`, `config_leakage_comparison.csv`,
`continuum_peak_lags.csv`, and `figures/atlas_{common,blue_leading,red_leading}.pdf`
(40 pages each, one page per epoch, set and target).

## Main findings

- Leakage is lowest at the shortest lag (36 of 40 curves, two bins) and highest at
  20 days or more (35 of 40): the target's own history dominates at short lags.
- Epoch mean leakage ranks like the grid tuple count, as the categorical
  calibration predicts, so it is not a physical comparison between epochs.
- In the F5100 set the unique continuum component of the line targets peaks at
  16 to 29 days in 13 of 15 curves. Median peak lags 23, 27, 20.5, 17, 23 days
  compare with the published 19.7, 18.6, 15.9, 11.0, 13.0 days (Lu 2022).
- The red peak comes before the blue peak in 1989, 1991, 1992 and 1993, and the
  order is reversed in 1990: red-leads-blue in four of five seasons, the
  signature Xi et al. 2025 report for 2019 and 2021.
- Moving one wing five days back changes mean leakage by at most 0.08 but
  changes the component allocation by a median total variation of 0.30.
- Two and three bins differ by a median total variation of 0.27 to 0.58.

No significance level or flux uncertainty propagation has been computed.
