# Lag curves with identical target dates

The comparison was completed on 7 October 2026 using the published flux
inputs and the unchanged authors' SURD code. Within each campaign, subset,
target, and scenario, target dates are retained only when every predictor
is available at every declared positive base lag from 1 through 30 days.
The same target dates are then used for both bin settings and every lag.
There is no shortened lag range or reduced tuple guard.

The ten day interpolation gap limit, stipulated scenario offsets, and native
campaign quantile boundaries are retained from the earlier scans. Bin edges
are not refitted at different lags or to the retained target subset. Zero lag
is excluded from this positive lag comparison. Configuration was saved before
execution. Native inputs and the authors' source hashes were verified.

Only four of 120 curves retain at least 32 tuples across the full lag range.
The other 116 curves are recorded as insufficient support and have no
decomposition in this run. All surviving curves have future F5100 as the
target and past F5100, blue, core, and red as predictors. Their target direction
is continuum prediction from its own history and past line components.

| Campaign and scenario | Fixed tuples | Median total variation between bin settings | Median absolute leakage change | Components with disjoint peak sets |
| :--- | :--- | :--- | :--- | :--- |
| 1992 common | 36 | 0.531073 | 0.124605 | 22 of 26 |
| 1993 common | 90 | 0.594429 | 0.122200 | 21 of 26 |
| 1993 blue earlier | 78 | 0.609769 | 0.143941 | 25 of 26 |
| 1993 red earlier | 78 | 0.522899 | 0.169819 | 21 of 26 |

The blue earlier and red earlier scenarios retain the stipulated additional
five day predictor shifts. These are analysis choices, not measured delays
or the expanded predictor construction of equation 3.18. Support is fixed
within each scenario, but differs between scenarios. Scenario differences
therefore cannot be attributed solely to a predictor delay change.

The campaign ranges of retained tuple counts are 1 through 7 in 1989,
1 through 7 in 1990, zero in 1991, 22 through 36 in 1992, and 25 through 90
in 1993. Every future line target falls below the guard. Passing the guard
is a numerical condition and does not establish estimator reliability.

The sensitivity remains large after removing changes in target dates across
lags. At a fixed bin count, replacing the earlier varying support with the
fixed support also changes the components substantially. For the 1993 common
curve, the median total variation from this support change is 0.474836 with
two bins and 0.390653 with three bins. Across all four curves and both bin
settings, these medians range from 0.272211 through 0.474836. This comparison
measures sensitivity to the retained observations, including possible temporal
variation; it does not isolate a statistical bias or prove nonstationarity.

Total variation uses the 26 normalized information components, while leakage
uses its separate target entropy denominator. Peak sets retain ties within
an absolute tolerance of 1e-12. All components are counted, including weak
components. Disjoint peak sets are descriptive diagnostics and do not measure
statistical significance or uncertainty in a physical delay.

There are 240 usable decompositions, from four curves, 30 lags, and two bin
settings, producing 6240 component rows. All have defined normalization.
Maximum absolute deviation of the component sum from one is 1.22e-15.
That agreement verifies arithmetic only. The 1993 common curve has median
occupied cell counts of 14.5 and 25 for two and three bins, respectively,
with 90 retained tuples. Temporal dependence, interpolation, and flux errors
have not been calibrated by these counts. No null test or uncertainty
propagation is reported.

Neither bin setting is selected as reliable. This review does not establish
a stable physical delay or causal detection and does not justify expanding
the predictor vector. Probability estimation calibration remains the next
methodological task. A calibration must check recovery against a known
probability law at comparable support and separately address temporal
dependence and interpolation before it can support physical interpretation.

Reproduce after the saved published scans are available:

```sh
python -m collaborator_rebuild.fixed_support
python -m pytest -q collaborator_rebuild/test_fixed_support.py collaborator_rebuild/test_campaign_review.py collaborator_rebuild/test_lag_design.py collaborator_rebuild/test_total_flux.py
```

The outputs are local under `published_results/fixed_support`. They include
the frozen configuration, curve support inventory, retained target dates,
support and leakage table, all individual components, bin comparisons,
component peak sets, comparisons with the previous varying support estimates,
four figures, and a provenance manifest. Earlier scans remain intact.

Twelve focused tests passed. Checks cover date intersections, interpolation
without extrapolation, unequal offsets, application of the tuple guard to
the intersection, degenerate targets, missing atoms, normalization, and tied
peaks. Output verification confirms complete 30 lag grids and identical date
hashes and tuple counts within every executed curve. All four figures were
rendered and visually inspected.
