# Implementation progress

Current active preparation now uses all four corrected published line columns.
`prepare_historical.py` verifies 242 instrument and date records and exports
five series per campaign. Two and three bin descriptive scans are complete
in a separate local output area. Their strong sensitivity is documented in
`BINNING_REVIEW.md`; no setting is promoted as a reliable physical result.

The campaign review is complete in `CAMPAIGN_REVIEW.md`. It verifies identical
native inputs and tuple support between the saved bin settings. All 720
positive lag configurations in 1993 pass the numerical guard, but median
component total variation is 0.515113 and 87.02 percent of individual component
curves have disjoint peak sets. This remains descriptive sensitivity, not a
physical delay result. Eight focused tests pass. Support is matched between
bin settings at each lag, but still changes across lags within a curve.

The subsequent fixed support scan is complete in `FIXED_SUPPORT_REVIEW.md`.
Intersecting target dates across the declared 1 through 30 day lag range
leaves four of 120 curves above the unchanged 32 tuple guard. All four have
future F5100 as the target. Both bin settings were run, yielding 240 usable
decompositions and 6240 component rows. Median component total variation
ranges from 0.522899 through 0.609769. Twelve focused tests pass, and all four
figures were visually checked. No stable physical delay is established.

The separate LaTeX results section is now drafted under `results_section`.
It contains three tables and four individual component figures, including
descriptive extrema and a limited contextual comparison with Lu and Xi.
The nine page PDF preview compiles without layout or unresolved reference
warnings, and all pages have been visually checked. The section reports
completed diagnostics and explicitly identifies unfinished scientific work.
It has not been incorporated into the previous full manuscript.

The first estimator calibration is complete in `ESTIMATOR_CALIBRATION.md`.
Eight exact categorical laws and 6400 sampled decompositions show substantial
finite sample information bias under independence. Four saved witnesses
also show component allocation changes when predictor axes and their labels
are reordered together, while total information stays unchanged to numerical
precision. Tracing the unchanged authors' routine records resets after
differences around machine precision. All 22 active tests pass. An axis order
audit of the four retained astronomical curves is now the immediate next
step. Observation matched calibration and a revision of the descriptive
results section remain pending.

The additional total integration audit is implemented in `audit_total_flux.py`.
See `TOTAL_FLUX_PROVENANCE.md`. The source section interval agrees much more
closely with archive broad line daily means than the current velocity sum.
No fitted calibration or scientific input replacement has been applied.

Published campaign intervals now replace the calendar convention in the runner.
Existing observation labels are unchanged. Existing result files remain intact;
their saved configuration still documents the original run.

The revised archive broad H beta source is saved with a checksum. All 224
adopted spectral dates have a match. The source contains 312 measurements and
cannot be uniquely matched to individual exposures from dates alone. Daily
means correlate with the profile sum at 0.998730, but flux ratios vary from
0.764439 to 0.974908. Correlation does not establish identical integration or
calibration. The source remains a candidate, not an automatic replacement.

The equation 3.18 construction is implemented and tested. It keeps the base
and an additional shifted copy of each predictor. A five day exploratory
spacing is used only for the support audit, not as a claimed physical delay.
Three tests pass for campaign membership, shifted values, and no extrapolation.

Positive lag configurations passing the existing 32 tuple numerical guard:

| Construction | Bins | Configurations | Tuple range | Possible cells | Median occupied cells | Median singleton cells |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Four predictors | 2 | 562 | 32 to 164 | 32 | 12 | 5 |
| Four predictors | 3 | 562 | 32 to 164 | 243 | 22 | 12 |
| Eight predictors | 2 | 382 | 32 to 155 | 512 | 24 | 15 |
| Eight predictors | 3 | 382 | 32 to 155 | 19683 | 32 | 26 |

These counts do not validate an estimator. Sparse support, particularly for
the expanded vector, prevents treating a routine execution as scientific
completion. No expanded SURD decomposition is reported as an AGN finding.

The reviewed figure renderer preserves the original descriptive atom estimates
but adds markers, distinguishable component styles, unsupported lag shading,
tuple counts, and explicit panels for cases with no supported lags. It renders
all 120 requested combinations into a separate directory.

The published total column supplies the current total flux input; exact
exposure and integration reconstruction from the profile archive remains
unresolved. Support held fixed across lags within each curve has now been
evaluated. Remaining work includes calibration and verification of probability
estimation settings with adequate support, expanded scientific runs where
justified, and scientific review and revision of the drafted results section.
The additional historical light curve is publicly available; using modern
campaign data is not required to perform these checks.
