# Implementation progress

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

Remaining work is the exposure and integration provenance of total H beta,
selection and verification of probability estimation settings with adequate
support, expanded scientific runs where justified, and the results section.
The additional historical light curve is publicly available; using modern
campaign data is not required to perform these checks.
