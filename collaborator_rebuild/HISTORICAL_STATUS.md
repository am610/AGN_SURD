# Historical implementation status

The user authorized the available historical observations as the first study.
The Lu and Xi observations are a possible later extension. This resolves the
data choice for implementation without asserting collaborator approval.

Run `python -m collaborator_rebuild.run_historical` from the repository root.
The runner prepares all five native series for five project epochs, named
1989 through 1993, now defined by the published intervals in Peterson et al.
2004 Table 2. The initial saved outputs used a November through October
convention, which assigns every current observation to the same campaign.
No lag pairing crosses campaign boundaries.

The total is the sum of the adopted broad profile windows. It differs from
the separately archived integrated H beta series. Adopted velocity windows
retain the source data preparation convention at boundaries of minus 6000,
minus 2000, plus 2000, and plus 6000 kilometres per second. Recalculation with
the narrower candidate windows remains a possible sensitivity check.

Every target uses all four past predictors, including itself. The unchanged
authors' SURD code calculates 26 individual components. Three bins use fixed
native marginal quantile boundaries within each epoch. Targets remain on
their native dates; predictors interpolate across gaps no greater than ten
days. A minimum of 32 native target tuples is a numerical guard. It does not
establish reliable estimation of a histogram with up to 243 cells.

The base lag scan spans zero through thirty observed days. Zero lag includes
the identical target as a predictor and therefore serves as an algebraic
control. Zero lag results must not be interpreted as temporal causality.

Two exploratory unequal delay configurations shift blue or red five days
earlier while leaving the other predictors at the base lag. These are
stipulated sensitivity settings, not measured physical delays. They do not
implement the expanded predictor vector in thesis equation 3.18.

All usable rows have 26 components whose fractions sum to one within numerical
precision. The 1991 epoch lacks enough support for either unequal setting.
Figures show every individual component in separate unique, redundancy, and
synergy panels, followed by a leakage panel. They are descriptive estimates.

Equation 3.18 construction is now implemented and tested, and its sampling
support has been audited without claiming a validated expanded decomposition.
Reviewed figures show all cases and unsupported lags. Remaining work includes
total flux provenance, estimation reliability, justified expanded runs, and a
LaTeX results section with source checked astrophysical comparison. The current
implementation is an initial scoped run, not completion of all feedback.
