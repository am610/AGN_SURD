# Historical run review

Reviewed on 6 October 2026. This is a diagnostic review, not an approved final
scientific analysis. No scientific settings or released outputs were changed.

## Verdict

The official SURD decomposition and numerical normalization pass. The current
run remains an exploratory implementation and does not complete the feedback.

## Confirmed checks

The implementation uses the unchanged authors' SURD routine with source hash
checks. Each target includes all four predictors, including its own history.
Both requested subsets are present. All 1402 usable configurations contain
26 components. Their normalized fractions sum to one with maximum absolute
error 1.33e-15. This checks arithmetic, not statistical reliability.

Published historical intervals in Peterson et al. 2004 Table 2 are JD offsets
47509 through 47809, 47861 through 48179, 48225 through 48534, 48623 through
48898, and 48954 through 49255. All adopted line and continuum observations
fall inside these intervals, and their existing epoch labels agree exactly.
Thus the November convention has not misassigned these observations. Explicit
published intervals should nevertheless replace the convention for provenance
and future input handling.

Source: https://arxiv.org/abs/astro-ph/0407299

## Issues requiring attention

1. Total currently means the sum of blue, core, and red within the adopted
   velocity limits. It is not an independently measured full broad H beta
   light curve. This creates an exact algebraic relationship in the total
   subset. The AGN Watch archive provides a revised broad H beta series with
   narrow emission removed and the published factor of two error corrected.
   Its calibration, date matching, and integration definition must be checked
   against the adopted profiles before selecting it as the total input.

2. Positive lag coverage is incomplete. In the common scenario, line targets
   in 1989 support only zero through three days; those in 1990 and 1991 support
   only zero days under the current numerical guard. Continuum targets have
   different sampling and support. There are no usable unequal scenarios in
   1991. Results cannot be presented as complete lag curves for every epoch.

3. A minimum of 32 tuples does not validate a histogram with up to 243 cells.
   Many configurations are sparsely populated. Bin and sample support review
   is required before interpreting components. Structural dependence can also
   produce empty cells, so occupancy alone is not a validity test.

4. Exact contemporaneous target history makes zero lag leakage essentially
   zero by construction. This is an algebraic control, not a physical finding.

5. The blue and red five day offsets are stipulated exploratory settings.
   They are not measured delays or a complete implementation of thesis
   equation 3.18, which requires multiple shifted copies of predictors.

6. Plotting needs correction. Single point configurations appear empty because
   there are no markers. Colors repeat across component labels. Unsupported
   lag ranges and missing configurations are not clearly displayed. There
   are 76 plots out of 120 possible epoch, subset, target, scenario cases.
   All 76 were inspected in overview contact sheets, with a representative
   full size plot inspected separately; this is not full resolution review
   of every plot.

7. The requested source grounded LaTeX results section is not yet written.

Archive source:
https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/
Candidate revised broad line data:
https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/iw-hb.lcv

## Next implementation order

Resolve total flux provenance and encode published campaign intervals. Review
usable positive lag coverage and bin support without silently relaxing guards.
Implement the explicit expanded predictor construction from equation 3.18,
with its own support accounting. Repair figures to show observed points,
unsupported regions, distinguishable components, and tuple counts. Only then
write the scoped LaTeX results section with conservative interpretation.

The historical data can remain the first study. Modern campaign observations
are a separate optional extension, not a prerequisite established by this
review. Collaborator acceptance cannot be guaranteed from numerical checks.
