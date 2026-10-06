# Exact analysis specification

The September feedback governs this rebuild. Preparation of observed time
series precedes the scientific decomposition and interpretation.

## Variables and targets

Prepare F5100, integrated broad H beta, blue H beta, core H beta, and red H beta
with native dates, flux uncertainties, units, and provenance for every selected
epoch. Epoch membership must refer to an observing campaign definition.

Analyze the subset total, blue, core, red, and the subset F5100, blue, core, red.
Each of the four variables becomes the future target in turn. All four past
variables remain predictors, including the target's own history.

For four predictors there are four unique components, eleven redundancy
components, and eleven synergy components. Thus the explicitly listed
combinations produce 26 components. The feedback calls them 25, but removing
a component to match that count would contradict its listed combinations.

## Normalization and figures

Equation 3.16 averages specific contributions over future target states.
Equation 3.17 defines the expected order used to arrange plotted components;
it does not define an additional normalization.

Following thesis Figure 4.1, divide every information component by the joint
mutual information between the future target and the complete predictor vector.
These 26 fractions sum to one when joint information is nonzero. Divide the
conditional entropy in equation 3.6 by future target entropy for normalized
leakage. Leakage has a separate denominator and separate plot.

Plot individual unique, redundant, and synergistic components against lag.
Separate panels may keep the curves readable while preserving every component.
Do not substitute total unique or total synergy curves for the individual atoms.

## Multiple lags

Equation 3.18 expands the predictor vector to include additional time shifted
copies. Assigning only one distinct delay to each variable is a different
predictor construction. The additional copies and their delays must be specified
before executing this part. Atom counts change when the predictor vector grows.

## Parameters awaiting specification

The September feedback leaves placeholders for common lag ranges and predictor
specific lags. The local task appendix supplies recommended common ranges
for six seasons and candidate velocity bins; see `DATA_SEARCH.md`. Exact
predictor specific delays remain absent. Time frame,
sampling alignment, and probability estimation settings remain unspecified. Save
these settings before any scientific scan. The former aggregate configuration
does not supply collaborator approval for these choices.

## Current implementation

`official_decomposition.py` calls the unchanged authors' code and verifies
source hashes before importing it. It exposes every component and the two
normalizations. Analytic verification is a software check; it is not an AGN
result. The authorized historical study has an exploratory run, with published
campaign intervals now encoded in `lag_design.py`. Equation 3.18 construction
is implemented and tested, but only its support audit has run. Total flux
provenance and estimation reliability remain open; see `PROGRESS.md`.
