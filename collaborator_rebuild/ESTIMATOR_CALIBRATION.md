# First stage estimator calibration

Completed on 7 October 2026 using the unchanged authors' SURD routine.
The experiment uses four predictors and one target, two or three known
categorical states, 200 independent repetitions per setting, and a fixed
random seed of 20261007. Sample sizes of 36, 78, and 90 match the retained
historical curves. An additional 512 tuple sample is an artificial reference.
There are eight exact probability laws and 6400 sampled decompositions.
All sampled targets have defined normalization in this run.

## Known probability laws

All variables are uniform in their allowed states. The predictor labels
are continuum history, blue, core, and red, with future continuum as target.
The labels are bookkeeping for artificial variables, not a physical AGN model.

| Model | Construction | Exact information structure |
| :--- | :--- | :--- |
| Independent | Target and every predictor are independent | Zero joint information; leakage one |
| Unique | Target equals blue; predictors are mutually independent | All information in U(blue) |
| Redundant | Target equals blue equals core; other predictors independent | All information in R(blue,core) |
| Synergistic | Target equals blue plus core modulo the number of states | All information in S(blue,core) |

The nonindependent models have exact joint information equal to the base
two logarithm of the state count and zero leakage. The exact laws pass
analytic checks through the official routine in the declared predictor order.
In the independent exact law, normalized components are undefined because
joint information is zero. They are not assigned artificial zero fractions.

## Finite sample bias under independence

The estimated information should be zero in this benchmark. Sampled
histograms nevertheless produce substantial positive information, particularly
with three states. Estimated leakage is also substantially below its true
value of one.

| Independent tuples | Median information with two states, bits | Median information with three states, bits | Median leakage with two states | Median leakage with three states |
| :--- | :--- | :--- | :--- | :--- |
| 36 | 0.352703 | 1.303841 | 0.640860 | 0.157221 |
| 78 | 0.158212 | 1.075949 | 0.839102 | 0.312568 |
| 90 | 0.136304 | 1.015025 | 0.862211 | 0.353165 |
| 512 | 0.019037 | 0.268462 | 0.980927 | 0.830381 |

These values come from the independent categorical law only. They are not
thresholds or p values for the observed AGN data, which have different
temporal and predictor dependence. The 32 tuple guard cannot establish
estimator reliability by itself. More samples reduce this benchmark bias,
but the 512 tuple reference is not an observational sample or a selected
minimum sample size.

## Recovery of perfect signals

At 90 tuples, the median fraction allocated to the known component is
0.982385 and 0.964760 for the unique model with two and three states,
0.991510 and 0.970957 for the redundant model, and 0.955916 and 0.861969
for the synergistic model. These are exceptionally strong deterministic
signals. Their recovery does not establish sensitivity to weak physical
relationships, nor does a favorable median describe every repetition.

There are substantial allocation failures in some binary controls. For
example, the fifth percentile of the expected unique component fraction
at 90 tuples is 0.505330. In the artificial 512 tuple binary unique case,
the expected component is the sole largest in only 77 percent of repetitions,
despite its median fraction being 0.996266. This nonmonotonic behavior
requires a numerical audit rather than an interpretation based on sample
size alone.

## Reproducing numerical allocation sensitivity

Four saved witnesses record original integer histogram counts, their
decomposition, and all 24 reorderings of predictor axes. Each reordering
moves the predictor labels with its axes. Thus it represents the same
distribution and the same named variables. Component labels are compared
after sorting the names within each component.

| Witness | Repetition | Maximum component total variation under reordering | Maximum absolute joint information change, bits |
| :--- | :--- | :--- | :--- |
| Binary unique, 90 tuples | 4 | 0.989753 | 2.22e-16 |
| Binary unique, 512 tuples | 1 | 0.993109 | 1.11e-16 |
| Binary synergy, 90 tuples | 41 | 0.511119 | 2.22e-16 |
| Binary synergy, 512 tuples | 16 | 0.506673 | 4.44e-16 |

An execution trace observes the unchanged routine at its strict comparison
and reset statement, `I1[inds_] = 0`. In these witnesses, positive specific
information values are reset after differences of 1.11e-16 through
4.44e-16 from the comparison maximum. These observed resets and the
axis order sensitivity identify a numerical implementation concern in
these artificial examples. The total component sum still normalizes.
Normalization therefore cannot validate the component allocation.

The authors' source files and the active wrapper have not been modified.
No altered numerical routine is promoted as a replacement. These examples
do not establish that this issue explains the observed AGN bin sensitivity
or affects every input. Its prevalence and effect on the saved astronomical
histograms remain to be assessed.

## Scope and next work

This completes an initial calibration with independent categorical draws.
It does not calibrate temporal dependence, irregular observing dates,
interpolation, quantile boundary estimation, measurement errors, lag scans,
or weak physical signals. The two and three state laws are separate
benchmarks, not different binning of the same continuous observations.
No physical delay recovery or AGN significance test is claimed.

The subsequent astronomical ordering audit is complete in
`PREDICTOR_ORDER_REVIEW.md`. The revised results section incorporates both
diagnostics. Analysis is now frozen for collaborator review; realistic
observation calibration remains unperformed. See `REVIEW_HANDOFF.md`.

Reproduce from the repository root:

```sh
python -m collaborator_rebuild.calibrate_estimator
python -m pytest -q collaborator_rebuild/test_estimator_calibration.py
```

Full repetition tables remain local under
`published_results/estimator_calibration`. Compact configuration, summary,
exact laws, witness histograms, figure, and hashes are saved under
`calibration_review`. Ten new analytic tests passed, and all 22 active
collaborator tests pass. The calibration figure was visually inspected.
The revised ten page PDF includes these diagnostics and the astronomical
ordering audit.
