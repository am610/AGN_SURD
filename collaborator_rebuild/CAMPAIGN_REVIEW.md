# Campaign support and component sensitivity

The published flux scans were reviewed on 7 October 2026. All 25 native input
files agree byte for byte between the two bin settings. Saved settings agree
apart from bin count, saved runner hashes agree, and all 3720 configurations
have identical tuple counts and support status. Configuration hashes agree
with their run manifests. The 1326 supported positive lag configurations each
contain the same 26 components with finite fractions summing to one.

The comparison holds native support fixed between bin settings at each lag.
It does not hold the target dates fixed across different lags. Consequently
lag maxima can also reflect changing sampling support. No new SURD scan was
needed for this review, and no probability estimation setting was selected.

| Campaign | Supported positive lag configurations out of 720 | Median total variation | Median absolute leakage change | Component curves with disjoint peak sets |
| :--- | :--- | :--- | :--- | :--- |
| 1989 | 111 | 0.465970 | 0.252885 | 35.38 percent of 260 |
| 1990 | 68 | 0.573625 | 0.109833 | 67.95 percent of 78 |
| 1991 | 1 | 0.473586 | 0.163987 | Not assessable |
| 1992 | 426 | 0.403103 | 0.130625 | 68.11 percent of 624 |
| 1993 | 720 | 0.515113 | 0.209931 | 87.02 percent of 624 |

Total variation is half the absolute difference summed over all normalized
components at a matched configuration. The leakage comparison uses its own
target entropy normalization. Peak sets retain all ties within an absolute
tolerance of 1e-12. The last column includes only individual component curves
with at least two supported positive lags. All components are included without
a minimum amplitude criterion, so the peak diagnostic also counts weak
components. It is a descriptive count, not a detection statistic or an
estimate of delay uncertainty.

The 1993 campaign passes the 32 tuple numerical guard for every positive lag
configuration. Even there, every total variation exceeds 0.2. Its median
tuples per occupied histogram cell are 3.071429 for two bins and 1.811254 for
three bins. These occupancy ratios summarize counts, not independent
observations or effective sample sizes. Interpolation and time dependence
remain relevant. More complete numerical support has not resolved the
sensitivity to bin choice.

The 1991 campaign has only one supported positive lag configuration. It cannot
supply a component lag curve under the retained guard. The 1989 and 1990
campaigns support mostly the continuum target. The 1992 and 1993 campaigns
are the practical starting points for further estimation checks, while their
present curves remain strongly sensitive to bin choice.

The current evidence supports reporting sampling limitations and estimator
sensitivity. It does not establish a stable positive astrophysical result,
select either bin setting as reliable, or justify expanding the predictor
vector. Before a physical interpretation, the next review should compare
lags on an intersection of target dates within each curve and assess whether
that stricter support remains adequate. Estimator calibration would still be
needed even if that comparison were stable.

That subsequent comparison is now complete in `FIXED_SUPPORT_REVIEW.md`.
Only four continuum target curves retain at least 32 dates across all
positive lags, and substantial sensitivity to bin choice persists.

Reproduce from the repository root:

```sh
python -m collaborator_rebuild.review_campaigns
python -m pytest -q collaborator_rebuild/test_campaign_review.py collaborator_rebuild/test_lag_design.py collaborator_rebuild/test_total_flux.py
```

The saved scans must already exist. Derived tables, the sensitivity figure,
and the review manifest are under `published_results/campaign_review`.
The manifest records the reviewer source hash and hashes of the reviewed
inputs. Original scans remain intact. Eight focused tests passed, including
rejection of missing atoms, invalid normalization, and preservation of tied
peak lags. The figure was rendered and visually inspected.
