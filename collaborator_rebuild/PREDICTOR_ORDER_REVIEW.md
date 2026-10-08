# Predictor ordering review

The bounded audit preserves histogram counts, target dates, lag range, bin
settings, and predictor identities for the four retained astronomical curves.
It checks all 24 predictor permutations across 240 saved configurations,
for 5760 evaluations. Labels are permuted with axes and mapped back before
comparing individual components. The numerical comparison tolerance is 1e-8;
it is not a statistical significance threshold.

Eleven configurations exceed that tolerance. Maximum total variation is
0.1878142285. Five affected configurations use 1992 common with two bins,
three use 1992 common with three bins, two use 1993 blue earlier with two
bins, and one uses 1993 red earlier with three bins. Total information changes
by at most 6.67e-16 bits and leakage by at most 1.84e-15.

The 1993 common curves are stable under predictor ordering. Their substantial
bin sensitivity therefore cannot be explained solely by ordering. No preferred
permutation is selected. The unchanged authors' source is retained and its
hashes are recorded. This audit establishes computational sensitivity, not a
physical effect or a correction to the algorithm.

Compact evidence, configurations, witnesses, and hashes are in `order_review`.
The full order table remains local in `published_results/predictor_order_audit`.
The revised section includes the curve summary as Table 5.

Reproduction from the repository root, with the saved scans available:

```sh
python -m collaborator_rebuild.audit_predictor_order
python -m pytest -q collaborator_rebuild/test_predictor_order.py
```

This completes the agreed audit. Analysis is frozen as specified in
`REVIEW_HANDOFF.md`.
