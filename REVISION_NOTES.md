# Changes warranted by the adopted pipeline

The local working draft `final_draft.tex` still describes earlier calculations and contains definitions that disagree with the verified loader. The historical manuscript sources already in this repository are not a submission release. The new pilot does not provide calibrated significance for the adopted sample.

## Update now

1. Replace MJD labels with JD less 2400000 where the plotted coordinates are 47512 to 49255. The filename rule produces that JD offset, not MJD. The revised continuum file instead needs 40000 added to its dates, as checked against corrected Table 1.
2. Correct the adopted velocity boundaries to minus 6000, minus 2000, plus 2000 and plus 6000 km/s. Document observed wavelengths, the project reference redshift 0.017175 and reference wavelength 4861.33 Å. These are explicit project choices and do not reproduce the exact bins or redshift used by Wanders and Peterson.
3. State that 247 archived spectra occur on 226 dates, while the adopted sample uses 242 spectra on 224 dates after excluding five profiles without exact date and instrument matches in the corrected table. The original continuum contributes 557 observations. List the five exclusions and preserve their provenance limitation.
4. Explain integration over 2 Å pixels, equal averaging on repeated dates, propagated errors, and missing calibration and pixel covariance. The integrated component sum is not the archive total Hβ series with its different subtraction conventions.
5. Replace statements treating 1744 interpolated values as independent samples. Report observed dates and lagged tuples separately. Native support endpoints differ between line and continuum inputs, and the sampler enforces each support at the requested lag without extrapolation.
6. Describe the four sampling methods and new exploratory settings. Keep the new results distinct from the earlier eight bin scans: the current pilot uses four marginal quantile bins and only the three line targets. It excludes each target component from its predictors at every lag.
7. Correct the history conditioning caption: 0.4043 bits is described by the underlying historical analysis as the median surrogate maximum, not a 95 percent global threshold. Nonsignificance under a particular null does not establish that the data are random noise or identify a unique bias mechanism.

The accompanying `reconciliation_update.tex` is an insertable methods and status section with tables generated from the adopted output files. It is a review fragment, not a rewritten complete manuscript. Integrate it while replacing the conflicting legacy methods text. Merely appending it would leave contradictions elsewhere. The local draft and live Overleaf project have not been replaced by this fragment.

## Wait for calibration before updating numerical conclusions

Do not transfer historical p values, ICCF delays, conditional results or estimator comparisons to the new 242 spectrum selection. Recompute the analyses that remain in the paper using the adopted inputs, or explicitly identify them as historical diagnostic calculations.

Do not describe the pilot peaks as detected lags. The independent simulations also produce substantial information maxima, and native histograms are sparse. The eight simulation runs measure runtime and check the execution path; they do not measure a reliable false positive rate or detection power.

## Work remaining beyond sampling

1. Define the scientific null and the family of scanned statistics, targets, lags and estimator settings. Calibrate the sampling comparison using separate calibration and evaluation seeds.
2. Measure false positive rates, detection power, lag errors and intended atom identification rates with uncertainty intervals. Vary response strength, width, noise and process timescale. Include a continuum only baseline to assess continuum response delays.
3. Test whether another velocity component adds information beyond continuum and target history. Controls must preserve the dependencies allowed by that null. Independent channel shifts do not automatically test this incremental question.
4. Assess estimator sensitivity, including histogram sparsity and exclusion of nearby times in continuous neighbour estimators. Any headline prediction claim needs improvement on observing seasons held out from fitting against simple baselines.
5. Repeat retained comparison analyses on the adopted sample. Treat Mrk 817 as an external dependence check unless the same decomposition and scientific test are actually replicated.
6. Update the abstract, conclusions, figures and captions only after the relevant results are validated. Choose one authoritative manuscript, map every displayed number to a saved output, compile and visually review the final paper, and verify the reproduction instructions from a clean checkout.

The immediate computational step is the null specification and continuum baseline, followed by a modest calibration batch on the laptop. The benchmark supports local execution for that scale. A larger simulation grid or nested calibration can then move to institutional computing resources.

## Evidence

The adopted manifests and data table are in `agn_surd_project/processed/reconciled/adopted`. Sampling counts are in `agn_surd_project/processed/reconciled/sampling`. The executed pilot is documented in `pipeline/BENCHMARK.md`. Archive conventions are checked against the [AGN Watch data descriptions](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/) and the [Wanders and Peterson source paper](https://articles.adsabs.harvard.edu/pdf/1996ApJ...466..174W), especially sections 3.1 and 3.2 and corrected Table 1.
