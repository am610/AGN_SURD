# Data reconciliation

## Adopted execution path

Run these commands from the project root:

```sh
python pipeline/reconcile_data.py
python pipeline/adopted_data.py
python -m unittest pipeline.test_sampling
python pipeline/sampling.py
```

The first command preserves the earlier candidate reduction. The second writes the adopted sample to `agn_surd_project/processed/reconciled/adopted`. Use `load_adopted` from `pipeline.adopted_data` for subsequent analysis. The final command writes sampling support counts and the exact configuration to `agn_surd_project/processed/reconciled/sampling`.

The adopted sample excludes five spectra without an exact date and instrument match in corrected Table 1. It contains 242 spectra on 224 dates. All 247 original profiles remain available in the candidate output. This is a documented provenance exclusion, not a claim that the five measurements are invalid. Their dates also occur among the five additional revised continuum rows absent from the table.

The primary continuum remains the original archive series, with 557 observations in the selected window. This preserves continuity with the existing analysis and retains denser native sampling. The revised continuum is exported separately as a sensitivity input, with 242 measurements on 224 dates. Its native flux units differ by a factor of ten from the original continuum units; the two series also have different calibration and galaxy subtraction conventions. They must not be silently substituted or concatenated.

The revised continuum file actually uses JD less 2440000 despite the archive landing page stating JD less 2400000. The loader adds 40000, then verifies the dates and continuum fluxes against every corrected table row. The flux comparison allows 0.00051 in the revised file units because the table prints three decimal places and one discrepancy is 0.000504. It requires a unique match without reusing any row.

Wanders and Peterson section 3.2 specifies observed wavelengths. It gives redshifts 0.0165 for the broad peak, 0.0174 for narrow lines, and 0.0172 for the radio reference. The existing value 0.017175 is retained as an explicit numerical convention, not attributed to that paper. The project velocity bins are likewise a continuity choice. The paper gives wavelength intervals from 4870 to 4920, 4920 to 4960, and 4960 to 5010 Å in section 3.2. Its Figure 1 caption differs in the first lower boundary. Section 3.1 describes the revised continuum window as centred at observed 5188 Å, although the table and archive label it 5150 Å. These discrepancies are preserved in the manifest.

Each series retains its actual support endpoints. Pairing enforces those supports at the requested lag. Repeatedly cropping all series to the last retained measurement would incorrectly discard valid delayed targets.

## Sampling status

The saved sampling configuration specifies daily interpolation, a 30 day maximum interpolation gap, seasons divided at gaps exceeding 60 days in either input, and native nearest observation pairing with a one day tolerance. Ties select the earlier observation. Native pairing uses at most one tuple per target date but can reuse predictor dates; the report counts that reuse. Positive lags place the target after its predictors. Seasonal analysis also applies the 30 day interpolation rule.

These settings precede new significance calculations but follow the existing cadence audit, so this is an exploratory comparison. The common sampler accepts replacement fluxes on the same observation times for simulation calibration. The current output contains support counts only. Information scans and a small simulation runtime pilot are now implemented. Null calibration, recovery rates and uncertainty intervals remain open. At zero lag the target component also appears in the predictor columns; inference code must remove the contemporaneous target predictor before estimating dependence.

## Candidate reduction retained for comparison

Run `python pipeline/reconcile_data.py` from the project root. Outputs are written to `agn_surd_project/processed/reconciled`. Original inputs and existing research outputs are preserved.

The candidate reduction retains the velocity intervals that reproduce the saved analysis: blue from −6000 to −2000, core from −2000 to 2000, and red from 2000 to 6000 km/s. These are an explicit continuity choice, not a claim that they match the bins in another study. Pixel centres determine membership. The adopted wavelength reference is 4861.33 Å with redshift 0.017175, inherited from the existing audit.

The archive labels continuum dates as JD less 2400000. Calling this MJD introduces a half day coordinate error. The spectral filename date rule is inherited from the original reduction and still requires confirmation against observation logs. The precise time scale and any heliocentric correction are unspecified by the archive landing page.

The Hβ profile archive supplies arbitrary flux density units. Integrating each selected pixel over its 2 Å width gives arbitrary integrated flux units. The old pipeline omitted this constant width. The new component fluxes and errors therefore equal twice the saved values. This constant scaling alone does not change histogram information estimates with correspondingly rescaled bin edges. Duplicate handling can change the estimates.

Repeated dates use arithmetic means with the same weights for every component. This preserves the identity between the component sum and `profile_total`. Errors use quadrature propagation under independent pixel and exposure assumptions. Shared calibration errors and pixel covariance are unavailable and remain a limitation. The component sum is not the archive total Hβ series, whose narrow line contribution was retained.

The native tables are not interpolated. `data_table.csv` reports actual unique dates and source observation counts. Continuum and line support endpoints differ because continuum filtering retains only actual measurements inside the overlap bounds; subsequent pairing must enforce both supports.

`manifest.json` records source hashes, conventions, dependencies, and the residual spectral timing issue. The runnable checks compare every spectrum against the previous reduction and check component sum preservation after date aggregation.

Archive documentation:

[Optical light curves](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/)

[Hβ profiles](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/spectra/)


## Information scans and local runtime pilot

Run from the project root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest pipeline.test_sampling pipeline.test_information
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pipeline.benchmark_information --workers 2
```

The runner reads `sampling_config.json` and `information_config.json`. It evaluates three line targets with continuum and the other two line components as predictors. The target component never predicts itself, even at positive lags. This is a dependence scan without target history conditioning. It is not the incremental prediction test in Task 4 and does not include a continuum target.

The exploratory estimator uses four marginal quantile bins based on native observations. Edges remain fixed across sampling methods and lags within each dataset and are recomputed by the identical rule for each simulation. Duplicate quantiles reduce the available bins. Histograms with fewer than 32 tuples are skipped; this threshold is a numerical guard, not an assertion of adequate sample size. Constant targets are skipped. Every usable row records joint mutual information, individual SURD atoms, total synergy in bits, synergy divided by joint information, target entropy and histogram occupancy. Information accounting is checked for every decomposition. The ratio is undefined when joint information is essentially zero.

The pilot uses four seeds in each of two stipulated process models: independent correlated signals and a shared continuum driver with a smoothed delayed line response plus independent line variability. Signals are sampled on the adopted native dates and perturbed using the existing propagated errors before entering the common sampler. Noise is applied after date aggregation. Shared calibration covariance and fitted process uncertainty are absent. This is a correctness and runtime pilot, not a fitted null model or power measurement. The shared driver also correlates the line predictors at zero lag, so a joint information peak need not recover the continuum response delay.

Each run writes to a directory identified by source and configuration hashes inside `agn_surd_project/processed/reconciled/information_pilot`. Outputs include observed curves and peaks, curves and peaks for every simulated dataset, and `benchmark.json`. Completed simulations have a completion marker with output hashes and can be reused on restart. A restart writes `benchmark_resume.json` without replacing the original fresh timing report. Input hashes are checked against the adopted manifest before execution. Changing code or settings starts a separate output directory.

The runtime projection covers repetitions of this exact scan only. It does not include a larger parameter grid, another estimator, or an outer simulation loop with a separately calibrated null distribution for each dataset. Those extensions multiply the workload. See `BENCHMARK.md` for the measured local resource recommendation.
