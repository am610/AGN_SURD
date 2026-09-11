# Local information benchmark

Measured on the current laptop using Python 3.11.5 and two worker processes.

The observed scan evaluated four sampling methods, three line targets and 201 lags, producing 2412 rows in 2.00 seconds. Eight simulated datasets completed in 8.82 seconds of elapsed time. Four seeds were evaluated for independent signals and four for a delayed shared driver.

The parent process peaked at 147.5 MiB. The largest worker peak was 160.0 MiB. These are separate process maxima, not a simultaneous system memory measurement.

At the measured throughput, 1000 repetitions of this exact scan would take approximately 18.4 minutes. This is an extrapolation from a small batch. Thermal effects, competing laptop work, different bin counts and larger estimator families may change it. An outer loop of 1000 datasets with 1000 null scans per dataset would instead imply approximately 306 hours at this throughput before further overhead.

Recommendation: keep the next modest calibration batch on the laptop with two workers. Use Midway or NERSC once the declared calibration and power grid requires substantially more repetitions or a nested null calculation. No remote job or cloud resource was started.

The output directory is `agn_surd_project/processed/reconciled/information_pilot/04dde4b289c6ac42`. `benchmark.json` contains fresh timings, source hashes, settings and simulation seeds. `benchmark_resume.json` confirms all eight completed simulations were reused after output hash verification. Ten tests passed, including independent variables, XOR synergy, duplicated predictor redundancy, direct mutual information agreement, constant targets, exclusion of the target predictor, simulation reproducibility and sampling boundaries.

## Scientific status

These curves are exploratory estimator outputs. No calibrated p values, detection rates or astrophysical lag detections are reported. Independent simulated signals also produce substantial raw information peaks. Native sampling has a median empty histogram fraction of 0.848, and its median fraction of tuples in singleton cells is 0.389. Calibration must reproduce this sparse support.

The shared driver case is not a pure test of a 20 day joint information peak: other line predictors share the same delayed response and can predict one another at zero lag. A continuum only baseline should be included when measuring continuum lag recovery.

Next: specify the null and statistic family, add a continuum only lag recovery baseline, and separate calibration seeds from evaluation seeds. Then measure false positive rates and power with uncertainty intervals under the same sampling rules. Expand noise, response strength, response width and process timescale settings before drawing conclusions about the real observations.
