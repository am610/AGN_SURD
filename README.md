# Information decomposition for AGN reverberation mapping

This repository investigates the behavior of Synergistic Unique Redundant Decomposition in irregular observations of NGC 5548. The current work validates data provenance and quantifies how sampling and information estimation affect the results. It does not establish a new astrophysical lag or causal relationship.

## Current analysis

The adopted pipeline uses 242 Hβ spectra on 224 dates and 557 original continuum observations. Five profiles without exact date and instrument matches in the corrected observation table are excluded from the adopted sample and preserved in the archive comparison. The revised continuum is available as a sensitivity input.

Four sampling methods have been implemented across 201 lags. Observed information scans and a small simulation pilot have run successfully. Significance calibration, false positive rates, detection power and incremental prediction tests remain open. Historical scripts, notebooks and manuscript versions elsewhere in the repository use earlier selections and settings; their results must not be attributed to the adopted pipeline.

## Reproduce the current pipeline

Use Python 3.11 in a separate environment. The requirements pin the versions used for the published pilot.

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r pipeline/requirements.txt
python pipeline/fetch_inputs.py
python pipeline/reconcile_data.py
python pipeline/adopted_data.py
python -m unittest pipeline.test_sampling pipeline.test_information
python pipeline/sampling.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pipeline.benchmark_information --workers 2
```

Input retrieval downloads public archive data and the upstream SURD implementation at commit `79dbdea85e6754ec2b5457b3e37204c5d53d1815`. It verifies file hashes, refuses to overwrite differing local files, and extracts only the expected spectra. Raw archive inputs and the source paper are retrieved from their providers rather than distributed here. The corrected table text required by the loader is included under `audit` with its provenance documented.

## Documentation and saved evidence

1. [Pipeline definitions and limitations](pipeline/README.md)
2. [Execution status and next steps](pipeline/PROGRESS.md)
3. [Measured laptop benchmark](pipeline/BENCHMARK.md)
4. [Paper revision notes](overleaf_draft/REVISION_NOTES.md)
5. [Adopted data and manifests](agn_surd_project/processed/reconciled/adopted)
6. [Sampling support audit](agn_surd_project/processed/reconciled/sampling)
7. [Executed information pilot](agn_surd_project/processed/reconciled/information_pilot/04dde4b289c6ac42)

The information pilot includes curves and peak summaries for all eight simulated datasets, observed curves and peaks, timings, settings and source hashes. Its peaks are descriptive; no calibrated significance is claimed. The benchmark estimates about 18 minutes for 1000 repetitions of that exact scan on the tested laptop with two workers. Expanded parameter grids and nested calibration can require much more computation.

## Remaining research

1. Complete the sampling comparison with a declared null, statistic family and significance calibration.
2. Measure false positive rates, detection power, lag recovery and atom identification with uncertainty intervals.
3. Test whether another velocity component adds information beyond continuum and target history, using controls that preserve the dependencies allowed by that null.
4. Check estimator settings and temporal neighbour exclusion. Support any predictive claim with evaluation on observing seasons held out from fitting.
5. Align manuscript text, tables and figures with the final validated pipeline and verify the complete release.

## Sources

The decomposition is provided by the [upstream SURD project](https://github.com/Computational-Turbulence-Group/SURD), under its own license. Archive observations are supplied by [AGN Watch](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/spectra/). Source papers and data attribution are recorded in the pipeline manifests. This repository supplies the astronomical validation workflow and retains earlier exploratory analyses for context.
