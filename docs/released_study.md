# Information decomposition for AGN reverberation mapping

This repository investigates the behavior of Synergistic Unique Redundant Decomposition in irregular observations of NGC 5548. The current work validates data provenance and quantifies how sampling and information estimation affect the results. It does not establish a new astrophysical lag or causal relationship.

The main result is deliberately narrow. SURD provides a richer partition of predictive information than conventional lag estimators, but no information atom in these observations is both globally significant and robust enough to establish additional astrophysical knowledge. The prompt relationships are recovered more directly by ICCF and simpler information measures. The demonstrated contribution is therefore a validation framework and a warning against interpreting uncalibrated multivariate information structure as AGN physics.

## Current analysis

The adopted pipeline uses 242 Hβ spectra on 224 dates and 557 original continuum observations. Five profiles without exact date and instrument matches in the corrected observation table are excluded from the adopted sample and preserved in the archive comparison. The revised continuum is available as a sensitivity input.

Four sampling methods have been implemented across 201 lags. Common-lag and unequal-lag decomposition now cover every target in both requested four-variable sets. The continuum baseline uses 999 scan-wide null surrogates and controls the three-target family. The power grid uses separate null calibration and evaluation samples, reports uncertainty intervals, and includes one-factor width, timescale, and noise sensitivities. Incremental scans use standardized KSG inputs and 100 surrogates. Their two-segment prediction check is diagnostic only.

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
python -m pipeline.continuum_baseline
python -m pipeline.calibrate_significance --surrogates 999 --workers 4
python -m pipeline.decompose_all_atoms
python -m pipeline.asymmetric_lags
python -m pipeline.detection_power_grid --reps 200 --workers 4
python -m pipeline.incremental_information --surrogates 100 --workers 4
python -m pipeline.final_closure_checks
```

Input retrieval downloads public archive data and the upstream SURD implementation at commit `79dbdea85e6754ec2b5457b3e37204c5d53d1815`. It verifies file hashes, refuses to overwrite differing local files, and extracts only the expected spectra. Raw archive inputs and the source paper are retrieved from their providers rather than distributed here. The corrected table text required by the loader is included under `audit` with its provenance documented.

## Documentation and saved evidence

1. [Pipeline definitions and limitations](../pipeline/README.md)
2. [Execution status and next steps](../pipeline/PROGRESS.md)
3. [Measured laptop benchmark](../pipeline/BENCHMARK.md)
4. [Paper revision notes](../overleaf_draft/REVISION_NOTES.md)
5. [Adopted data and manifests](../agn_surd_project/processed/reconciled/adopted)
6. [Sampling support audit](../agn_surd_project/processed/reconciled/sampling)
7. [Executed information pilot](../agn_surd_project/processed/reconciled/information_pilot/04dde4b289c6ac42)

## Authoritative release files

`main.tex` is the authoritative manuscript source and `main.pdf` is its compiled
release artifact.  The matching files under `overleaf_draft` are synchronized
copies for Overleaf packaging, not independent manuscript versions.  Historical
notebooks and versioned PDFs are exploratory records and are not inputs to the
adopted pipeline.

The release is considered reproducible only after the commands above have run
from a clean checkout, the unit tests pass, regenerated numerical outputs have
been compared with the committed evidence, and the compiled PDF has received a
complete visual review.

The information pilot includes curves and peak summaries for all eight simulated datasets, observed curves and peaks, timings, settings and source hashes. Its peaks are descriptive; no calibrated significance is claimed. The benchmark estimates about 18 minutes for 1000 repetitions of that exact scan on the tested laptop with two workers. Expanded parameter grids and nested calibration can require much more computation.

## Final closure status

The five release tasks are now addressed. The central power grid was rerun with 200 calibration and 200 evaluation null realizations per run, for 3400 central realizations. Intended atom identification is recorded in the 100 realization synthetic validation. Temporal neighbour exclusion was tested at 1, 3, and 5 days. A four block contiguous date validation was added, with two blocks meeting the strict overlap requirement at the 15 day evaluation lag. The clean checkout workflow and 22 page PDF review have passed.

No further exploratory search for a positive SURD result is recommended with this dataset. A credible test of added astrophysical value requires a predefined analysis using denser simultaneous optical, ultraviolet, X ray, and velocity resolved observations, ideally with an estimator independent decomposition. This is future work rather than a release blocker.

## Sources

The decomposition is provided by the [upstream SURD project](https://github.com/Computational-Turbulence-Group/SURD), under its own license. Archive observations are supplied by [AGN Watch](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/spectra/). Source papers and data attribution are recorded in the pipeline manifests. This repository supplies the astronomical validation workflow and retains earlier exploratory analyses for context.
