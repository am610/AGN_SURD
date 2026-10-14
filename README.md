# SURD for AGN reverberation mapping

The active study follows the September collaborator feedback: prepare five
observed series by campaign, use the authors' SURD implementation, show every
information component and leakage, and deliver a scoped LaTeX results section.
Historical NGC 5548 observations are the authorized first dataset. Lu and Xi
campaigns are a possible later extension.

## Current status

The September instructions are answered by the per season analysis in
[collaborator_rebuild/EPOCH_ANALYSIS.md](collaborator_rebuild/EPOCH_ANALYSIS.md): five
historical seasons, both four variable sets, every target, 26 components and
leakage against lag, common and uncommon predictor lags, 7200 decompositions and
120 figures. The results section is `collaborator_rebuild/results_section/results_section.tex`
and its compiled review copy is `output/pdf/surd_epoch_results.pdf`. The replies to
each instruction are in `output/docx/SURD_September_Instructions_With_Replies.docx`.

## Earlier review boundary (native date analysis, 7 October)

The bounded diagnostic work is complete and the analysis is frozen for
collaborator review. Initial categorical calibration and the predictor ordering
audit are complete. The separate results section now has five tables, four
figures, and a ten page PDF. All pages were visually inspected.

Ordering changes component allocation in 11 of 240 configurations. The maximum
total variation is 0.187814. Total information and leakage remain unchanged to
numerical precision. The 1993 common curves are stable under ordering but remain
strongly sensitive to bin choice. No stable physical delay is established.

Read [review handoff](collaborator_rebuild/REVIEW_HANDOFF.md) for completed requirements, remaining
gaps, and the stopping rule. Further scientific runs or changes to the authors'
routine require a new agreed scope. Earlier progress entries below are history.

Start with [current progress](collaborator_rebuild/PROGRESS.md),
[analysis specification](collaborator_rebuild/SPECIFICATION.md), and
[review findings](collaborator_rebuild/REVIEW.md).

## Repository layout

| Location | Role |
| :--- | :--- |
| `collaborator_rebuild/` | Active study code, specification, and status |
| `pipeline/` | Shared data provenance tools and previous validation pipeline |
| `agn_surd_project/` | Existing data workspace and saved release evidence |
| `audit/` | Corrected observation table and provenance evidence |
| `docs/` | Repository guidance and previous release documentation |
| `history/` | Earlier exploratory notebooks, documents, figures, and utilities |
| `SURD/` | Locally fetched, unchanged official SURD source, ignored by Git |
| `main.tex`, `main.pdf`, `refs.bib` | The paper: per season SURD analysis of NGC 5548. Compile `main.tex` (MNRAS class). It inputs `collaborator_rebuild/results_section` and the figure atlases in `collaborator_rebuild/epoch_results/figures`. The earlier validation study is the git tag `validation-study-main-tex` |

## Setup and verification

Use Python 3.11 and an isolated environment. Existing numerical dependencies
are recorded in `pipeline/requirements.txt`; active test dependencies are in
`requirements-dev.txt`.

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python -m pipeline.fetch_inputs
python -m unittest pipeline.test_sampling pipeline.test_information
python -m pytest -q collaborator_rebuild/test_lag_design.py
```

Input retrieval verifies hashes and preserves differing local files. Raw data
and third party source documents are retrieved from providers, not republished
as new project code. The authors' SURD implementation is pinned at commit
`79dbdea85e6754ec2b5457b3e37204c5d53d1815`.

## Active commands

Run from the repository root. These commands write local derived products.

```sh
python -m collaborator_rebuild.check_support
python -m collaborator_rebuild.run_historical
python -m collaborator_rebuild.plot_reviewed
```

The historical runner is descriptive. It does not resolve total flux provenance
or establish reliable information estimates. The support checker builds the
expanded predictor vector but does not report an expanded scientific result.
The plot renderer uses the existing historical run, so rerun the historical
runner first if its inputs or settings change.

See [repository policy](docs/REPOSITORY.md) for tracked and local artifacts.
See [previous release documentation](docs/released_study.md) for the older
study and its reproduction commands. Its completion statements apply to that
release only, not to the active collaborator study.
