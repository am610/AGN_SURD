# SURD for AGN reverberation mapping

The active study follows the September collaborator feedback: prepare five
observed series by campaign, use the authors' SURD implementation, show every
information component and leakage, and deliver a scoped LaTeX results section.
Historical NGC 5548 observations are the authorized first dataset. Lu and Xi
campaigns are a possible later extension.

## Current status

The historical run is exploratory, not a completed astrophysical result.
Campaign assignments and decomposition arithmetic have been checked. Total
H beta now uses corrected published flux columns; exact profile integration
reconstruction and estimation reliability remain unresolved. Equation 3.18 predictor
construction is implemented and tested; expanded scientific results are not yet
validated. A separate descriptive results section is drafted under
`collaborator_rebuild/results_section`, with a compiled preview at
`output/pdf/surd_campaign_results.pdf`. Scientific review and calibration
remain pending.

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
| `main.tex`, `main.pdf`, `refs.bib` | Previous released manuscript, not the new deliverable |

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
