# Research Summary: AGN Reverberation Mapping Datasets for SURD Analysis

This document summarizes the datasets prepared for the Synergistic–Unique–Redundant Decomposition (SURD) causal information analysis of Active Galactic Nuclei (AGN).

## 1. Datasets Overview

We have prepared three primary datasets from well-known AGN monitoring campaigns, focusing on high-quality light curves suitable for measuring information flow between continuum and emission lines.

### A. NGC 5548 (AGN Watch)
- **Campaign:** International AGN Watch (1988–2001)
- **Continuum:** Optical 5100 Å
- **Emission Line:** Integrated Hβ flux
- **Characteristics:** Long baseline (13 years), high signal-to-noise ratio.
- **Data Source:** [AGN Watch Archive](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/)

### B. NGC 5548 (STORM 1)
- **Campaign:** Space Telescope and Optical Reverberation Mapping (2014)
- **Continuum:** UV (1158 Å, 1367 Å, 1469 Å, 1745 Å)
- **Emission Lines:** Deblended Lyα, N V, Si IV, C IV, He II
- **Characteristics:** Exceptional daily cadence over 6 months, UV-focused.
- **Data Source:** [MAST HLSP STORM](https://archive.stsci.edu/hlsp/storm)

### C. Mrk 817 (STORM 2)
- **Campaign:** STORM 2 (2020–2022)
- **Continuum:** UV 1180 Å
- **Emission Lines:** Lyα, C IV, He II
- **Characteristics:** Recent high-cadence monitoring of a highly variable AGN.
- **Data Source:** [MAST HLSP STORM 2](https://archive.stsci.edu/hlsp/storm2)

## 2. File Structure and Formats

All datasets have been processed into clean CSV files with standardized column names:

| Dataset | File Name | Columns | Units |
|---------|-----------|---------|-------|
| NGC 5548 Watch | `ngc5548_agnwatch_cont_clean.csv` | `mjd`, `flux`, `error` | Flux: 10⁻¹⁵ erg/s/cm²/Å |
| NGC 5548 Watch | `ngc5548_agnwatch_hb_clean.csv` | `mjd`, `flux`, `error` | Flux: 10⁻¹³ erg/s/cm² |
| STORM 1 | `ngc5548_storm1_clean.csv` | `mjd`, `f_lya`, `e_f_lya`, ... | Flux: 10⁻¹³ erg/s/cm² |
| STORM 2 | `mrk817_storm2_continuum-1180_clean.csv` | `hjd-2400000`, `flux`, `error` | Flux: 10⁻¹³ erg/s/cm²/Å |

## 3. Velocity-Resolved Capability

For NGC 5548, we have downloaded the **Hβ Profiles (1989–1996)**. These consist of 247 high-quality spectra that can be integrated over specific velocity windows:
- **Blue Wing:** -6000 to -2000 km/s
- **Core:** -2000 to +2000 km/s
- **Red Wing:** +2000 to +6000 km/s

This allows for the analysis of information flow as a function of gas kinematics in the Broad Line Region (BLR).

## 4. Usage Instructions

1. **Load Data:** Use `pd.read_csv()` or the provided `agn_helpers.py`.
2. **Resample:** Use `prepare_surd_input` to interpolate continuum and line data onto a common grid.
3. **Standardize:** The helpers automatically apply zero-mean unit-variance scaling.
4. **SURD Analysis:** Feed the resulting `X = np.vstack([cont, line])` into your SURD implementation.
