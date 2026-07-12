# AGN Reverberation Mapping Data Research Notes

## NGC 5548 (AGN Watch)
- **Continuum (5100 Å):** [https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/OPT1988/c5100.dat](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/OPT1988/c5100.dat)
- **H-beta (Integrated):** [https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/OPT1988/hb.dat](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/OPT1988/hb.dat)
- **UV Continuum (1350 Å):** [https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/IUE1989/c1350.sips](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/IUE1989/c1350.sips)
- **Format:** 
    - Col 1: JD - 2440000
    - Col 2: Flux (Continuum: 10^-15 ergs/s/cm^2/A, Lines: 10^-13 ergs/s/cm^2)
    - Col 3: Uncertainty
- **Velocity-resolved:** Need to check if velocity-resolved data is available for NGC 5548 in AGN Watch or STORM.

## AGN STORM (NGC 5548)
- **STORM Data:** Usually hosted on MAST or specific project sites.
- **STORM 2:** [https://archive.stsci.edu/hlsp/storm2](https://archive.stsci.edu/hlsp/storm2) (FITS format light curves).

## SDSS-RM
- **DR18/DR19:** [https://www.sdss.org/dr18/bhm/programs/rm/](https://www.sdss.org/dr18/bhm/programs/rm/)
- **Data Release:** Shen et al. 2024 (final data).

## OzDES
- **OzDES RM:** Usually part of DES releases.

## AGN STORM 1 (NGC 5548)
- **Emission Lines Table (deblended):** [http://archive.stsci.edu/hlsps/storm/hlsp_storm_hst_cos_ngc-5548_g130m-g160m_v1_emission-lines.txt](http://archive.stsci.edu/hlsps/storm/hlsp_storm_hst_cos_ngc-5548_g130m-g160m_v1_emission-lines.txt)
- **Description:** Contains deblended emission line fluxes (Ly alpha, N V, Si IV, C IV, He II) for each visit.

## AGN STORM 2 (Mrk 817)
- **Light Curve FITS:** [https://archive.stsci.edu/hlsps/storm2/hlsp_storm2_hst_cos_mrk-817-go16196_g130m-g160m_v1_data-files.tar.gz](https://archive.stsci.edu/hlsps/storm2/hlsp_storm2_hst_cos_mrk-817-go16196_g130m-g160m_v1_data-files.tar.gz)
- **Description:** Contains a single FITS file with multiple extensions for continuum and emission lines (Ly alpha, N V, Si IV, C IV, He II).
- **Extensions:**
    - ext1: CONTINUUM-1180
    - ext2: LYA
    - ext7: CIV
    - ext8: HEII

## NGC 5548 Velocity-Resolved Data (AGN Watch)
- **H-beta Profiles (1989-1996):** [https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/spectra/Profiles/all.tar.gz](https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/spectra/Profiles/all.tar.gz)
- **Description:** 247 high-quality optical spectra with continuum and narrow lines subtracted. 3-column ASCII (Wavelength, Flux, Uncertainty).
- **Note:** Suitable for splitting into blue wing, core, and red wing.
