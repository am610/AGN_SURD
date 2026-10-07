# Published flux run and bin sensitivity

The corrected published Table 1 supplies total, blue, core, and red line
measurements for every one of the 242 adopted spectra. These are now prepared
on 224 native dates, alongside the original F5100 continuum. No fitted flux
correction or R squared analysis was introduced. The previous release and
earlier exploratory run remain unchanged.

Source paper section 3.2 on printed page 180 explicitly gives observed
wavelength intervals 4870 through 4920, 4920 through 4960, and 4960 through
5010 Angstrom. Its figure caption differs on the blue lower boundary. The
preparation uses the corrected published measurements rather than claiming
an exact numerical reconstruction from the profile pixels. Shared spectra
and the near sum relationship between total and components remain important
interpretive limitations.

Both two bin and three bin runs use the unchanged authors' SURD routine and
the same campaign intervals, target times, ten day maximum interpolation gap,
zero through thirty day base lags, and 32 tuple numerical guard. Unequal five
day blue and red shifts remain stipulated sensitivity scenarios, not measured
delays and not the expanded predictor construction of equation 3.18.

Each run contains 3720 configurations, of which 1402 pass the numerical guard.
Each usable configuration has 26 components. Maximum normalization errors are
9.99e-16 for two bins and 1.44e-15 for three bins. These are software checks.

For 1326 matched positive lag configurations, the median total variation
between normalized component distributions is 0.483783 and the maximum is
0.851933. Total variation here is one half of the sum of the absolute changes
across the 26 components. About 98.87 percent of these configurations exceed
0.2. Median absolute normalized leakage change is 0.181997. These values are
descriptive diagnostics, not a new hypothesis test or confidence statement.

Consequently neither bin setting is selected as physically reliable merely
because it normalizes or gives an interesting component maximum. The next
scientific step is to review component stability and sampling support by
campaign before interpreting maxima or executing expanded predictor scans.
A defensible results section may need to report estimator sensitivity and
sampling limitations rather than a new astrophysical inference.

Reproduce from the repository root:

```sh
python -m collaborator_rebuild.prepare_historical
python -m collaborator_rebuild.run_historical --published --bins 2 --no-plots
python -m collaborator_rebuild.run_historical --published --bins 3 --no-plots
python -m collaborator_rebuild.compare_bins
```

Prepared inputs and derived tables are local in `historical_prepared` and
`published_results`. Original continuum retains host starlight, and duplicate
uncertainties omit unavailable shared calibration covariance. Those limitations
are preserved in the preparation manifest.
