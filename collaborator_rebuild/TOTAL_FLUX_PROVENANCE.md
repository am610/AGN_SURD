# Total broad H beta provenance check

The collaborator feedback requests SURD analysis, not an R squared regression
or a fitted flux recalibration. No additional regression algorithm or fitted
calibration was introduced. The user authorized the most appropriate next
step, interpreted as following the collaborator's SURD scope.

The current adopted total sums the three velocity windows. A diagnostic now
compares that sum with integration over the interval recorded from source paper
section 3.2, 4870 through 5010 Angstrom, its caption alternative beginning at
4840 Angstrom, and the complete profile file. These intervals are not selected
by optimizing agreement with the archive.

For all 224 adopted dates, the section 3.2 interval has correlation 0.999967
with the revised archive broad H beta daily means. Its maximum absolute
relative difference is 1.91 percent. The current velocity window sum has
correlation 0.998730 but a varying scale discrepancy. Integrating the entire
file worsens correspondence, so the file span cannot simply be equated with
the intended broad line integration interval.

These observations support further investigation of the source paper interval,
not an automatic flux rescaling. The archive contains extra measurements and
cannot be matched to individual instruments using dates alone. The profile
archive labels its units arbitrary, whereas the integrated line archive gives
physical units. Close correspondence does not establish absolute calibration.

Source section and figure caption disagree on the lower boundary. Confirm the
printed definition before adopting a revised integration policy. Decide total
and component boundaries together; a deliberately partitioned total can have
an exact algebraic relation to its components, which must be acknowledged in
the interpretation rather than treated as independent evidence.

The printed section 3.2 definition is now visually confirmed on page 180.
More importantly, the corrected table already supplies all four measured line
columns. `prepare_historical.py` uses those published values directly, with
instrument and date checks, rather than fitting a scale or selecting an
integration interval by numerical agreement. See `BINNING_REVIEW.md` for
the resulting descriptive scans and their substantial bin sensitivity.

Run `python -m collaborator_rebuild.check_support` to retrieve the candidate
archive series, then `python -m collaborator_rebuild.audit_total_flux` to
generate per campaign diagnostics and provenance. No active scientific inputs
or previous release data are replaced by either check.

Primary archive descriptions:
https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/spectra/
https://www.asc.ohio-state.edu/astronomy/agnwatch/n5548/lcv/
