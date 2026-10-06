"""Prepare only the public Lu et al. native light curves by observing season."""

import hashlib
import json
import io
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "collaborator_rebuild"
SOURCE = BASE / "sources" / "lu_2022_table2.dat"
OUT = BASE / "prepared"
EXPECTED_SHA256 = "5e35a557eeacd0f11242bf0edfad384030f20219a7ca3571c94f2adedadd72f7"
EXPECTED_SEASONS = [2015, 2018, 2019, 2020, 2021]
CAMPAIGNS = {
    2015: ('2015-01-07', '2015-08-01', 62),
    2018: ('2018-03-12', '2018-06-18', 40),
    2019: ('2018-11-28', '2019-06-20', 81),
    2020: ('2020-01-11', '2020-06-21', 52),
    2021: ('2020-12-24', '2021-08-06', 80),
}
XI_SOURCE = BASE / 'sources' / 'xi_2025_table1.txt'
XI_SHA256 = '04bf0f8940d11d5f30878458fcf9b04cb5899bd0f8931f0e034fc50f131eb99a'


def read_lu_table(path=SOURCE):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError(f"Unexpected Lu table hash: {digest}")
    names = [
        "jd",
        "f5100",
        "f5100_error",
        "fheii",
        "fheii_error",
        "fhei",
        "fhei_error",
        "fhgamma",
        "fhgamma_error",
        "fhbeta",
        "fhbeta_error",
    ]
    frame = pd.read_csv(path, sep=r"\s+", names=names)
    frame = frame[names].apply(pd.to_numeric, errors="raise")
    dates = pd.to_datetime(frame.jd, unit='D', origin='julian').dt.normalize()
    frame['season'] = 0
    for season, (start, end, expected) in CAMPAIGNS.items():
        selected = dates.between(pd.Timestamp(start), pd.Timestamp(end))
        if selected.sum() != expected:
            raise ValueError(f'Campaign {season} count differs from Lu Table 1')
        frame.loc[selected, 'season'] = season
    if (frame.season == 0).any():
        raise ValueError('Observation outside published Lu campaign boundaries')
    return frame


def validate(frame):
    if len(frame) != 315:
        raise ValueError(f"Expected 315 observations, found {len(frame)}")
    if sorted(frame.season.unique().tolist()) != EXPECTED_SEASONS:
        raise ValueError(f"Unexpected observing seasons: {sorted(frame.season.unique())}")
    if frame.jd.duplicated().any() or not frame.jd.is_monotonic_increasing:
        raise ValueError("Julian Dates must be unique and increasing")
    if (frame[["f5100_error", "fhbeta_error"]] <= 0).any().any():
        raise ValueError("Uncertainties must be positive")


def main():
    frame = read_lu_table()
    validate(frame)
    if hashlib.sha256(XI_SOURCE.read_bytes()).hexdigest() != XI_SHA256:
        raise ValueError('Unexpected Xi table hash')
    data_lines = [line for line in XI_SOURCE.read_text().splitlines()
                  if line.startswith('24')]
    xi = pd.read_csv(io.StringIO('\n'.join(data_lines)), sep=r'\s+',
                     names=frame.columns[:-1].tolist())
    if len(xi) != 74 or xi.isna().any().any():
        raise ValueError('Expected 74 complete Xi observations')
    if xi.jd.duplicated().any() or not xi.jd.is_monotonic_increasing:
        raise ValueError('Invalid Xi dates')
    if (xi[['f5100_error', 'fhbeta_error']] <= 0).any().any():
        raise ValueError('Invalid Xi uncertainties')
    xi['season'] = 2023
    frame = pd.concat([frame, xi], ignore_index=True)
    seasons = EXPECTED_SEASONS + [2023]
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for season in seasons:
        selected = frame.loc[frame.season == season]
        season_dir = OUT / str(season)
        season_dir.mkdir(parents=True, exist_ok=True)
        for name, value, error in [
            ("f5100", "f5100", "f5100_error"),
            ("hbeta_total", "fhbeta", "fhbeta_error"),
        ]:
            output = selected[["jd", value, error]].copy()
            output.columns = ["jd", "value", "error"]
            output.to_csv(season_dir / f"{name}.csv", index=False)
            rows.append({"season": season, "series": name, "observations": len(output),
                         "first_jd": float(output.jd.min()),
                         "last_jd": float(output.jd.max())})
    pd.DataFrame(rows).to_csv(OUT / "index.csv", index=False)
    manifest = {
        "source": "Lu et al. 2022 Table 2 and Xi et al. 2025 Table 1",
        "doi": ["10.3847/1538-4365/ac94d3", "10.3847/1538-4357/ae1ccb"],
        "source_sha256": EXPECTED_SHA256,
        "xi_source_sha256": XI_SHA256,
        "native_observations": 389,
        "seasons": seasons,
        "season_2023_definition": "All 74 observations in Xi Table 1, including December 2022",
        "lu_campaign_boundaries": CAMPAIGNS,
        "units": {"f5100": "1e-15 erg/s/cm^2/Angstrom",
                  "hbeta_total": "1e-13 erg/s/cm^2"},
        "interpolation": "none",
        "available_series": ["f5100", "hbeta_total"],
        "missing_required_series": ["hbeta_blue", "hbeta_core", "hbeta_red"],
        "status": "incomplete for requested five variable SURD analysis",
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
