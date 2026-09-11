"""Run observed scans and a reproducible, resumable local simulation pilot."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import time

import numpy as np
import pandas as pd

from .adopted_data import ROOT, load_adopted
from .information import scan, peaks

HERE = Path(__file__).resolve().parent


def peak_memory_mib():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / (1024 ** 2 if platform.system() == 'Darwin' else 1024)


def ou(rng, times, tau):
    values = np.empty(len(times))
    values[0] = rng.normal()
    coefficients = np.exp(-np.diff(times) / tau)
    innovations = rng.normal(size=len(times) - 1) * np.sqrt(1 - coefficients ** 2)
    for i, (coefficient, innovation) in enumerate(zip(coefficients, innovations), 1):
        values[i] = coefficient * values[i - 1] + innovation
    return values


def simulate(lines, cont, config, case, seed):
    rng = np.random.default_rng(seed)
    start = min(lines.jd_offset.min(), cont.jd_offset.min()) - config['simulation_burn_days']
    stop = max(lines.jd_offset.max(), cont.jd_offset.max())
    step = config['simulation_grid_days']
    grid = np.arange(start, stop + step, step)
    driver = ou(rng, grid, config['driver_timescale_days'])
    sim_lines, sim_cont = lines.copy(), cont.copy()
    continuum_scale = float(cont.continuum.std(ddof=1))
    sim_cont['continuum'] = (cont.continuum.mean() + continuum_scale * np.interp(cont.jd_offset, grid, driver)
                             + rng.normal(size=len(cont)) * cont.continuum_error.to_numpy())
    response_lags = np.arange(config['response_delay_days'] - config['response_width_days'] / 2,
                              config['response_delay_days'] + config['response_width_days'] / 2 + step / 2, step)
    if response_lags.min() < 0 or config['simulation_burn_days'] < response_lags.max():
        raise ValueError('Simulation history does not cover response')
    response = np.mean([np.interp(lines.jd_offset.to_numpy() - delay, grid, driver) for delay in response_lags], axis=0)
    for name in ('blue', 'core', 'red'):
        independent = np.interp(lines.jd_offset, grid, ou(rng, grid, config['driver_timescale_days']))
        if case == 'independent':
            signal = independent
        elif case == 'delayed_shared_driver':
            signal = config['response_amplitude'] * response + config['independent_line_amplitude'] * independent
        else:
            raise ValueError('Unknown simulation case')
        sim_lines[name] = (lines[name].mean() + lines[name].std(ddof=1) * signal
                           + rng.normal(size=len(lines)) * lines[name + '_error'].to_numpy())
    sim_lines['profile_total'] = sim_lines[['blue', 'core', 'red']].sum(axis=1)
    return sim_lines, sim_cont


def worker(job):
    case, seed, sampling, config, output = job
    stem = Path(output) / f'{case}_{seed}'
    marker = stem.with_suffix('.json')
    if marker.exists():
        saved = json.loads(marker.read_text())
        for suffix, digest in saved['output_sha256'].items():
            if hashlib.sha256(Path(str(stem) + suffix).read_bytes()).hexdigest() != digest:
                raise ValueError('Checkpoint output hash mismatch')
        return saved
    start = time.perf_counter()
    lines, cont = load_adopted(sampling['continuum'])
    simulated = simulate(lines, cont, config, case, seed)
    curves = scan(*simulated, sampling, config)
    summary = peaks(curves)
    curves.to_csv(str(stem) + '_curves.csv', index=False)
    summary.to_csv(str(stem) + '_peaks.csv', index=False)
    result = dict(case=case, seed=seed, elapsed_seconds=time.perf_counter() - start,
                  process_peak_memory_mib=peak_memory_mib(),
                  output_sha256={s: hashlib.sha256(Path(str(stem) + s).read_bytes()).hexdigest()
                                 for s in ['_curves.csv', '_peaks.csv']})
    marker.write_text(json.dumps(result, indent=2) + '\n')
    print(f'Completed {case} seed {seed}: {result["elapsed_seconds"]:.2f} seconds', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('workers must be positive')
    sampling = json.loads((HERE / 'sampling_config.json').read_text())
    config = json.loads((HERE / 'information_config.json').read_text())
    sources = [HERE / n for n in ['sampling_config.json', 'information_config.json', 'sampling.py',
                                 'information.py', 'benchmark_information.py', 'adopted_data.py', 'reconcile_data.py']]
    sources += [ROOT / 'SURD/utils/surd.py', ROOT / 'SURD/utils/it_tools.py',
                ROOT / 'agn_surd_project/processed/reconciled/adopted/manifest.json']
    fingerprints = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    adopted = json.loads(sources[-1].read_text())
    for name, expected in adopted['sha256'].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Adopted input changed: {name}; regenerate adopted products')
    identity = hashlib.sha256(json.dumps(fingerprints, sort_keys=True).encode()).hexdigest()[:16]
    out = ROOT / 'agn_surd_project/processed/reconciled/information_pilot' / identity
    out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    lines, cont = load_adopted(sampling['continuum'])
    observed = scan(lines, cont, sampling, config)
    observed.to_csv(out / 'observed_curves.csv', index=False)
    peaks(observed).to_csv(out / 'observed_peaks.csv', index=False)
    observed_seconds = time.perf_counter() - start
    print(f'Observed scan: {len(observed)} rows in {observed_seconds:.2f} seconds', flush=True)
    jobs = [(case, seed, sampling, config, str(out)) for case in config['simulation_cases'] for seed in config['pilot_seeds']]
    cached = sum((out / f'{case}_{seed}.json').exists() for case, seed, *_ in jobs)
    start = time.perf_counter()
    if args.workers == 1:
        results = [worker(job) for job in jobs]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            results = list(executor.map(worker, jobs))
    wall = time.perf_counter() - start
    elapsed = np.array([r['elapsed_seconds'] for r in results])
    report = dict(status='Runtime and correctness pilot only; not calibrated significance or power',
                  sampling=sampling, information=config, sha256=fingerprints,
                  platform=platform.platform(), python=platform.python_version(),
                  numpy=np.__version__, pandas=pd.__version__, logical_cpus=os.cpu_count(),
                  workers=args.workers, observed_seconds=observed_seconds, simulation_wall_seconds=wall,
                  parent_peak_memory_mib=peak_memory_mib(),
                  maximum_worker_peak_memory_mib=max(r['process_peak_memory_mib'] for r in results),
                  simulations=len(results), cached_simulations=cached,
                  median_simulation_seconds=float(np.median(elapsed)),
                  throughput_seconds_per_simulation=wall / len(results) if cached == 0 else None,
                  projected_1000_simulations_hours=1000 * wall / len(results) / 3600 if cached == 0 else None,
                  projection_limit='Same cases, estimator, scan extent and worker count only; excludes nested calibration and expanded parameter grids',
                  simulation_limit='Noise added at aggregated native dates using propagated errors; shared calibration covariance omitted; process parameters not fitted',
                  results=results)
    report_name = 'benchmark.json' if cached == 0 else 'benchmark_resume.json'
    (out / report_name).write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['observed_seconds', 'simulation_wall_seconds', 'simulations',
                                            'cached_simulations', 'projected_1000_simulations_hours']}, indent=2))
    print(f'Outputs: {out}')


if __name__ == '__main__':
    main()
