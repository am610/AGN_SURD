"""
SURD AGN starter script
=======================

Purpose
-------
This script collects the main steps testing
the SURD repository for an AGN-style reverberation mapping toy problem.

What this script does
---------------------
1. Assumes the SURD repository has already been cloned in Colab or locally.
2. Adds the SURD `utils` directory to Python's import path.
3. Imports the core SURD modules.
4. Builds a synthetic AGN-like time-series data set:
      - hidden UV driver
      - optical continuum proxy
      - X-ray continuum proxy
      - delayed Hβ emission-line response
5. Runs the SURD decomposition for a single lag.
6. Scans across many lags and stores metrics for the final target signal.
7. Provides a helper that returns results for *all* target signals.
8. Shows how to run one of the built-in analytic benchmark systems.

Important notes
---------------
- This script intentionally avoids `transport_map.py`, because that path
  requires the optional `mpart` / MParT dependency.
- The SURD `run()` helper prints and plots results for all targets, but it
  returns only the decomposition for the *last* target processed. For clean
  analysis, this script includes a custom `run_collect()` wrapper.
- The lag parameter `nlag` is in samples, not physical days.

Expected environment
--------------------
This script is written for the same setup that worked in Colab:
    pip install pymp-pypi numpy scipy matplotlib pandas scikit-learn

and then:
    sys.path.append("/content/SURD/utils")

Adapt the SURD_UTILS_PATH below if your cloned repo lives elsewhere.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Any

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# 1. Configure path to SURD repo utilities
# =============================================================================
# Change this if your SURD clone is somewhere else.
SURD_UTILS_PATH = Path("/content/SURD/utils")

if not SURD_UTILS_PATH.exists():
    raise FileNotFoundError(
        f"Could not find SURD utils directory at {SURD_UTILS_PATH}. "
        "Update SURD_UTILS_PATH to match your local or Colab clone."
    )

sys.path.append(str(SURD_UTILS_PATH))

# Core SURD imports
import surd            # noqa: E402
import it_tools        # noqa: E402  # imported for completeness / inspection
import analytic_eqs    # noqa: E402


# =============================================================================
# 2. Utility helpers
# =============================================================================
def zscore(x: np.ndarray) -> np.ndarray:
    """
    Standardise an array to zero mean and unit variance.

    Parameters
    ----------
    x : np.ndarray
        Input array.

    Returns
    -------
    np.ndarray
        Standardised array.
    """
    x = np.asarray(x, dtype=float)
    return (x - np.mean(x)) / np.std(x)


def run_collect(
    X: np.ndarray,
    nvars: int,
    nlag: int,
    nbins: int,
    axs: np.ndarray | None = None,
    print_results: bool = False,
) -> Dict[int, Dict[str, Any]]:
    """
    Run SURD for each target signal and return a clean dictionary of results.

    Why this helper exists
    ----------------------
    The repo's `surd.run(...)` function loops over all targets and plots all of
    them, but returns only the decomposition of the final target in the loop.
    For analysis, it is much more useful to collect results for every target.

    Parameters
    ----------
    X : np.ndarray
        Array of shape (nvars, ntimes). Each row is one signal.
    nvars : int
        Number of signals.
    nlag : int
        Lag in samples.
    nbins : int
        Number of histogram bins used by np.histogramdd.
    axs : np.ndarray or None
        Optional axes array of shape (nvars, 2) for SURD plotting.
    print_results : bool
        If True, print the SURD decomposition for each target.

    Returns
    -------
    dict
        Dictionary keyed by target index, with values containing:
            - "I_R"
            - "I_S"
            - "MI"
            - "info_leak"
    """
    results = {}

    for i in range(nvars):
        # Future of target i versus present/past of all variables
        Y = np.vstack([X[i, nlag:], X[:, :-nlag]])

        # Multi-dimensional histogram approximation to the joint distribution
        hist, _ = np.histogramdd(Y.T, nbins)

        # SURD decomposition
        I_R, I_S, MI, info_leak = surd.surd(hist)

        results[i] = {
            "I_R": I_R,
            "I_S": I_S,
            "MI": MI,
            "info_leak": info_leak,
        }

        if print_results:
            print(f"\nSURD CAUSALITY FOR TARGET SIGNAL {i + 1}")
            surd.nice_print(I_R, I_S, MI, info_leak)

        if axs is not None:
            surd.plot(I_R, I_S, info_leak, axs[i, :], nvars, threshold=-0.01)

    return results


# =============================================================================
# 3. Synthetic AGN-style toy data
# =============================================================================
def make_synthetic_agn(
    N: int = 3000,
    tau: int = 8,
    seed: int = 42,
    random_walk_driver: bool = True,
) -> np.ndarray:
    """
    Build a synthetic AGN-like system with three signals.

    Signal definitions
    ------------------
    1. optical continuum: noisy proxy for hidden UV driver
    2. X-ray continuum:  noisy proxy for hidden UV driver
    3. Hβ emission line: delayed response to the UV driver

    Parameters
    ----------
    N : int
        Number of time samples.
    tau : int
        Reverberation lag in samples.
    seed : int
        Random seed for reproducibility.
    random_walk_driver : bool
        If True, use a random-walk hidden driver (long memory).
        If False, use white noise hidden driver (sharper lag localisation).

    Returns
    -------
    np.ndarray
        Stacked array X of shape (3, N).
    """
    rng = np.random.default_rng(seed)

    # Hidden UV driver
    if random_walk_driver:
        uv = np.cumsum(rng.normal(0, 0.15, N))
    else:
        uv = rng.normal(0, 1.0, N)
    uv = zscore(uv)

    # Two observed continuum proxies
    optical = 0.9 * uv + 0.25 * rng.normal(size=N)
    xray = 0.7 * uv + 0.35 * rng.normal(size=N)

    # Delayed broad-line response
    hbeta = np.zeros(N)
    hbeta[tau:] = 0.8 * uv[:-tau] + 0.25 * rng.normal(size=N - tau)
    hbeta[:tau] = 0.25 * rng.normal(size=tau)

    optical = zscore(optical)
    xray = zscore(xray)
    hbeta = zscore(hbeta)

    return np.vstack([optical, xray, hbeta])


# =============================================================================
# 4. Demonstration: single-lag SURD run
# =============================================================================
def demo_single_lag() -> None:
    """
    Run a single-lag SURD decomposition on the synthetic AGN toy system.
    """
    tau = 8
    X = make_synthetic_agn(N=3000, tau=tau, seed=42, random_walk_driver=True)
    nvars = X.shape[0]

    fig, axs = plt.subplots(nvars, 2, figsize=(12, 4 * nvars), squeeze=False)

    # This is the repo's built-in wrapper.
    # It prints and plots all targets, but returns only the final target result.
    surd.run(X=X, nvars=nvars, nlag=tau, nbins=8, axs=axs)

    plt.tight_layout()
    plt.show()


# =============================================================================
# 5. Demonstration: lag scan for target Signal 3
# =============================================================================
def demo_lag_scan() -> None:
    """
    Scan over lag values and track selected information measures.

    We focus on the final target (Signal 3, Hβ) by using run_collect().
    """
    tau_true = 8
    X = make_synthetic_agn(N=3000, tau=tau_true, seed=42, random_walk_driver=True)
    nvars = X.shape[0]

    lags = np.arange(1, 21)
    mi1, mi2, syn12, leak = [], [], [], []

    for lag in lags:
        results = run_collect(X=X, nvars=nvars, nlag=lag, nbins=8)

        # Target index 2 = Signal 3 in 1-based human counting
        res = results[2]

        # Notes:
        # MI keys refer to variable combinations in the SURD histogram indexing.
        # We keep the same convention used during the Colab exploration.
        mi1.append(res["MI"].get((1,), np.nan))
        mi2.append(res["MI"].get((2,), np.nan))
        syn12.append(res["I_S"].get((1, 2), np.nan))
        leak.append(res["info_leak"])

    plt.figure(figsize=(8, 5))
    plt.plot(lags, mi1, marker="o", label="MI(1)")
    plt.plot(lags, mi2, marker="o", label="MI(2)")
    plt.plot(lags, syn12, marker="o", label="Synergy(1,2)")
    plt.axvline(tau_true, linestyle="--", label="True lag")
    plt.xlabel("Lag (samples)")
    plt.ylabel("Information")
    plt.title("Lag scan for target Signal 3 (Hβ)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(7, 4))
    plt.plot(lags, leak, marker="o")
    plt.axvline(tau_true, linestyle="--")
    plt.xlabel("Lag (samples)")
    plt.ylabel("Information leak")
    plt.title("Leak vs lag for target Signal 3 (Hβ)")
    plt.tight_layout()
    plt.show()

    print("\nInterpretation note:")
    print(
        "If you use a random-walk hidden driver, the information curves often "
        "look broad rather than sharply peaked, because the driver has strong "
        "temporal autocorrelation."
    )
    print(
        "For a cleaner lag-localisation test, rerun with "
        "random_walk_driver=False in make_synthetic_agn()."
    )


# =============================================================================
# 6. Demonstration: compare random-walk vs white-noise hidden driver
# =============================================================================
def demo_driver_comparison() -> None:
    """
    Compare how lag localisation changes when the hidden driver has long memory
    versus no memory.
    """
    tau_true = 8
    lags = np.arange(1, 21)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4), squeeze=False)

    for ax, use_rw, title in zip(
        axs[0],
        [True, False],
        ["Random-walk driver", "White-noise driver"],
        strict=True,
    ):
        X = make_synthetic_agn(
            N=3000,
            tau=tau_true,
            seed=42,
            random_walk_driver=use_rw,
        )

        leak = []
        for lag in lags:
            results = run_collect(X=X, nvars=3, nlag=lag, nbins=8)
            leak.append(results[2]["info_leak"])

        ax.plot(lags, leak, marker="o")
        ax.axvline(tau_true, linestyle="--")
        ax.set_xlabel("Lag (samples)")
        ax.set_ylabel("Information leak")
        ax.set_title(title)

    plt.tight_layout()
    plt.show()


# =============================================================================
# 7. Demonstration: built-in analytic benchmark
# =============================================================================
def demo_analytic_benchmark() -> None:
    """
    Run one of the built-in synthetic benchmark systems from analytic_eqs.

    Available functions discovered during the Colab session:
        - confounder
        - mediator
        - redundant_collider
        - synergistic_collider
    """
    q1, q2, q3 = analytic_eqs.confounder(3000)
    X = np.vstack([q1, q2, q3])

    fig, axs = plt.subplots(3, 2, figsize=(10, 12), squeeze=False)
    surd.run(X=X, nvars=3, nlag=1, nbins=8, axs=axs)
    plt.tight_layout()
    plt.show()


# =============================================================================
# 8. Main execution block
# =============================================================================
if __name__ == "__main__":
    print("Core SURD AGN starter script")
    print("-" * 40)

    # Uncomment whichever demos you want to run.

    # 1) Single-lag toy example
    demo_single_lag()

    # 2) Lag scan for the Hβ target
    demo_lag_scan()

    # 3) Compare long-memory and short-memory hidden drivers
    # demo_driver_comparison()

    # 4) Run a built-in benchmark motif
    # demo_analytic_benchmark()

    print("\nDone.")
