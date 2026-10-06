"""Expose all atoms from the unchanged authors' SURD implementation."""

import hashlib
import itertools
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
UTILS = ROOT / 'SURD' / 'utils'
HASHES = {
    'surd.py': 'c843b337d6e1fc4437d26ee1551980c21da973cf0fb1f3f5f1f82455773bd2bf',
    'it_tools.py': '3a6c752e37961774cc951dbaf85d15e2fc39486c3d42efa8e07d0ede90554e18',
}


def official_module():
    for name, expected in HASHES.items():
        if hashlib.sha256((UTILS / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Official source hash mismatch: {name}')
    sys.path.insert(0, str(UTILS))
    import surd
    if Path(surd.__file__).resolve() != (UTILS / 'surd.py').resolve():
        raise ValueError('Unexpected SURD import location')
    return surd


def decompose_counts(counts, predictor_names):
    """Return every individual atom and normalized leakage for a joint histogram.

    Axis zero is the future target. Subsequent axes are past predictors,
    including the target history when specified by the caller. Atom fractions
    use joint mutual information; leakage uses future target entropy.
    """
    counts = np.asarray(counts, dtype=float)
    names = tuple(predictor_names)
    if counts.ndim != len(names) + 1 or len(set(names)) != len(names):
        raise ValueError('Histogram axes must match distinct predictor labels')
    if not np.isfinite(counts).all() or (counts < 0).any() or counts.sum() <= 0:
        raise ValueError('Invalid histogram counts')
    if np.count_nonzero(counts.sum(axis=tuple(range(1, counts.ndim)))) < 2:
        raise ValueError('Target must have at least two observed states')
    redundancy, synergy, mi, leakage = official_module().surd(counts.copy())
    joint = float(mi[tuple(range(1, counts.ndim))])
    atoms = []
    for order in range(1, len(names) + 1):
        for indices in itertools.combinations(range(1, len(names) + 1), order):
            labels = [names[i - 1] for i in indices]
            kinds = [('U' if order == 1 else 'R', redundancy)]
            if order > 1:
                kinds.append(('S', synergy))
            for kind, values in kinds:
                bits = float(values[indices])
                atoms.append(dict(kind=kind, predictors=labels, bits=bits,
                                  fraction=bits / joint if joint > 1e-10 else None))
    if not np.isclose(sum(a['bits'] for a in atoms), joint, atol=1e-8, rtol=1e-8):
        raise ValueError('Atom sum differs from joint mutual information')
    return dict(atoms=atoms, joint_mi_bits=joint,
                normalized_leakage=float(leakage),
                normalization_defined=joint > 1e-10)
