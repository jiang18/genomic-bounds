"""
compute_table1_eigenspectrum.py
-------------------------------
Generates the eigenvalue spectral-parameter tables.

Outputs (all CSV, in tables/):

  table1_eigenspectrum.csv           [Main Table 1]
      Me = 4*Ne*L (theoretical).  Long form: one row per (Ne, h2)
      combination, so 9 rows.  Columns:
          Ne, M, h2, Me, rho, N_A_star, N_B_star, ratio_NB_over_NA
      Eigenvalues are computed at the largest available N for each Ne
      (N=100,000 for Ne=20 and 50; N=200,000 for Ne=100) so that the
      spectrum is well-resolved.

  tableS6_eig98.csv                  [Supp]
      Same long-form layout, Me = number of top eigenvalues capturing
      98% of total spectral mass.  Adds an 'EIG98_over_4NeL' column
      showing the ratio against the theoretical Me.

  tableS7_eig99.csv                  [Supp]
      Same as above with 99% mass cutoff.

Usage:
    python compute_table1_eigenspectrum.py \
        --eigen_dir ../eigen_out --out_dir tables

Eigenvalue files: eigen_out/Ne{Ne}/eigenvalues_N{N}.npy
Each .npy: 1-D float64 array, length = min(N, M), sorted descending,
           eigenvalues of the LD matrix R = Z'Z / N.
"""

import argparse
import os
import numpy as np
import pandas as pd

from style import (NE_VALUES, me_theory, compute_rho,
                   compute_NA_star, compute_NB_star)

H2_VALUES = [0.1, 0.3, 0.5]


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------
def discover_files(eigen_dir):
    data = {}
    for ne in NE_VALUES:
        ne_dir = os.path.join(eigen_dir, f'Ne{ne}')
        if not os.path.isdir(ne_dir):
            print(f'WARNING: directory not found: {ne_dir}')
            continue
        data[ne] = {}
        for fname in sorted(os.listdir(ne_dir)):
            if fname.startswith('eigenvalues_N') and fname.endswith('.npy'):
                n_str = fname.replace('eigenvalues_N', '').replace('.npy', '')
                try:
                    data[ne][int(n_str)] = os.path.join(ne_dir, fname)
                except ValueError:
                    pass
    return data


def infer_M(files_by_ne):
    M_by_ne = {}
    for ne, files in files_by_ne.items():
        largest_N = max(files.keys())
        eigs = np.load(files[largest_N])
        M_by_ne[ne] = len(eigs)
    return M_by_ne


def eig_threshold(eigenvalues, frac):
    total = eigenvalues.sum()
    if total == 0:
        return 0
    cumsum = np.cumsum(eigenvalues)
    idx = np.searchsorted(cumsum, frac * total)
    return int(idx) + 1


# ---------------------------------------------------------------------------
# Build long-form rows
# ---------------------------------------------------------------------------
def build_long_rows(ne, M, eigs, Me_value, extra_cols=None):
    """
    Returns 3 rows (one per h2) with Ne, M, h2, Me, rho, N_A*, N_B*, ratio.
    extra_cols: dict to merge into every row (e.g. ratio of EIG vs 4NeL).
    """
    Me = min(Me_value, len(eigs))
    rho  = compute_rho(eigs, Me)
    rows = []
    for h2 in H2_VALUES:
        NA = compute_NA_star(Me, h2, rho)
        NB = compute_NB_star(M, Me, h2, rho)
        ratio = NB / NA if (NA > 0 and np.isfinite(NB)) else np.inf
        row = {
            'Ne':             ne,
            'M':              M,
            'h2':             h2,
            'Me':          Me,
            'rho':            round(rho, 4),
            'N_A_star':       int(round(NA)),
            'N_B_star':       (int(round(NB)) if np.isfinite(NB) else 'inf'),
            'ratio_NB_NA':    (int(round(ratio)) if np.isfinite(ratio) else 'inf'),
        }
        if extra_cols:
            row.update(extra_cols)
        rows.append(row)
    return rows


def main(eigen_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    files_by_ne = discover_files(eigen_dir)
    M_by_ne     = infer_M(files_by_ne)

    rows_main = []
    rows_e98  = []
    rows_e99  = []

    for ne in NE_VALUES:
        if ne not in files_by_ne:
            continue
        M         = M_by_ne[ne]
        largest_N = max(files_by_ne[ne].keys())
        eigs      = np.load(files_by_ne[ne][largest_N])
        eigs      = np.maximum(eigs, 0.0)

        # Theory-driven Me = 4*Ne*L
        Me_th = min(me_theory(ne), M)
        rows_main.extend(
            build_long_rows(ne, M, eigs, Me_th)
        )

        # Data-driven: EIG_98%
        Me_98 = eig_threshold(eigs, 0.98)
        rows_e98.extend(
            build_long_rows(ne, M, eigs, Me_98,
                            extra_cols={'EIG98_over_4NeL': round(Me_98 / Me_th, 3)})
        )

        # Data-driven: EIG_99%
        Me_99 = eig_threshold(eigs, 0.99)
        rows_e99.extend(
            build_long_rows(ne, M, eigs, Me_99,
                            extra_cols={'EIG99_over_4NeL': round(Me_99 / Me_th, 3)})
        )

    for rows, fname, label in [
        (rows_main, 'table1_eigenspectrum.csv',
         'Main Table 1 (Me = 4*Ne*L, long form)'),
        (rows_e98,  'tableS6_eig98.csv',
         'Supp Table S6 (Me = EIG_98%, long form)'),
        (rows_e99,  'tableS7_eig99.csv',
         'Supp Table S7 (Me = EIG_99%, long form)'),
    ]:
        df   = pd.DataFrame(rows)
        path = os.path.join(out_dir, fname)
        df.to_csv(path, index=False)
        print(f'\n=== {label} ===')
        print(df.to_string(index=False))
        print(f'Saved: {path}')

    print('\n[Footnote] Eigenvalues computed at the largest simulated N '
          '(100,000 for Ne=20 and 50; 200,000 for Ne=100) so that the '
          'spectrum is well-resolved.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute Table 1 (long form) and supplementary '
                    'EIG98/EIG99 tables.')
    parser.add_argument('--eigen_dir', required=True,
                        help='Path to eigen_out/ directory')
    parser.add_argument('--out_dir', default='tables',
                        help='Output directory for CSV files (default: tables)')
    args = parser.parse_args()
    main(args.eigen_dir, args.out_dir)
