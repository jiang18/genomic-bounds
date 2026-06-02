"""
compute_tableS3_sigmoid_sum.py
------------------------------
Generates Supplementary Table S3: S(N) values underlying Figure 3.

The table has two blocks distinguished by the column ``fig3_row``:

  fig3_row = 'top'  (panels A, B, C of Figure 3)
      h2 in {0.1, 0.3, 0.5} -- the three curves plotted in the top row.
      For each (Ne, h2), evaluate S(N) at every available N in the
      simulation grid (the full N-trajectory shown in the top-row
      curves).

  fig3_row = 'bottom'  (panels D, E, F of Figure 3)
      N in {10,000; 20,000; 50,000} -- the three curves plotted in
      the bottom row.  For each (Ne, N), evaluate S(N) at a dense
      grid of h2 values spanning the plotted range [0.05, 0.95].
      This is the dense h2 sampling needed to verify the shape of
      the bottom-row curves.

For each row, the table reports:
    fig3_row, Ne, M, Me, h2, N, Mlambda, S, S_A, S_B, S_over_Me

Usage:
    python compute_tableS3_sigmoid_sum.py \
        --eigen_dir ../eigen_out --out_dir tables
"""

import argparse
import os
import numpy as np
import pandas as pd

from style import (NE_VALUES, me_theory, compute_S_decomposed)

# Top-row: three heritability curves plotted at all simulated N.
H2_TOP = [0.1, 0.3, 0.5]

# Bottom-row: three fixed N values plotted at a dense h2 grid.
N_BOTTOM = [10_000, 20_000, 50_000]
H2_BOTTOM = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50,
             0.60, 0.70, 0.80, 0.90, 0.95]


def discover_files(eigen_dir):
    data = {}
    for ne in NE_VALUES:
        ne_dir = os.path.join(eigen_dir, f'Ne{ne}')
        if not os.path.isdir(ne_dir):
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


def _load_eigs(path, M):
    ev = np.load(path)
    ev = np.maximum(ev, 0.0)
    if len(ev) < M:
        ev = np.concatenate([ev, np.zeros(M - len(ev))])
    return ev


def _row(fig3_row, ne, M, Me, h2, n_val, ev):
    M_lambda = M * (1.0 - h2) / h2
    S, S_A, S_B, _, _ = compute_S_decomposed(ev, n_val, h2, M, Me)
    return {
        'fig3_row':     fig3_row,
        'Ne':           ne,
        'M':            M,
        'Me':        Me,
        'h2':           h2,
        'N':            n_val,
        'Mlambda':      round(M_lambda, 2),
        'S':            round(float(S), 3),
        'S_A':          round(float(S_A), 3),
        'S_B':          round(float(S_B), 3),
        'S_over_Me':  round(float(S) / Me, 4),
    }


def main(eigen_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    files_by_ne = discover_files(eigen_dir)

    rows = []
    for ne in NE_VALUES:
        if ne not in files_by_ne:
            continue
        # Use the largest eigenvalue array to determine M.
        largest_N = max(files_by_ne[ne].keys())
        M = len(np.load(files_by_ne[ne][largest_N]))
        Me = min(me_theory(ne), M)

        # ---- Top row: full N grid x H2_TOP ----
        for n_val in sorted(files_by_ne[ne].keys()):
            ev = _load_eigs(files_by_ne[ne][n_val], M)
            for h2 in H2_TOP:
                rows.append(_row('top', ne, M, Me, h2, n_val, ev))

        # ---- Bottom row: N_BOTTOM x H2_BOTTOM ----
        for n_val in N_BOTTOM:
            if n_val not in files_by_ne[ne]:
                continue
            ev = _load_eigs(files_by_ne[ne][n_val], M)
            for h2 in H2_BOTTOM:
                rows.append(_row('bottom', ne, M, Me, h2, n_val, ev))

    df = pd.DataFrame(rows)
    # Sort: top block first (by Ne, h2, N), then bottom block (by Ne, N, h2).
    df_top = (df[df['fig3_row'] == 'top']
              .sort_values(['Ne', 'h2', 'N']))
    df_bot = (df[df['fig3_row'] == 'bottom']
              .sort_values(['Ne', 'N', 'h2']))
    df_out = pd.concat([df_top, df_bot], ignore_index=True)

    out_path = os.path.join(out_dir, 'tableS3_sigmoid_sum.csv')
    df_out.to_csv(out_path, index=False)
    print('\n=== Table S3: S(N) values underlying Figure 3 ===')
    print(f'Top-row rows (panels A, B, C):    {len(df_top):4d}')
    print(f'Bottom-row rows (panels D, E, F): {len(df_bot):4d}')
    print(f'Total rows:                       {len(df_out):4d}')
    print(f'\nSaved: {out_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Table S3: sigmoid sum S(N) underlying Figure 3.')
    parser.add_argument('--eigen_dir', required=True,
                        help='Path to eigen_out/ directory')
    parser.add_argument('--out_dir', default='tables',
                        help='Output directory for CSV files')
    args = parser.parse_args()
    main(args.eigen_dir, args.out_dir)
