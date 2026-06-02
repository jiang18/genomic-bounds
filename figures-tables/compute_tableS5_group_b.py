"""
compute_tableS5_group_b.py
--------------------------
Generates Supplementary Table S5: Group B properties underlying
Figures S1 and S2.

For each (Ne, h2, N) combination, reports:
  - S_A(N), S_B(N): Group A / Group B sigmoid sums
  - S_B / S_A: the small-fraction ratio plotted in Figure S1
  - U_B(N): Jensen upper bound on S_B
  - gap = U_B - S_B and rel_gap = (U_B - S_B) / U_B: the Jensen
    tightness diagnostic plotted in Figure S2

The single table supports both supplementary figures since they share
the same (Ne, h2, N) evaluation grid.

Usage:
    python compute_tableS5_group_b.py \
        --eigen_dir ../eigen_out --out_dir tables
"""

import argparse
import os
import numpy as np
import pandas as pd

from style import (NE_VALUES, me_theory, compute_S_decomposed)

H2_VALUES = [0.1, 0.3, 0.5]


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


def main(eigen_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    files_by_ne = discover_files(eigen_dir)

    rows = []
    for ne in NE_VALUES:
        if ne not in files_by_ne:
            continue
        largest_N = max(files_by_ne[ne].keys())
        M = len(np.load(files_by_ne[ne][largest_N]))
        Me = min(me_theory(ne), M)

        for n_val in sorted(files_by_ne[ne].keys()):
            ev = np.load(files_by_ne[ne][n_val])
            ev = np.maximum(ev, 0.0)
            if len(ev) < M:
                ev = np.concatenate([ev, np.zeros(M - len(ev))])

            for h2 in H2_VALUES:
                _, S_A, S_B, _, U_B = compute_S_decomposed(ev, n_val, h2, M, Me)
                ratio   = float(S_B / S_A) if S_A > 0 else np.nan
                gap     = float(U_B - S_B)
                rel_gap = (gap / U_B) if U_B > 0 else np.nan
                rows.append({
                    'Ne':         ne,
                    'M':          M,
                    'Me':        Me,
                    'h2':         h2,
                    'N':          n_val,
                    'S_A':        round(float(S_A), 3),
                    'S_B':        round(float(S_B), 3),
                    'SB_over_SA': round(ratio, 5),
                    'U_B':        round(float(U_B), 3),
                    'gap':        round(gap, 3),
                    'rel_gap':    round(rel_gap, 4),
                })

    df = pd.DataFrame(rows).sort_values(['Ne', 'h2', 'N']).reset_index(drop=True)
    out_path = os.path.join(out_dir, 'tableS5_group_b.csv')
    df.to_csv(out_path, index=False)

    print('\n=== Table S5: Group B properties underlying Figures S1 and S2 ===')
    print(df.to_string(index=False))
    print(f'\nMax S_B/S_A ratio:                {df["SB_over_SA"].max():.4f}')
    print(f'Max relative Jensen gap (U_B-S_B)/U_B: {df["rel_gap"].max():.4f}')
    print('\nNote: rows with S_B = U_B = 0 occur when N <= Me (no '
          'eigenvalues below the Me cutoff in that regime).')
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Table S5: Group B properties underlying Figures S1 and S2.')
    parser.add_argument('--eigen_dir', required=True,
                        help='Path to eigen_out/ directory')
    parser.add_argument('--out_dir', default='tables',
                        help='Output directory for CSV files')
    args = parser.parse_args()
    main(args.eigen_dir, args.out_dir)
