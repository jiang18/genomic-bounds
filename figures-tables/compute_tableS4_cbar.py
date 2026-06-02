"""
compute_tableS4_cbar.py
-----------------------
Generates Supplementary Table S4: average GRAMMAR-Gamma coefficients
(c-bar full-GRM and LOCO) underlying Figure 4.

Reads the three pre-computed cbar CSV tables in ../msprime/ and
stacks them into one publication-ready table.

Usage:
    python compute_tableS4_cbar.py --cbar_dir ../msprime --out_dir tables
"""

import argparse
import os
import pandas as pd

NE_VALUES = [20, 50, 100]


def main(cbar_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    frames = []
    for ne in NE_VALUES:
        path = os.path.join(cbar_dir, f'cbar_Ne{ne}.csv')
        if not os.path.isfile(path):
            print(f'WARNING: {path} not found, skipping Ne={ne}')
            continue
        df = pd.read_csv(path)
        df.insert(0, 'Ne', ne)
        frames.append(df)

    if not frames:
        raise FileNotFoundError(f'No cbar CSV files found in {cbar_dir}')

    df = pd.concat(frames, ignore_index=True)
    # Reorder/relabel columns for the publication table
    keep_cols = ['Ne', 'N', 'M', 'M_loco', 'h2', 'h2_loco',
                 'cbar_full', 'cbar_loco', 'loco_upper_bound']
    keep_cols = [c for c in keep_cols if c in df.columns]
    df = df[keep_cols]

    # Round numeric columns for readability
    for c in df.columns:
        if df[c].dtype.kind == 'f':
            df[c] = df[c].round(4)

    df = df.sort_values(['Ne', 'h2', 'N']).reset_index(drop=True)

    out_path = os.path.join(out_dir, 'tableS4_cbar.csv')
    df.to_csv(out_path, index=False)
    print('\n=== Table S4: c-bar (full-GRM and LOCO) underlying Figure 4 ===')
    print(df.to_string(index=False))
    print(f'\nSaved: {out_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Table S4: c-bar values underlying Figure 4.')
    parser.add_argument('--cbar_dir', default='../msprime',
                        help='Directory containing cbar_Ne{20,50,100}.csv')
    parser.add_argument('--out_dir', default='tables',
                        help='Output directory for CSV files')
    args = parser.parse_args()
    main(args.cbar_dir, args.out_dir)
