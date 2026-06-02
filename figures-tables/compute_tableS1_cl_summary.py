"""
compute_tableS1_cl_summary.py
-----------------------------
Generates Supplementary Table S1: detailed per-SNP GRAMMAR-Gamma coefficient
summary statistics for the three livestock chip-panel datasets (full-GRM
and LOCO).  Underlies Figure 1.

For each (dataset, mode, h2) combination, reports:
    N, M, Mean, SD, CV, Min, P10, Median, P90, Max

Heritability values: h2 = 0.1, 0.3, 0.5 (matches the paper's figures).

Inputs:
    Per-SNP coefficient TSV files in
        ../grammar-gamma-coeffs/analysis_outputs/grammar_gamma/results/
    {dataset}_coeffs.tsv          (full-GRM, columns c_h2_<h2>)
    {dataset}_loco_coeffs.tsv     (LOCO,     columns c_loco_h2_<h2>)

Output:
    tables/tableS1_cl_summary.csv

Usage:
    python compute_tableS1_cl_summary.py \
        --coeff_dir ../grammar-gamma-coeffs/analysis_outputs/grammar_gamma/results \
        --out_dir tables
"""

import argparse
import os
import numpy as np
import pandas as pd

H2_VALUES = [0.1, 0.3, 0.5]

DATASETS = [
    {
        'key':    'chinese_holstein',
        'label':  'Chinese Holstein',
        'full':   'chinese_holstein_coeffs.tsv',
        'loco':   'chinese_holstein_loco_coeffs.tsv',
    },
    {
        'key':    'karacabey_merino',
        'label':  'Karacabey Merino',
        'full':   'karacabey_merino_coeffs.tsv',
        'loco':   'karacabey_merino_loco_coeffs.tsv',
    },
    {
        'key':    'german_holstein',
        'label':  'German Holstein',
        'full':   'german_holstein_coeffs.tsv',
        'loco':   None,
    },
]

# Sample sizes (post-QC) — taken from the GRAMMAR-Gamma analysis
DATASET_N = {
    'chinese_holstein': 2510,
    'karacabey_merino': 734,
    'german_holstein':  5024,
}


def summarise(values):
    """Return mean, SD, CV, min, P10, median, P90, max for a vector."""
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    mean = float(v.mean())
    sd   = float(v.std(ddof=1))
    cv   = sd / mean if mean != 0 else np.nan
    return {
        'N_SNPs':   int(v.size),
        'Mean':     round(mean, 4),
        'SD':       round(sd, 4),
        'CV':       round(cv, 4),
        'Min':      round(float(v.min()), 4),
        'P10':      round(float(np.percentile(v, 10)), 4),
        'Median':   round(float(np.median(v)), 4),
        'P90':      round(float(np.percentile(v, 90)), 4),
        'Max':      round(float(v.max()), 4),
    }


def main(coeff_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    for ds in DATASETS:
        # full-GRM
        full_path = os.path.join(coeff_dir, ds['full'])
        full_df   = pd.read_csv(full_path, sep='\t')
        for h2 in H2_VALUES:
            col = f'c_h2_{h2}'
            stats = summarise(full_df[col].values)
            row = {
                'Dataset':  ds['label'],
                'Mode':     'full-GRM',
                'N':        DATASET_N[ds['key']],
                'h2':       h2,
            }
            row.update(stats)
            rows.append(row)

        # LOCO (if available)
        if ds['loco']:
            loco_path = os.path.join(coeff_dir, ds['loco'])
            loco_df   = pd.read_csv(loco_path, sep='\t')
            for h2 in H2_VALUES:
                col = f'c_loco_h2_{h2}'
                stats = summarise(loco_df[col].values)
                row = {
                    'Dataset':  ds['label'],
                    'Mode':     'LOCO',
                    'N':        DATASET_N[ds['key']],
                    'h2':       h2,
                }
                row.update(stats)
                rows.append(row)

    df = pd.DataFrame(rows)

    out_path = os.path.join(out_dir, 'tableS1_cl_summary.csv')
    df.to_csv(out_path, index=False)

    print('\n=== Table S1: per-SNP c_l summary (full-GRM and LOCO) ===')
    print(df.to_string(index=False))
    print(f'\nSaved: {out_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Table S1: per-SNP GRAMMAR-Gamma coefficient summary.')
    parser.add_argument(
        '--coeff_dir',
        default='../grammar-gamma-coeffs/analysis_outputs/grammar_gamma/results',
        help='Directory of per-SNP coefficient TSV files')
    parser.add_argument(
        '--out_dir', default='tables',
        help='Output directory for CSV files')
    args = parser.parse_args()
    main(args.coeff_dir, args.out_dir)
