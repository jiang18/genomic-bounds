"""
compute_table2_thresholds.py
----------------------------
Generates Table 2 of the manuscript: detection and fine-mapping thresholds
for livestock and human GWAS.

The output is *transposed* relative to a row-per-species layout: parameters
are rows and species are columns, so the table reads vertically by species.

All values are computed from closed-form formulas (no external data needed).
Species parameters follow Supplementary Methods S14:
    Cattle  (Ne=100,    L=25 Morgans, Me=10,000)
    Pig     (Ne=50,     L=20 Morgans, Me= 4,000)
    Chicken (Ne=50,     L=30 Morgans, Me= 6,000)
    Human   (Ne=10,000, L=35 Morgans, Me=1,400,000)

Output:
    tables/table2_thresholds.csv

Usage:
    python compute_table2_thresholds.py [--out_dir tables]
"""

import argparse
import os
import pandas as pd

# ---------------------------------------------------------------------------
# Species parameters (Supplementary Methods S14)
# ---------------------------------------------------------------------------
SPECIES = [
    {'name': 'Cattle',  'Ne': 100,   'L': 25},
    {'name': 'Pig',     'Ne': 50,    'L': 20},
    {'name': 'Chicken', 'Ne': 50,    'L': 30},
    {'name': 'Human',   'Ne': 10000, 'L': 35},
]

H2_VALUES = [0.1, 0.3, 0.5]

# Proximal-contamination mass fraction (>0.998 for all simulated populations)
RHO = 1.0

# Fine-mapping parameters
FM_N = 100_000
FM_THETA = 3
FM_DISTANCES_KB = [10, 100]    # 1 kb omitted: full-GRM ceiling at 1 kb is
                                # biologically infeasible in livestock (e.g.,
                                # 22.5% PVE for pig at h2=0.3) and 1 kb is
                                # below the resolution of practical interest.
FM_DISTANCES_M  = [d * 1e-5 for d in FM_DISTANCES_KB]


# ---------------------------------------------------------------------------
# Formulas
# ---------------------------------------------------------------------------
def meff(Ne, L):
    return 4 * Ne * L


def q_min(h2, Me):
    return 30.0 * h2 / Me


def N_A_star(Me, h2, rho):
    return Me * (1.0 - h2) / (h2 * rho)


def fm_min_pve(Ne, d_morgans, N, theta):
    """Simple regression / LOCO fine-mapping threshold.

    Equation 13 in main text: beta^2_min = 2 theta / (N (1 - r^2))
    with 1 - r^2 = 4 Ne d / (1 + 4 Ne d).
    """
    x = 4.0 * Ne * d_morgans
    one_minus_r2 = x / (1.0 + x)
    if one_minus_r2 == 0:
        return float('inf')
    return 2.0 * theta / (N * one_minus_r2)


def fm_min_pve_full_ceiling(Ne, L, h2, d_morgans, theta):
    """Full-GRM fine-mapping ceiling (Section S11.3).

    beta^2_min = 2 theta h^2 / (Me (1 - r^2))
              = theta h^2 (1 + 4 Ne d) / (8 Ne^2 L d)
    Independent of N (asymptotic limit as N -> infinity).
    """
    Me = meff(Ne, L)
    x = 4.0 * Ne * d_morgans
    one_minus_r2 = x / (1.0 + x)
    if one_minus_r2 == 0 or Me == 0:
        return float('inf')
    return 2.0 * theta * h2 / (Me * one_minus_r2)


# ---------------------------------------------------------------------------
# Pretty formatters
# ---------------------------------------------------------------------------
def fmt_int(x):
    return f'{int(round(x)):,}'


def fmt_pct(x):
    """Format a percentage value, switching to scientific notation when small."""
    if x >= 0.01:
        return f'{x:.3f}'
    elif x >= 1e-4:
        return f'{x:.4f}'
    else:
        return f'{x:.2e}'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # Compute everything we need per species
    rows_data = {}
    for sp in SPECIES:
        Ne, L = sp['Ne'], sp['L']
        Me = meff(Ne, L)
        d = {
            'Ne':       Ne,
            'L':        L,
            'Me':    Me,
        }
        for h2 in H2_VALUES:
            d[f'q_h2_{h2}']  = q_min(h2, Me) * 100.0          # %
            d[f'NA_h2_{h2}'] = N_A_star(Me, h2, RHO)
        # Fine-mapping under simple regression / LOCO (no h2 dependence)
        for d_kb, d_m in zip(FM_DISTANCES_KB, FM_DISTANCES_M):
            d[f'FM_{d_kb}kb'] = fm_min_pve(Ne, d_m, FM_N, FM_THETA) * 100.0  # %
        # Fine-mapping ceiling under full-GRM (at three h2 values, independent of N)
        for h2 in H2_VALUES:
            for d_kb, d_m in zip(FM_DISTANCES_KB, FM_DISTANCES_M):
                d[f'FM_full_h2_{h2}_{d_kb}kb'] = fm_min_pve_full_ceiling(
                    Ne, L, h2, d_m, FM_THETA) * 100.0
        rows_data[sp['name']] = d

    # ------------------------------------------------------------------
    # Build the transposed table (rows = parameters, cols = species)
    # ------------------------------------------------------------------
    species_names = [sp['name'] for sp in SPECIES]

    table_rows = []

    def add_row(label, values, fmt=str):
        row = {'Parameter': label}
        for sp_name, val in zip(species_names, values):
            row[sp_name] = fmt(val) if fmt else val
        table_rows.append(row)

    # Demographics
    add_row('Ne',
            [rows_data[s]['Ne'] for s in species_names], fmt_int)
    add_row('L (Morgans)',
            [rows_data[s]['L'] for s in species_names], fmt_int)
    add_row('Me',
            [rows_data[s]['Me'] for s in species_names], fmt_int)

    # Detection thresholds q_min (%)
    for h2 in H2_VALUES:
        add_row(f'q_min (%) at h2={h2}',
                [rows_data[s][f'q_h2_{h2}'] for s in species_names], fmt_pct)

    # Transition scale N_A*
    for h2 in H2_VALUES:
        add_row(f'N_A* at h2={h2}',
                [rows_data[s][f'NA_h2_{h2}'] for s in species_names], fmt_int)

    # Fine-mapping under simple regression at N=100k (no h2 dependence)
    for d_kb in FM_DISTANCES_KB:
        add_row(f'Min PVE (%), simple regression at d={d_kb} kb',
                [rows_data[s][f'FM_{d_kb}kb'] for s in species_names], fmt_pct)

    # Fine-mapping ceiling under full-GRM at three h2 values (independent of N)
    for h2 in H2_VALUES:
        for d_kb in FM_DISTANCES_KB:
            add_row(f'Min PVE (%), full-GRM ceiling at d={d_kb} kb (h2={h2})',
                    [rows_data[s][f'FM_full_h2_{h2}_{d_kb}kb']
                     for s in species_names], fmt_pct)

    df = pd.DataFrame(table_rows)

    # Save
    out_path = os.path.join(out_dir, 'table2_thresholds.csv')
    df.to_csv(out_path, index=False)

    print('\n=== Table 2: Detection and Fine-Mapping Thresholds (transposed) ===\n')
    print(df.to_string(index=False))
    print(f'\nSimple-regression fine-mapping thresholds use '
          f'N = {FM_N:,}, theta = {FM_THETA}.')
    print(f'Full-GRM fine-mapping ceiling is independent of N '
          f'(asymptotic limit at theta = {FM_THETA}).')
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute Table 2: detection and fine-mapping thresholds.')
    parser.add_argument(
        '--out_dir', default='tables',
        help='Output directory for CSV files (default: tables)')
    args = parser.parse_args()
    main(args.out_dir)
