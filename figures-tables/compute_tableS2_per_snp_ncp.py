"""
compute_tableS2_per_snp_ncp.py
------------------------------
Generates Supplementary Table S2: per-SNP NCP underlying Figure 2.

The table has 360 rows = 3 Ne x 3 h2 x 2 beta2 x 4 N x 5 focal SNPs.
Columns: Ne, h2, beta2, N, focal_SNP, chr, MAF, c_l,
         theoretical_NCP, empirical_NCP, ratio.

MAF is read from each population's focal_snps_info.csv (computed using all
individuals; one MAF per focal SNP, independent of subsample N).

Usage:
    python compute_tableS2_per_snp_ncp.py \
        --plink_dir ../plink_out --out_dir tables
"""

import argparse
import os
import re
import numpy as np
import pandas as pd

RE_NE   = re.compile(r'Ne(\d+)')
RE_SCEN = re.compile(r'h2_([\d.]+)__beta2_([\d.]+)')
RE_N    = re.compile(r'N_(\d+)')


def load_focal_info(plink_dir):
    """Return dict {Ne: DataFrame from focal_snps_info.csv}."""
    info = {}
    for ne_dir in sorted(os.listdir(plink_dir)):
        m_ne = RE_NE.match(ne_dir)
        if not m_ne:
            continue
        ne = int(m_ne.group(1))
        path = os.path.join(plink_dir, ne_dir, 'focal_snps_info.csv')
        if os.path.isfile(path):
            info[ne] = pd.read_csv(path)
    return info


def collect(plink_dir, focal_info):
    rows = []
    for ne_dir in sorted(os.listdir(plink_dir)):
        m_ne = RE_NE.match(ne_dir)
        if not m_ne:
            continue
        ne = int(m_ne.group(1))
        ne_path = os.path.join(plink_dir, ne_dir)
        if not os.path.isdir(ne_path):
            continue

        for scen_dir in sorted(os.listdir(ne_path)):
            m_sc = RE_SCEN.match(scen_dir)
            if not m_sc:
                continue
            h2    = float(m_sc.group(1))
            beta2 = float(m_sc.group(2))
            sc_path = os.path.join(ne_path, scen_dir)

            for n_dir in sorted(os.listdir(sc_path)):
                m_n = RE_N.match(n_dir)
                if not m_n:
                    continue
                N = int(m_n.group(1))
                csv_path = os.path.join(sc_path, n_dir, 'summary.csv')
                if not os.path.isfile(csv_path):
                    continue

                df = pd.read_csv(csv_path)

                # First column 'covar' is the focal SNP id
                for _, r in df.iterrows():
                    snp_id = r['covar']
                    # Look up MAF + chr from focal_snps_info
                    info_df = focal_info.get(ne)
                    chr_val = np.nan
                    maf_val = np.nan
                    if info_df is not None:
                        match = info_df[info_df['SNP'] == snp_id]
                        if not match.empty:
                            chr_val = int(match.iloc[0]['CHR'])
                            maf_val = float(match.iloc[0]['MAF'])

                    rows.append({
                        'Ne':              ne,
                        'h2':              h2,
                        'beta2':           beta2,
                        'N':               N,
                        'focal_SNP':       snp_id,
                        'chr':             chr_val,
                        'MAF':             maf_val,
                        'c_l':             round(float(r['c_l']),
                                                 6) if not pd.isna(r['c_l']) else np.nan,
                        'theoretical_NCP': round(float(r['ncp_theoretical_cl']),
                                                 4) if not pd.isna(r['ncp_theoretical_cl']) else np.nan,
                        'empirical_NCP':   round(float(r['ncp_empirical']),
                                                 4) if not pd.isna(r['ncp_empirical']) else np.nan,
                        'ratio':           round(float(r['ncp_ratio']),
                                                 4) if not pd.isna(r['ncp_ratio']) else np.nan,
                    })
    return pd.DataFrame(rows)


def main(plink_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    focal_info = load_focal_info(plink_dir)
    df = collect(plink_dir, focal_info)
    df = df.sort_values(['Ne', 'h2', 'beta2', 'N', 'chr',
                         'focal_SNP']).reset_index(drop=True)

    out_path = os.path.join(out_dir, 'tableS2_per_snp_ncp.csv')
    df.to_csv(out_path, index=False)

    print('\n=== Table S2: per-SNP NCP underlying Figure 2 ===')
    print(f'Total rows: {len(df)}  '
          f'(expected 3 Ne x 3 h2 x 2 beta2 x 4 N x 5 SNPs = 360)')
    print(df.head(8).to_string(index=False))
    print('...')
    print(f'\nSaved: {out_path}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Table S2: per-SNP NCP underlying Figure 2.')
    parser.add_argument('--plink_dir', required=True,
                        help='Path to plink_out/ directory')
    parser.add_argument('--out_dir', default='tables',
                        help='Output directory for CSV files')
    args = parser.parse_args()
    main(args.plink_dir, args.out_dir)
