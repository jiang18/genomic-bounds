#!/usr/bin/env python3
"""
collect_results.py — Summarize focal SNP test statistics across replicates.

Usage:
    python3 collect_results.py <simdir> <eigen_dir> <h2> <beta2>

Example:
    python3 collect_results.py ./plink_out/Ne20/h2_0.3__beta2_0.01/N_5000 ./eigen_out/Ne20 0.3 0.01

Reads rep_*/reml.reml.blue.csv, extracts BLUE/SE/pval for focal SNPs,
computes chi-squared = (BLUE/SE)^2, and reports E[chi2] - 1 (empirical NCP).

Reads focal_proj_N{N}.npz from eigen_dir (saved by project_focal_snps.py)
to compute per-focal-SNP c_l and theoretical NCP:
    c_l = (1/N) * sum_k a_{lk}^2 / (d_k * h2 + 1 - h2)
    NCP_l = (c_l * N * h2 / M) * Wtilde
where Wss = beta2*M/h2, Wtilde = Wss - 1.

Note: SLEMM reports BLUE relative to A1 allele (--recode A coding), while
sim.snp.csv uses A2 allele convention. We negate BLUE before unbiasedness
checks so that the sign convention matches the true effect.

Output:
    <simdir>/summary.csv
    <simdir>/per_replicate_chi2.csv
    <simdir>/blue_bias_check.csv
"""

import pandas as pd
import numpy as np
import sys, os, glob, re

if len(sys.argv) < 5:
    print("Usage: python3 collect_results.py <simdir> <eigen_dir> <h2> <beta2>")
    sys.exit(1)

simdir = sys.argv[1]
eigen_dir = sys.argv[2]
h2 = float(sys.argv[3])
beta2 = float(sys.argv[4])

# --- Parse N from simdir path ---
N_match = re.search(r'N_(\d+)', simdir)
if not N_match:
    print("ERROR: Cannot parse N from simdir path")
    sys.exit(1)
N = int(N_match.group(1))
print(f"  N={N}, h2={h2}, beta2={beta2}")

# --- Load precomputed projections and compute c_l ---
proj_file = os.path.join(eigen_dir, f"focal_proj_N{N}.npz")
cl_data = None
cbar = S_N = Wss = Wtilde = M = None

if os.path.exists(proj_file):
    proj = np.load(proj_file, allow_pickle=True)
    d = proj['d']            # GRM eigenvalues, length N
    a_matrix = proj['a']     # n_focal x N projections
    focal_proj_ids = list(proj['focal_ids'])
    M = int(proj['M'])

    # Compute c_l for this h2
    inv_w = 1.0 / (d * h2 + (1 - h2))
    cbar_N = np.sum(d * inv_w)
    cbar = cbar_N / N
    S_N = h2 * cbar_N
    alpha = M * (1 - h2) / h2
    Wss = beta2 * M / h2
    Wtilde = Wss - 1

    cl_data = {}
    for i, sid in enumerate(focal_proj_ids):
        a2 = a_matrix[i, :] ** 2
        cl = np.sum(a2 * inv_w) / N
        cl_data[sid] = {
            'c_l': cl,
            'c_l_N': cl * N,
            'ratio_cl_cbar': cl / cbar,
        }

    print(f"  Loaded projections for {len(focal_proj_ids)} focal SNPs at N={N}")
    print(f"  M={M}, h2={h2}, Wss={Wss:.2f}, Wtilde={Wtilde:.2f}")
    print(f"  c_bar={cbar:.6f}, S(N)={S_N:.2f}, alpha={alpha:.2f}")
else:
    print(f"  WARNING: {proj_file} not found — theoretical NCP will not be computed")

# Find all replicate BLUE files
blue_files = sorted(glob.glob(f"{simdir}/rep_*/reml.reml.blue.csv"))
nrep = len(blue_files)
if nrep == 0:
    print(f"ERROR: No reml.reml.blue.csv files found in {simdir}/rep_*/")
    sys.exit(1)

print(f"=== Collecting results from {simdir} ({nrep} replicates) ===")

# --- Load true effects from sim.snp.csv files ---
focal_snp_ids = set()
true_effects_per_rep = {}
for rep_idx in range(1, nrep + 1):
    snp_file = f"{simdir}/rep_{rep_idx}/sim.snp.csv"
    if os.path.exists(snp_file):
        snp_df = pd.read_csv(snp_file)
        max_weight = snp_df['Weight'].max()
        focal_rows = snp_df[snp_df['Weight'] == max_weight]
        for _, row in focal_rows.iterrows():
            focal_snp_ids.add(row['SNP'])
        true_effects_per_rep[rep_idx] = dict(zip(focal_rows['SNP'], focal_rows['Effect']))

print(f"  Focal SNPs identified: {sorted(focal_snp_ids)}")

# Collect per-replicate chi-squared for each focal SNP
records = []
for f in blue_files:
    rep = int(f.split('/rep_')[1].split('/')[0])
    df = pd.read_csv(f)
    focal = df[df['covar'] != 'intercept'][['covar', 'blue', 'se', 'pval']].copy()
    focal['blue'] = -focal['blue']  # negate: SLEMM uses A1 coding, sim uses A2
    focal['chi2'] = (focal['blue'] / focal['se']) ** 2
    focal['rep'] = rep

    if rep in true_effects_per_rep:
        focal['true_effect'] = focal['covar'].map(true_effects_per_rep[rep])
    else:
        focal['true_effect'] = np.nan

    records.append(focal)

all_df = pd.concat(records, ignore_index=True)

# --- BLUE unbiasedness check ---
print(f"\n{'='*80}")
print("BLUE UNBIASEDNESS CHECK")
print(f"{'='*80}")

bias_records = []
for snp in sorted(all_df['covar'].unique()):
    sub = all_df[all_df['covar'] == snp]
    true_eff = sub['true_effect'].iloc[0]
    mean_blue = sub['blue'].mean()
    sd_blue = sub['blue'].std()
    se_of_mean = sd_blue / np.sqrt(len(sub))
    ratio = mean_blue / true_eff if abs(true_eff) > 1e-10 else np.nan
    t_stat = (mean_blue - true_eff) / se_of_mean if se_of_mean > 0 else np.nan

    bias_records.append({
        'covar': snp,
        'true_effect': true_eff,
        'mean_blue': mean_blue,
        'sd_blue': sd_blue,
        'se_of_mean': se_of_mean,
        'ratio_blue_over_true': ratio,
        't_stat_bias': t_stat,
        'nrep': len(sub),
    })

    print(f"  {snp}:")
    print(f"    true_effect  = {true_eff:12.4f}")
    print(f"    mean(BLUE)   = {mean_blue:12.4f} +/- {se_of_mean:.4f}")
    print(f"    ratio        = {ratio:12.4f}  (should be ~1.0 if unbiased)")
    print(f"    t-stat(bias) = {t_stat:12.4f}  (|t|>2 suggests bias)")

bias_df = pd.DataFrame(bias_records)
bias_df.to_csv(f"{simdir}/blue_bias_check.csv", index=False)

# --- Summarize per focal SNP ---
summary = all_df.groupby('covar').agg(
    mean_blue=('blue', 'mean'),
    sd_blue=('blue', 'std'),
    mean_se=('se', 'mean'),
    mean_chi2=('chi2', 'mean'),
    sd_chi2=('chi2', 'std'),
    median_chi2=('chi2', 'median'),
    max_chi2=('chi2', 'max'),
    min_chi2=('chi2', 'min'),
    ncp_empirical=('chi2', lambda x: x.mean() - 1),
    mean_pval=('pval', 'mean'),
    median_pval=('pval', 'median'),
    n_sig_5e8=('pval', lambda x: (x < 5e-8).sum()),
    nrep=('rep', 'count'),
).reset_index()

# --- Merge theoretical NCP ---
if cl_data is not None:
    summary['c_l'] = summary['covar'].map(lambda s: cl_data.get(s, {}).get('c_l', np.nan))
    summary['c_l_N'] = summary['covar'].map(lambda s: cl_data.get(s, {}).get('c_l_N', np.nan))
    summary['ratio_cl_cbar'] = summary['covar'].map(lambda s: cl_data.get(s, {}).get('ratio_cl_cbar', np.nan))
    summary['c_bar'] = cbar
    summary['S_N'] = S_N
    summary['Wss'] = Wss
    summary['Wtilde'] = Wtilde
    summary['M'] = M
    summary['ncp_theoretical_cl'] = summary['c_l'] * N * h2 / M * Wtilde
    summary['ncp_theoretical_cbar'] = S_N / M * Wtilde
    summary['ncp_ratio'] = summary['ncp_empirical'] / summary['ncp_theoretical_cl']

summary.to_csv(f"{simdir}/summary.csv", index=False)

# --- Print comparison table ---
print(f"\n{'='*100}")
print("EMPIRICAL vs THEORETICAL NCP")
print(f"{'='*100}")

if cl_data is not None:
    print(f"  M={M}, h2={h2}, beta2={beta2}, Wss={Wss:.2f}, Wtilde={Wtilde:.2f}")
    print(f"  S(N)={S_N:.2f}, c_bar={cbar:.6f}")
    print(f"  NCP(c_bar) = S(N)/M * Wtilde = {S_N/M * Wtilde:.2f}")
    print()
    print(f"  {'SNP':<25} {'cl/cbar':>8} {'NCP(cl)':>10} {'NCP(cbar)':>10} "
          f"{'NCP(emp)':>10} {'emp/cl':>8} {'power':>8}")
    print(f"  {'-'*90}")
    for _, row in summary.iterrows():
        ncp_cl = row.get('ncp_theoretical_cl', np.nan)
        ncp_cb = row.get('ncp_theoretical_cbar', np.nan)
        ratio = row.get('ncp_ratio', np.nan)
        print(f"  {row['covar']:<25} {row.get('ratio_cl_cbar', np.nan):8.4f} "
              f"{ncp_cl:10.2f} {ncp_cb:10.2f} {row['ncp_empirical']:10.2f} "
              f"{ratio:8.4f} {row['n_sig_5e8']:>5.0f}/{row['nrep']:.0f}")
    # Average across all focal SNPs
    print(f"  {'-'*90}")
    mean_ratio_cl = summary['ratio_cl_cbar'].mean()
    mean_ncp_cl = summary['ncp_theoretical_cl'].mean()
    mean_ncp_cb = summary['ncp_theoretical_cbar'].mean()
    mean_ncp_emp = summary['ncp_empirical'].mean()
    mean_ncp_ratio = summary['ncp_ratio'].mean()
    total_sig = summary['n_sig_5e8'].sum()
    total_n = summary['nrep'].sum()
    print(f"  {'MEAN':<25} {mean_ratio_cl:8.4f} "
          f"{mean_ncp_cl:10.2f} {mean_ncp_cb:10.2f} {mean_ncp_emp:10.2f} "
          f"{mean_ncp_ratio:8.4f} {total_sig:>5.0f}/{total_n:.0f}")
else:
    print(f"  {'SNP':<25} {'E[chi2]':>8} {'NCP':>8} {'SD':>8} {'power':>8}")
    print(f"  {'-'*60}")
    for _, row in summary.iterrows():
        print(f"  {row['covar']:<25} {row['mean_chi2']:8.2f} {row['ncp_empirical']:8.2f} "
              f"{row['sd_chi2']:8.2f} {row['n_sig_5e8']:>5.0f}/{row['nrep']:.0f}")

# --- Per-replicate detail ---
detail_path = f"{simdir}/per_replicate_chi2.csv"
pivot = all_df.pivot(index='rep', columns='covar', values='chi2').reset_index()
pivot.to_csv(detail_path, index=False)

print(f"\nSaved: {simdir}/summary.csv")
print(f"Saved: {simdir}/blue_bias_check.csv")
print(f"Saved: {detail_path}")

