#!/usr/bin/env python3
"""
project_focal_snps.py — Eigendecompose GRM and project focal SNPs.

For each sample size N (up to ~50k), performs the N×N GRM eigendecomposition
ONCE, projects each focal SNP onto the eigenvectors ONCE, and saves:

    d_k  : GRM eigenvalues (length N)
    a_lk : projections u_k' z_l for each focal SNP (n_focal × N)

These are purely geometric quantities independent of h2. The h2-dependent
c_l computation is deferred to collect_results.py.

Usage:
    python3 project_focal_snps.py <geno_dir> <eigen_dir> <plink_dir> [N1 N2 ...]

    If N values are omitted, uses all indices_N*.npy files with N <= 50000.

Example:
    python3 project_focal_snps.py sim_out/Ne20 eigen_out/Ne20 plink_out/Ne20
    python3 project_focal_snps.py sim_out/Ne20 eigen_out/Ne20 plink_out/Ne20 5000 10000 20000 50000
"""

import os
import sys
import glob
import numpy as np
import pandas as pd

MAX_N = 50_000  # maximum N for N×N eigendecomposition


def main():
    if len(sys.argv) < 4:
        print(__doc__.strip())
        sys.exit(2)

    geno_dir = sys.argv[1]
    eigen_dir = sys.argv[2]
    plink_dir = sys.argv[3]

    # Determine N values
    if len(sys.argv) > 4:
        N_values = [int(x) for x in sys.argv[4:]]
    else:
        idx_files = glob.glob(os.path.join(eigen_dir, "indices_N*.npy"))
        N_values = []
        for f in idx_files:
            n = int(os.path.basename(f).replace("indices_N", "").replace(".npy", ""))
            if n <= MAX_N:
                N_values.append(n)
        N_values = sorted(N_values)

    print(f"=== project_focal_snps ===")
    print(f"  geno_dir:  {geno_dir}")
    print(f"  eigen_dir: {eigen_dir}")
    print(f"  plink_dir: {plink_dir}")
    print(f"  N values:  {N_values}")

    # --- Load focal SNP info ---
    focal_info = pd.read_csv(os.path.join(plink_dir, "focal_snps_info.csv"))
    focal_ids = focal_info['SNP'].tolist()
    n_focal = len(focal_ids)
    print(f"  Focal SNPs ({n_focal}): {focal_ids}")

    # --- Reconstruct SNP IDs to find focal SNP column indices ---
    positions = np.load(os.path.join(geno_dir, "positions.npy"))
    chrom_ids = np.load(os.path.join(geno_dir, "chrom_ids.npy"))
    M = len(positions)
    snp_ids = [f"snp_{int(chrom_ids[j])}_{int(positions[j])}" for j in range(M)]

    focal_col_idx = []
    for sid in focal_ids:
        if sid in snp_ids:
            focal_col_idx.append(snp_ids.index(sid))
        else:
            print(f"  WARNING: {sid} not found in SNP IDs")
            focal_col_idx.append(None)

    print(f"  Focal column indices: {focal_col_idx}")

    # --- Load full genotype matrix ---
    print("\nLoading genotypes...", flush=True)
    G_full = np.load(os.path.join(geno_dir, "G.npy"))
    N_full = G_full.shape[0]
    print(f"  Full matrix: N={N_full}, M={M}")

    # --- Process each N ---
    for N in N_values:
        if N > MAX_N:
            print(f"\nSkipping N={N} (> MAX_N={MAX_N})")
            continue

        print(f"\n{'='*60}")
        print(f"N = {N}")
        print(f"{'='*60}")

        # Load indices
        idx_file = os.path.join(eigen_dir, f"indices_N{N}.npy")
        if not os.path.exists(idx_file):
            print(f"  WARNING: {idx_file} not found, skipping")
            continue
        idx = np.load(idx_file)

        # Subsample and standardize
        G = G_full[idx, :]
        mu = G.mean(axis=0)
        sd = G.std(axis=0)
        sd[sd == 0] = 1.0
        Z = (G - mu) / sd  # N x M
        del G

        # --- N×N GRM eigendecomposition ---
        print(f"  Computing GRM eigendecomposition ({N}x{N})...", flush=True)
        GRM = Z @ Z.T / M  # N x N
        d, U = np.linalg.eigh(GRM)
        d = d[::-1].copy()
        U = U[:, ::-1].copy()
        d = np.clip(d, 0, None)
        del GRM

        print(f"  Eigenvalues: min={d.min():.6f}, max={d.max():.4f}, "
              f"sum={d.sum():.2f}, n_nonzero={np.sum(d > 1e-10)}")

        # --- Project focal SNPs ---
        print(f"  Projecting {n_focal} focal SNPs...", flush=True)
        a_matrix = np.zeros((n_focal, N), dtype=np.float64)  # n_focal x N

        for i, sid in enumerate(focal_ids):
            ci = focal_col_idx[i]
            if ci is None:
                print(f"    {sid}: SKIPPED (not found)")
                continue
            z_l = Z[:, ci]
            a = U.T @ z_l  # length N
            a_matrix[i, :] = a

            # Sanity check: ||a||^2 should be ~N (since z_l'z_l = N)
            a_norm2 = np.sum(a**2)
            b2 = N - a_norm2
            print(f"    {sid}: ||a||^2={a_norm2:.2f}, ||b||^2={b2:.6f}")

        del Z, U

        # --- Save ---
        outfile = os.path.join(eigen_dir, f"focal_proj_N{N}.npz")
        np.savez(outfile,
                 d=d,
                 a=a_matrix,
                 focal_ids=np.array(focal_ids),
                 N=N, M=M)
        print(f"  Saved: {outfile}")
        print(f"    d: shape={d.shape}")
        print(f"    a: shape={a_matrix.shape}")

        del d

    del G_full
    print("\n=== Done ===")


if __name__ == "__main__":
    main()

