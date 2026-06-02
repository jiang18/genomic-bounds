#!/usr/bin/env python3
"""
compute_cbar_loco.py

Compute full-GRM and LOCO average GRAMMAR-Gamma coefficients (c_bar) as a
function of sample size N.  Produces CSV input for the divergence figure.

  c_bar_full  : from stored eigenvalues (no G needed)
  c_bar_LOCO  : Cholesky-based x^T V^{-1} x, leaving out chromosome 1

Dual trick for LOCO:
  N <= M_LOCO  :  Cholesky of  N x N  matrix V_LOCO
  N >  M_LOCO  :  Woodbury  -> Cholesky of  M_LOCO x M_LOCO  matrix

Usage
-----
  python3 compute_cbar_loco.py <sim_dir> <eigen_dir> <out_csv> [N1 N2 ...]

  sim_dir   : directory containing G.npy, chrom_ids.npy  (from simulate_genotypes.py)
  eigen_dir : directory containing eigenvalues_N*.npy, indices_N*.npy
  out_csv   : output CSV file

Example
-------
  python3 compute_cbar_loco.py sim_out/Ne50 eigen_out/Ne50 cbar_Ne50.csv \
      5000 10000 20000 50000 100000
"""

import os
import sys
import time
import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular


# ---------------------------------------------------------------------------
# Heritability values to evaluate
# ---------------------------------------------------------------------------
H2_VALUES = [0.1, 0.3, 0.5]

# Chromosome to leave out (1-based)
LOCO_CHR = 1


# ---------------------------------------------------------------------------
# c_bar_full from eigenvalues  (no G needed)
# ---------------------------------------------------------------------------
def cbar_full_from_eigenvalues(eigen_dir, N, M, h2_values):
    """
    c_bar_full = (1/N) * sum_k  d_k / (d_k * h2 + 1 - h2)

    where d_k = N * ell_k / M  are GRM eigenvalues.
    """
    ell = np.load(os.path.join(eigen_dir, f"eigenvalues_N{N}.npy"))
    d = ell * (N / M)                                   # GRM eigenvalues

    results = {}
    for h2 in h2_values:
        cbar = np.sum(d / (d * h2 + 1.0 - h2)) / N
        results[h2] = cbar
    return results


# ---------------------------------------------------------------------------
# c_bar_LOCO via Cholesky   (needs G)
# ---------------------------------------------------------------------------
def cbar_loco_cholesky(Z_loco, Z_chr1, N, h2_values):
    """
    Compute  c_bar_LOCO  for each h2 value using Cholesky.

    c_bar_LOCO = (1 / (M_chr1 * N)) * tr( Z_chr1^T  V_LOCO^{-1}  Z_chr1 )

    Dual trick:
      N <= M_LOCO  ->  Cholesky of N x N matrix V_LOCO
      N >  M_LOCO  ->  Woodbury   -> Cholesky of M_LOCO x M_LOCO matrix
    """
    M_loco = Z_loco.shape[1]
    M_chr1 = Z_chr1.shape[1]
    use_woodbury = (N > M_loco)

    results = {}

    if use_woodbury:
        # Pre-compute products shared across h2 values
        # ZtZ_loco : M_LOCO x M_LOCO
        # B        : M_LOCO x M_chr1
        print(f"    Woodbury path (N={N:,} > M_LOCO={M_loco:,})")
        print(f"    Computing Z_loco^T @ Z_loco  ({M_loco:,} x {M_loco:,}) ...",
              flush=True)
        ZtZ_loco = Z_loco.T @ Z_loco                    # M_LOCO x M_LOCO

        print(f"    Computing Z_loco^T @ Z_chr1  ({M_loco:,} x {M_chr1:,}) ...",
              flush=True)
        B = Z_loco.T @ Z_chr1                            # M_LOCO x M_chr1

        for h2 in h2_values:
            h2_loco = h2 * M_loco / (M_loco + M_chr1)
            gamma = h2_loco / (M_loco * (1.0 - h2_loco))

            # C = I + gamma * Z_loco^T Z_loco   (M_LOCO x M_LOCO)
            C = np.eye(M_loco) + gamma * ZtZ_loco
            print(f"    h2={h2}: Cholesky of {M_loco:,} x {M_loco:,} ...",
                  flush=True)
            t0 = time.time()
            L_C = np.linalg.cholesky(C)
            print(f"      Cholesky: {time.time()-t0:.1f}s", flush=True)

            # W = L_C^{-1} B    (M_LOCO x M_chr1)
            W = solve_triangular(L_C, B, lower=True)

            cbar = (1.0 / (1.0 - h2_loco)
                    - gamma / (M_chr1 * N * (1.0 - h2_loco))
                    * np.sum(W ** 2))
            results[h2] = cbar

        del ZtZ_loco, B

    else:
        # Direct path: Cholesky of  N x N  matrix V_LOCO
        print(f"    Direct path (N={N:,} <= M_LOCO={M_loco:,})")
        print(f"    Computing GRM_LOCO  ({N:,} x {N:,}) ...", flush=True)
        GRM_loco = Z_loco @ Z_loco.T / M_loco           # N x N

        for h2 in h2_values:
            h2_loco = h2 * M_loco / (M_loco + M_chr1)

            # V_LOCO = GRM_LOCO * h2_LOCO + I * (1 - h2_LOCO)
            V = GRM_loco * h2_loco
            np.fill_diagonal(V, V.diagonal() + (1.0 - h2_loco))

            print(f"    h2={h2}: Cholesky of {N:,} x {N:,} ...", flush=True)
            t0 = time.time()
            L = np.linalg.cholesky(V)
            print(f"      Cholesky: {time.time()-t0:.1f}s", flush=True)

            # W = L^{-1} Z_chr1   (N x M_chr1)
            W = solve_triangular(L, Z_chr1, lower=True)

            cbar = np.sum(W ** 2) / (M_chr1 * N)
            results[h2] = cbar

        del GRM_loco

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) < 4:
        print(__doc__.strip())
        sys.exit(2)

    sim_dir   = sys.argv[1]
    eigen_dir = sys.argv[2]
    out_csv   = sys.argv[3]
    N_values  = [int(x) for x in sys.argv[4:]] if len(sys.argv) > 4 \
                else [5_000, 10_000, 20_000, 50_000, 100_000]

    # --- Load full genotype matrix and chromosome ids ---
    print("Loading genotypes ...", flush=True)
    G_full    = np.load(os.path.join(sim_dir, "G.npy"))       # (N_max, M)
    chrom_ids = np.load(os.path.join(sim_dir, "chrom_ids.npy"))
    N_max, M  = G_full.shape
    print(f"  G shape: {N_max:,} x {M:,}")

    chr1_mask = (chrom_ids == LOCO_CHR)
    M_chr1    = int(chr1_mask.sum())
    M_loco    = M - M_chr1
    print(f"  Chr{LOCO_CHR}: {M_chr1:,} SNPs;  LOCO: {M_loco:,} SNPs")

    # Validate N values
    N_values = sorted(set(N for N in N_values if N <= N_max))
    print(f"  N values: {N_values}")

    # --- Process each N ---
    rows = []

    for N in N_values:
        print(f"\n{'='*60}")
        print(f"N = {N:,}")
        print(f"{'='*60}")

        # Load subsample indices (same as eigenvalue computation)
        idx_file = os.path.join(eigen_dir, f"indices_N{N}.npy")
        if os.path.exists(idx_file):
            idx = np.load(idx_file)
            print(f"  Loaded indices from {idx_file}")
        else:
            rng = np.random.default_rng(0)
            idx = np.sort(rng.choice(N_max, size=N, replace=False))
            print(f"  Generated fresh indices (seed=0)")

        G = G_full[idx, :]                                  # (N, M)

        # Re-standardise within subsample
        mu = G.mean(axis=0)
        sd = G.std(axis=0)
        sd[sd == 0] = 1.0
        G = (G - mu) / sd

        # Split into LOCO and chr1
        Z_chr1 = G[:, chr1_mask].copy()                    # (N, M_chr1)
        Z_loco = G[:, ~chr1_mask].copy()                   # (N, M_loco)
        del G

        # c_bar_full from eigenvalues
        print("  c_bar_full from eigenvalues ...", flush=True)
        cbar_full = cbar_full_from_eigenvalues(
            eigen_dir, N, M, H2_VALUES
        )

        # c_bar_LOCO via Cholesky
        print("  c_bar_LOCO via Cholesky ...", flush=True)
        cbar_loco = cbar_loco_cholesky(
            Z_loco, Z_chr1, N, H2_VALUES
        )
        del Z_loco, Z_chr1

        for h2 in H2_VALUES:
            h2_loco = h2 * M_loco / M
            row = {
                "N": N,
                "M": M,
                "M_loco": M_loco,
                "M_chr1": M_chr1,
                "h2": h2,
                "h2_loco": round(h2_loco, 6),
                "cbar_full": round(cbar_full[h2], 6),
                "cbar_loco": round(cbar_loco[h2], 6),
                "loco_upper_bound": round(1.0 / (1.0 - h2_loco), 6),
            }
            rows.append(row)
            print(f"  h2={h2:.1f}: c_bar_full={cbar_full[h2]:.4f}, "
                  f"c_bar_LOCO={cbar_loco[h2]:.4f}, "
                  f"bound={1/(1-h2_loco):.4f}")

    del G_full

    # --- Write CSV ---
    keys = list(rows[0].keys())
    with open(out_csv, "w") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in keys) + "\n")

    print(f"\nDone. Saved {len(rows)} rows to {out_csv}")


if __name__ == "__main__":
    main()
