#!/usr/bin/env python3
"""
compute_eigenvalues_by_N.py

For a fixed simulated genotype matrix (Ne, M fixed), compute LD eigenvalues
at each of several sample sizes N by subsampling rows of G.

Each N gets its own eigenvalue array, correctly reflecting the rank and
shrinkage at that sample size. This is the right input for a sigmoid sum
trajectory where eigenvalues are N-appropriate.

Strategy by matrix size:
  - min(N, M) <= MAX_EXACT_DIM : exact eigvalsh via dual trick
  - min(N, M) >  MAX_EXACT_DIM : randomized SVD (top N_COMPONENTS singular
                                  values), capturing all eigenvalues > ~0.01

For the randomized SVD path, summary statistics are computed correctly by
using M_total (= tr(R) = M) as the known total variance rather than the
truncated sum, so participation ratio and EIG thresholds remain accurate.

Individual indices for each subsample are saved as indices_N{n}.npy so that
downstream analyses (phenotype simulation, REML) use exactly the same
individuals.

Outputs (all in outdir/):
  eigenvalues_N{n}.npy   : LD eigenvalues at sample size n, descending
  indices_N{n}.npy       : sorted row indices into the full genotype matrix
  eigen_meta.csv         : one row per N with summary statistics

Usage:
  python3 compute_eigenvalues_by_N.py <geno_dir> <outdir> [N1 N2 ...]

  If N values are omitted, defaults to: 5000 10000 20000 50000 100000

Example:
  python3 compute_eigenvalues_by_N.py sim_out/Ne50_N200k eigen_out/Ne50 \
      5000 10000 20000 50000 100000 200000
"""

import os
import sys
import numpy as np

try:
    from sklearn.utils.extmath import randomized_svd
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

DEFAULT_N_VALUES = [5_000, 10_000, 20_000, 50_000, 100_000]

# Exact eigvalsh is feasible up to this matrix dimension.
# 50k x 50k float64 ~ 20 GB; 30k x 30k ~ 7 GB.
# Adjust downward if RAM is tight.
MAX_EXACT_DIM = 60_000

# Number of singular values to compute in randomized SVD.
# Must capture all eigenvalues meaningfully above 0.
# At N=50k with M=100k, EIG99.9% ~ 10,648; 20,000 gives ample margin.
N_COMPONENTS = 20_000


def compute_ld_eigenvalues(G, M_total, max_exact_dim=MAX_EXACT_DIM,
                           n_components=N_COMPONENTS):
    """
    Eigenvalues of the M x M LD matrix R = Z^T Z / N.

    For small min(N, M): exact computation via dual trick.
    For large min(N, M): randomized SVD returning top n_components values;
                         remaining eigenvalues are set to zero.

    Parameters
    ----------
    G        : (N, M) standardized genotype array
    M_total  : int, total number of markers (= tr(R), used to normalize
               summary statistics correctly when SVD truncates the spectrum)

    Returns
    -------
    ell      : descending array of length min(N, M), clipped to >= 0
    truncated: bool, True if randomized SVD was used
    """
    N, M = G.shape
    small_dim = min(N, M)

    if small_dim <= max_exact_dim:
        # --- Exact path via dual trick ---
        if N <= M:
            # Form N x N GRM (smaller matrix), then rescale eigenvalues
            GRM = (G @ G.T) / M
            d = np.linalg.eigvalsh(GRM)[::-1].copy()
            ell = d * (M / N)
        else:
            # Form M x M LD matrix directly
            R = (G.T @ G) / N
            ell = np.linalg.eigvalsh(R)[::-1].copy()
        return np.clip(ell, 0, None), False

    else:
        # --- Randomized SVD path ---
        if not HAS_SKLEARN:
            raise ImportError(
                "scikit-learn is required for large N/M. "
                "Install with: pip install scikit-learn"
            )
        n_comp = min(n_components, small_dim)
        print(f"  Using randomized SVD (n_components={n_comp}) ...", flush=True)

        # SVD of G / sqrt(M): singular values s_k satisfy ell_k = s_k^2 * M/N
        _, s, _ = randomized_svd(
            G,
            n_components=n_comp,
            n_iter=1,
            power_iteration_normalizer='QR',
            random_state=0
        )

        ell_top = np.clip(s**2 / N, 0, None)

        # Remaining eigenvalues are near zero; pad so len = small_dim
        # Their sum = M_total - ell_top.sum() distributed over (small_dim - n_comp) values
        ell = np.concatenate([ell_top, np.zeros(small_dim - n_comp)])
        return ell, True


def variance_thresholds(ell, M_total, thresholds=(0.90, 0.95, 0.98, 0.99, 0.999)):
    """
    Number of eigenvalues needed to explain fraction p of total variance.
    Uses M_total (the known trace) rather than ell.sum() so that truncated
    spectra give correct fractions.
    """
    cumvar = np.cumsum(ell) / M_total
    return {p: int(np.searchsorted(cumvar, p)) + 1 for p in thresholds}


def participation_ratio(ell, M_total):
    """
    PR = (sum ell_k)^2 / sum(ell_k^2).
    Uses M_total for the numerator so truncation does not bias the result.
    """
    return float(M_total)**2 / float(np.sum(ell**2))


def main():
    if len(sys.argv) < 3:
        print(__doc__.strip())
        sys.exit(2)

    geno_dir = sys.argv[1]
    outdir   = sys.argv[2]
    N_values = [int(x) for x in sys.argv[3:]] if len(sys.argv) > 3 else DEFAULT_N_VALUES
    os.makedirs(outdir, exist_ok=True)

    # --- Load full genotype matrix ---
    print("Loading genotypes ...", flush=True)
    G_full = np.load(os.path.join(geno_dir, "G.npy"))   # (N_full, M)
    N_full, M = G_full.shape
    print(f"  Full matrix: N={N_full:,}, M={M:,}")
    print(f"  MAX_EXACT_DIM = {MAX_EXACT_DIM:,}, N_COMPONENTS = {N_COMPONENTS:,}")

    # Validate requested N values
    for N in N_values:
        if N > N_full:
            print(f"  WARNING: requested N={N:,} > available N={N_full:,}; skipping.")
    N_values = sorted(set(N for N in N_values if N <= N_full))
    print(f"  N values to process: {N_values}")

    rng = np.random.default_rng(0)
    summary_rows = []

    for N in N_values:
        print(f"\nN = {N:,}", flush=True)

        # Subsample rows (individuals) without replacement
        idx = np.sort(rng.choice(N_full, size=N, replace=False))
        G = G_full[idx, :]                               # (N, M)

        # Save indices for downstream use (phenotype simulation, REML)
        idx_fname = os.path.join(outdir, f"indices_N{N}.npy")
        np.save(idx_fname, idx)
        print(f"  Saved individual indices -> {idx_fname}")

        # Re-standardise within the subsample
        mu = G.mean(axis=0)
        sd = G.std(axis=0)
        sd[sd == 0] = 1.0
        G = (G - mu) / sd

        small_dim = min(N, M)
        print(f"  Computing eigenvalues ({small_dim:,}x{small_dim:,} matrix) ...",
              flush=True)

        ell, truncated = compute_ld_eigenvalues(G, M_total=M)
        del G

        # Save full eigenvalue array
        fname = os.path.join(outdir, f"eigenvalues_N{N}.npy")
        np.save(fname, ell)

        # --- Summaries (all use M as the known total variance) ---
        n_nonzero  = int((ell > 1e-10).sum())
        pr         = participation_ratio(ell, M_total=M)
        vt         = variance_thresholds(ell, M_total=M)

        # Minimum eigenvalue that is meaningfully above zero
        meaningful  = ell[ell > 0.01]
        min_meaningful = float(meaningful.min()) if len(meaningful) > 0 else float('nan')

        print(f"  Truncated (randomized SVD) : {truncated}")
        print(f"  Nonzero eigenvalues (>1e-10): {n_nonzero:,}")
        print(f"  Eigenvalues > 0.01          : {len(meaningful):,}")
        print(f"  Min eigenvalue > 0.01       : {min_meaningful:.4f}")
        print(f"  Max eigenvalue              : {ell[0]:.2f}")
        print(f"  Sum (should = {M:,})        : {ell.sum():.1f}")
        print(f"  Participation ratio         : {pr:.1f}")
        for p, k in vt.items():
            print(f"  EIG{p*100:.1f}%                  : {k:,}")

        row = {
            "N": N, "M": M,
            "truncated": truncated,
            "n_nonzero": n_nonzero,
            "n_above_001": len(meaningful),
            "min_above_001": min_meaningful,
            "max_eigenvalue": ell[0],
            "participation_ratio": pr,
        }
        row.update({f"EIG{p:.3f}": k for p, k in vt.items()})
        summary_rows.append(row)

    del G_full  # intentional: free memory after all N values processed

    # --- Save summary ---
    meta_path = os.path.join(outdir, "eigen_meta.csv")
    if HAS_PANDAS:
        pd.DataFrame(summary_rows).to_csv(meta_path, index=False)
    else:
        keys = list(summary_rows[0].keys())
        with open(meta_path, "w") as f:
            f.write(",".join(keys) + "\n")
            for row in summary_rows:
                f.write(",".join(str(row[k]) for k in keys) + "\n")

    print(f"\nDone. Results saved to {outdir}/")
    print("  eigenvalues_N{{n}}.npy  for each N")
    print("  indices_N{{n}}.npy     for each N")
    print("  eigen_meta.csv")


if __name__ == "__main__":
    main()
