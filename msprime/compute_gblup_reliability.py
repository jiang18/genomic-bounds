#!/usr/bin/env python3
"""
compute_gblup_reliability.py

Empirical vs theoretical in-sample GBLUP reliability (validates Equation 17,
R_in^2 = 1 - lambda*S(N)/N = 1 - cbar_full*(1 - h2)).

Uses the same phenotype-simulation data that validated the per-SNP NCP. For the
polygenic-dominated scenario (beta2 = 0.001: five focal QTLs each 0.1% PVE +
polygenic background), for each (Ne, h2, N) and each replicate it computes:

  Empirical R_in^2 = cor(GEBV_hat, GV_true)^2  (for BLUP, = Var(GEBV)/Var(GV)),
       GV_true = rep_*/sim.gv.csv = TOTAL true GV (polygenic + 5 focal QTLs).
         R2_random : GEBV_hat = polygenic random prediction only
                     (slemm --pred on rep_*/reml.reml.snp.csv).
         R2_full   : GEBV_hat = random + fixed; the fixed part is the 5-QTL BLUE
                     (rep_*/reml.reml.blue.csv col 'blue') times the focal
                     genotypes (N_*/focal_snp_genotypes.csv). Because the 5 QTLs
                     are fixed-effect covariates in REML, the full prediction
                     matches the total GV -> R2_full is the apples-to-apples
                     in-sample reliability.

  Theoretical R_in^2 = 1 - cbar_full * (1 - h2),  cbar_full from cbar_<Ne>.csv at
       the matching (N, h2), at nominal h2 (effective GRM heritability is
       h2 - 5*beta2 = h2 - 0.005, close enough).

Replicates are averaged (mean +/- SD): R2_theory is fixed per scenario (fixed
genotypes + nominal h2); the empirical value varies across reps from the random
effect/residual draws and the per-rep REML h2 estimate. Per-replicate values are
written to --out_perrep so error bars / a scatter vs theory can be made later.

Alignment: gebv.csv, sim.gv.csv, focal_snp_genotypes.csv all derive from sub_N in
PLINK .fam order, so individuals are aligned by ROW ORDER (the files use
inconsistent IID labels). Lengths are checked.

The per-rep slemm --pred calls are independent and run concurrently (--jobs).

Usage
-----
  python3 compute_gblup_reliability.py --nrep 1                 # quick test
  python3 compute_gblup_reliability.py --nrep 100 --jobs 28     # full run

Run from the directory containing plink_out/ and cbar_Ne*.csv. Needs slemm on PATH.
"""

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd


def load_cbar(cbar_dir, ne):
    df = pd.read_csv(os.path.join(cbar_dir, f"cbar_{ne}.csv"))
    return {(int(r.N), float(r.h2)): float(r.cbar_full) for r in df.itertuples()}


def run_slemm_pred(slemm, bfile, snp_estimate, out):
    subprocess.run(
        [slemm, "--pred", "--bfile", bfile, "--snp_estimate", snp_estimate, "--out", out],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def values(path):
    """Second column as a row-ordered float array (the GEBV/GV value)."""
    return pd.read_csv(path).iloc[:, 1].to_numpy(dtype=float)


def fixed_prediction(blue_csv, focal_geno_csv):
    """Sum_s (focal genotype_s * BLUE_s) over the 5 focal QTLs, row-ordered.

    BLUE from reml.reml.blue.csv (col 'blue' indexed by 'covar'); focal genotypes
    (raw 0/1/2) from focal_snp_genotypes.csv. Intercept omitted (irrelevant to cor).
    """
    blue = pd.read_csv(blue_csv)
    bmap = {str(c): float(b) for c, b in zip(blue["covar"], blue["blue"])}
    geno = pd.read_csv(focal_geno_csv)
    pred = np.zeros(len(geno))
    for c in (col for col in geno.columns if col.startswith("snp_")):
        if c in bmap:
            pred += geno[c].to_numpy(dtype=float) * bmap[c]
        else:
            sys.stderr.write(f"    [warn] focal SNP {c} not in {blue_csv}\n")
    return pred


def r2(x, y):
    return float(np.corrcoef(x, y)[0, 1] ** 2)


def process_rep(task):
    """One (scenario, rep) work unit: slemm --pred + R2_random, R2_full. Thread-safe."""
    ne, h2, N, rep, bfile, simdir, focal_geno, slemm = task
    repdir = os.path.join(simdir, f"rep_{rep}")
    snp_est = os.path.join(repdir, "reml.reml.snp.csv")
    blue_csv = os.path.join(repdir, "reml.reml.blue.csv")
    gv_true = os.path.join(repdir, "sim.gv.csv")
    if not all(os.path.exists(p) for p in (snp_est, blue_csv, gv_true)):
        return None
    gebv = os.path.join(repdir, "gebv.csv")          # unique per rep -> safe concurrently
    try:
        run_slemm_pred(slemm, bfile, snp_est, gebv)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        sys.stderr.write(f"    [warn] slemm --pred failed ({repdir}): {e}\n")
        return None
    g_rand, g_tot, g_fix = values(gebv), values(gv_true), fixed_prediction(blue_csv, focal_geno)
    if not (len(g_rand) == len(g_tot) == len(g_fix)):
        sys.stderr.write(f"    [warn] length mismatch in {repdir}; skipping\n")
        return None
    return dict(Ne=ne, h2=h2, N=N, rep=rep,
                R2_random=r2(g_rand, g_tot), R2_full=r2(g_rand + g_fix, g_tot))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plink_dir", default="./plink_out")
    ap.add_argument("--cbar_dir", default=".")
    ap.add_argument("--ne", nargs="+", default=["Ne20", "Ne50", "Ne100"])
    ap.add_argument("--h2", nargs="+", type=float, default=[0.1, 0.3, 0.5])
    ap.add_argument("--N", nargs="+", type=int, default=[5000, 10000, 20000, 50000])
    ap.add_argument("--beta2", type=float, default=0.001)
    ap.add_argument("--nrep", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=8,
                    help="concurrent slemm --pred jobs (try 28 on the cluster)")
    ap.add_argument("--slemm", default="slemm")
    ap.add_argument("--out", default="gblup_reliability.csv")
    ap.add_argument("--out_perrep", default="gblup_reliability_perrep.csv")
    args = ap.parse_args()

    # --- build the flat (scenario, rep) task list; theory once per scenario ---
    theory, tasks = {}, []
    for ne in args.ne:
        try:
            cbar = load_cbar(args.cbar_dir, ne)
        except FileNotFoundError:
            sys.stderr.write(f"  [skip] cbar_{ne}.csv not found in {args.cbar_dir}\n")
            continue
        for h2 in args.h2:
            for N in args.N:
                if (N, h2) not in cbar:
                    sys.stderr.write(f"  [skip] no cbar_full for {ne} N={N} h2={h2}\n")
                    continue
                simdir = os.path.join(args.plink_dir, ne,
                                      f"h2_{h2}__beta2_{args.beta2}", f"N_{N}")
                bfile = os.path.join(simdir, f"sub_{N}")
                focal_geno = os.path.join(simdir, "focal_snp_genotypes.csv")
                if not (os.path.exists(bfile + ".bim") and os.path.exists(focal_geno)):
                    sys.stderr.write(f"  [skip] missing sub_{N}.bim / focal_snp_genotypes.csv in {simdir}\n")
                    continue
                theory[(ne, h2, N)] = 1.0 - cbar[(N, h2)] * (1.0 - h2)
                for rep in range(1, args.nrep + 1):
                    tasks.append((ne, h2, N, rep, bfile, simdir, focal_geno, args.slemm))

    if not tasks:
        print("No tasks (check paths / data).")
        return

    # --- run replicates concurrently ---
    print(f"Running {len(tasks)} (scenario, rep) jobs on {args.jobs} workers ...", flush=True)
    results, done = [], 0
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(process_rep, t) for t in tasks]
        for fut in as_completed(futs):
            r = fut.result()
            if r:
                results.append(r)
            done += 1
            if done % 25 == 0 or done == len(tasks):
                print(f"  {done}/{len(tasks)} done", flush=True)

    perrep = pd.DataFrame(results)
    perrep.to_csv(args.out_perrep, index=False)

    # --- per-scenario summary (mean +/- SD) in scenario order ---
    rows = []
    for (ne, h2, N), R2_th in theory.items():
        g = perrep[(perrep.Ne == ne) & (perrep.h2 == h2) & (perrep.N == N)]
        if len(g) == 0:
            continue
        rows.append(dict(
            Ne=ne, h2=h2, N=N, nrep=len(g),
            R2_random_mean=round(g.R2_random.mean(), 4),
            R2_full_mean=round(g.R2_full.mean(), 4),
            R2_full_sd=round(g.R2_full.std(ddof=1), 4) if len(g) > 1 else 0.0,
            R2_theory=round(R2_th, 4),
            ratio_full_over_theory=round(g.R2_full.mean() / R2_th, 3) if R2_th > 0 else float("nan"),
        ))
    summary = pd.DataFrame(rows)
    summary.to_csv(args.out, index=False)

    print(f"\n{'Ne':>5} {'h2':>4} {'N':>7} {'nrep':>4} "
          f"{'R2_random':>9} {'R2_full':>8} {'+/-SD':>7} {'R2_theory':>9} {'full/th':>7}")
    for r in rows:
        print(f"{r['Ne']:>5} {r['h2']:>4} {r['N']:>7} {r['nrep']:>4} "
              f"{r['R2_random_mean']:>9.4f} {r['R2_full_mean']:>8.4f} {r['R2_full_sd']:>7.4f} "
              f"{r['R2_theory']:>9.4f} {r['ratio_full_over_theory']:>7.3f}")
    print(f"\nSaved summary -> {args.out}  (per-replicate -> {args.out_perrep})")


if __name__ == "__main__":
    main()
