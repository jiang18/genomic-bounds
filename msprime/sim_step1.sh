#!/bin/bash
# sim_step1.sh — Subset individuals, simulate phenotypes
#
# Usage:
#   bash sim_step1.sh <N> <h2> <beta2> <nrep> <dir> <eigendir>
#
# Arguments:
#   N        — Number of individuals to subset
#   h2       — Total heritability
#   beta2    — Per-focal-SNP PVE (beta_i^2 on standardized scale)
#   nrep     — Number of replicates
#   dir      — Plink directory (e.g., ./plink_out/Ne20)
#   eigendir — Eigenvalue output directory (e.g., ./eigen_out/Ne20)
#              Must contain indices_N{N}.npy from compute_eigenvalues_by_N.py
#
# Requires: <dir>/focal_snps.txt (created by select_focal_snps.sh)
#           <dir>/genotypes.bed/bim/fam
#           <eigendir>/indices_N<N>.npy
#
# Output in <dir>/h2_<h2>__beta2_<beta2>/N_<N>/:
#   sub_<N>.bed/bim/fam       — Subset genotypes
#   sub_<N>.frq               — MAF from PLINK
#   focal_snp_genotypes.csv   — IID, intercept, focal SNP dosages (REML covariate)
#   rep_<r>/sim.snp.csv       — SNP effect file for SLEMM --pred
#   sim.pheno.txt             — FID IID pheno_1 ... pheno_<nrep>

set -euo pipefail

N=$1; h2=$2; beta2=$3; nrep=$4; dir=$5; eigendir=$6
bed_prefix=${dir}/genotypes
outdir="${dir}/h2_${h2}__beta2_${beta2}/N_${N}"
mkdir -p "${outdir}"

focal_snps_file="${dir}/focal_snps.txt"
if [ ! -f "${focal_snps_file}" ]; then
    echo "ERROR: ${focal_snps_file} not found. Run select_focal_snps.sh first."
    exit 1
fi

indices_file="${eigendir}/indices_N${N}.npy"
if [ ! -f "${indices_file}" ]; then
    echo "ERROR: ${indices_file} not found. Run compute_eigenvalues_by_N.py first."
    exit 1
fi

echo "=== sim_step1: N=${N} h2=${h2} beta2=${beta2} nrep=${nrep} ==="
echo "  outdir: ${outdir}"
echo "  indices: ${indices_file}"

# --- Step 1: Extract individuals matching eigenvalue analysis ---
echo "[1/4] Extracting ${N} individuals (matching eigenvalue indices)..."
python3 -c "
import numpy as np, pandas as pd
idx = np.load('${indices_file}')
fam = pd.read_csv('${bed_prefix}.fam', sep=r'\s+', header=None, engine='python')
keep = fam.iloc[idx, [0, 1]]
keep.to_csv('${outdir}/keep.txt', sep=' ', header=False, index=False)
print(f'  Wrote {len(keep)} individuals to keep.txt')
"
plink --chr-set 30 --bfile ${bed_prefix} --keep ${outdir}/keep.txt --make-bed \
      --out ${outdir}/sub_${N} --allow-no-sex --silent

# --- Step 2: Compute MAF ---
echo "[2/4] Computing MAF..."
plink --chr-set 30 --bfile ${outdir}/sub_${N} --freq --out ${outdir}/sub_${N} --allow-no-sex --silent

# --- Step 3: Extract focal SNP genotypes ---
echo "[3/4] Extracting focal SNP genotypes..."
plink --chr-set 30 --bfile ${outdir}/sub_${N} --extract ${focal_snps_file} --recode A \
      --out ${outdir}/focal_geno --allow-no-sex --silent

# --- Step 4: Simulate effects and phenotypes ---
echo "[4/4] Simulating effects and phenotypes..."

python3 - ${N} ${h2} ${beta2} ${nrep} "${outdir}" "${focal_snps_file}" << 'PYSIM'
import numpy as np
import pandas as pd
import os, sys

N = int(sys.argv[1])
h2 = float(sys.argv[2])
beta2 = float(sys.argv[3])
nrep = int(sys.argv[4])
outdir = sys.argv[5]
focal_snps_file = sys.argv[6]

# --- Read inputs ---
with open(focal_snps_file) as f:
    focal_ids = set(l.strip() for l in f if l.strip())
n_focal = len(focal_ids)

bim = pd.read_csv(f"{outdir}/sub_{N}.bim", sep='\t', header=None,
                   names=['Chr','SNP','cM','Pos','A1','A2'])
M = len(bim)

frq = pd.read_csv(f"{outdir}/sub_{N}.frq", sep=r'\s+', engine='python')
maf_map = dict(zip(frq['SNP'], frq['MAF']))

fam = pd.read_csv(f"{outdir}/sub_{N}.fam", sep=r'\s+', header=None,
                   names=['FID','IID','PID','MID','Sex','Pheno'], engine='python')

# --- Convert beta2 (PVE) to Wss ---
# beta2 = Wss * h2 / M  =>  Wss = beta2 * M / h2
Wss = beta2 * M / h2
print(f"  beta2={beta2}, M={M}, h2={h2} => Wss={Wss:.4f}")

# --- Weights: mean(W) = 1 ---
w_bg = (M - n_focal * Wss) / (M - n_focal)
assert w_bg > 0, f"w_bg={w_bg:.4f} <= 0; reduce beta2 (Wss={Wss:.2f} too large for M={M})"

# --- Variance ---
sigma2_alpha = 1.0
sigma2_e = M * (1.0 - h2) / h2

print(f"  M={M}, n_focal={n_focal}, Wss={Wss:.4f}, w_bg={w_bg:.6f}")
print(f"  sigma2_e={sigma2_e:.2f}, total_gvar={M:.0f}, total_pvar={M/h2:.0f}")

# --- Focal SNP signs (random, fixed across replicates) ---
np.random.seed(2025)
focal_signs = {sid: np.random.choice([-1, 1]) for sid in sorted(focal_ids)}

# --- Precompute SNP arrays ---
snp_snps = bim['SNP'].values
snp_chrs = bim['Chr'].values
snp_pos = bim['Pos'].values
snp_a1 = bim['A1'].values
snp_a2 = bim['A2'].values
snp_maf = np.array([maf_map.get(s, 0.0) for s in snp_snps])
snp_pq2 = np.maximum(2.0 * snp_maf * (1.0 - snp_maf), 1e-10)
snp_is_focal = np.array([s in focal_ids for s in snp_snps])
snp_weight = np.where(snp_is_focal, Wss, w_bg)

# Focal fixed effects (on raw 0/1/2 scale)
focal_eff_vec = np.zeros(M)
for i, s in enumerate(snp_snps):
    if s in focal_ids:
        focal_eff_vec[i] = focal_signs[s] * np.sqrt(Wss * sigma2_alpha) / np.sqrt(snp_pq2[i])

# --- Focal SNP genotype CSV (intercept + dosages for REML covariate) ---
raw = pd.read_csv(f"{outdir}/focal_geno.raw", sep=r'\s+', engine='python')
geno_out = pd.DataFrame({'IID': raw['IID'], 'intercept': 1})
snp_cols = [c for c in raw.columns if c not in ['FID','IID','PAT','MAT','SEX','PHENOTYPE']]
for col in snp_cols:
    clean = col.rsplit('_', 1)[0]
    geno_out[clean] = raw[col]
geno_out.to_csv(f"{outdir}/focal_snp_genotypes.csv", index=False)
print(f"  Focal genotype CSV: {len(snp_cols)} SNPs")

# --- Simulate replicates ---
np.random.seed(42)
pheno_cols = {}
maf_str = [f"{m:.6g}" for m in snp_maf]

for rep in range(1, nrep + 1):
    rep_dir = f"{outdir}/rep_{rep}"
    os.makedirs(rep_dir, exist_ok=True)

    # Background effects: N(0, w_bg) / sqrt(2pq); focal: fixed
    bg_alpha = np.random.normal(0, np.sqrt(w_bg * sigma2_alpha), size=M)
    effects = np.where(snp_is_focal, focal_eff_vec, bg_alpha / np.sqrt(snp_pq2))

    # Write SLEMM SNP effect CSV
    out_eff = pd.DataFrame({
        'SNP': snp_snps, 'Chr': snp_chrs, 'Pos': snp_pos,
        'Allele1': snp_a1, 'Allele2': snp_a2,
        'MAF': maf_str, 'HWE_Pval': 1, 'Group': 'NULL',
        'Weight': snp_weight, 'Effect': effects
    })
    out_eff.to_csv(f"{rep_dir}/sim.snp.csv", index=False)

    # Compute GEBV via SLEMM
    os.system(
        f"slemm --pred --bfile {outdir}/sub_{N} "
        f"--snp_estimate {rep_dir}/sim.snp.csv "
        f"--out {rep_dir}/sim.gv.csv 2>/dev/null"
    )

    # Read GV, add residual
    gv = pd.read_csv(f"{rep_dir}/sim.gv.csv")
    resid = np.random.normal(0, np.sqrt(sigma2_e), size=len(gv))
    pheno_cols[f"pheno_{rep}"] = gv['GEBV'].values + resid

    if rep % 20 == 0 or rep == 1 or rep == nrep:
        print(f"  Rep {rep}/{nrep} done")

# --- Write phenotype file ---
pheno_df = pd.concat([fam[['FID','IID']], pd.DataFrame(pheno_cols)], axis=1)
pheno_df.to_csv(f"{outdir}/sim.pheno.txt", sep=' ', index=False)
pheno_df.drop(columns=['FID']).to_csv(f"{outdir}/sim.pheno.csv", index=False)

print(f"\n  Phenotypes: {outdir}/sim.pheno.txt ({nrep} columns)")
print("=== sim_step1 complete ===")
PYSIM
