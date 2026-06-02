#!/bin/bash
# select_focal_snps.sh — Select 5 focal SNPs from the full dataset (run once)
#
# Usage:
#   bash select_focal_snps.sh ./plink_out/Ne20
#
# Computes MAF from the full dataset, selects 5 SNPs on chr1-5 with
# specified MAF bins, and writes results to the same directory.
# These files are reused by all subsequent simulation runs.
#
# Output:
#   <dir>/full.frq              — MAF from full dataset
#   <dir>/focal_snps.txt        — SNP IDs (one per line)
#   <dir>/focal_snps_info.csv   — SNP details (chr, MAF, bin)

set -euo pipefail

dir=$1
bed_prefix=${dir}/genotypes

echo "=== Selecting focal SNPs from full dataset ==="

# Compute MAF from full data
echo "[1/2] Computing MAF..."
plink --chr-set 30 --bfile ${bed_prefix} --freq --out ${dir}/full --allow-no-sex --silent

# Select focal SNPs
echo "[2/2] Selecting focal SNPs..."

python3 - "${dir}" << 'PYSELECT'
import pandas as pd
import sys

dir = sys.argv[1]
frq = pd.read_csv(f"{dir}/full.frq", sep=r'\s+', engine='python')

bins = [
    (1, 0.05, 0.10),
    (2, 0.15, 0.20),
    (3, 0.25, 0.30),
    (4, 0.35, 0.40),
    (5, 0.45, 0.50),
]

focal = []
for chrom, maf_lo, maf_hi in bins:
    candidates = frq[(frq['CHR'] == chrom) & (frq['MAF'] >= maf_lo) & (frq['MAF'] <= maf_hi)]
    if len(candidates) == 0:
        print(f"  WARNING: no SNP found for chr{chrom} MAF [{maf_lo},{maf_hi}]")
        continue
    mid = (maf_lo + maf_hi) / 2
    best = candidates.iloc[(candidates['MAF'] - mid).abs().argsort().iloc[0]]
    focal.append({'SNP': best['SNP'], 'CHR': int(best['CHR']),
                  'MAF': best['MAF'], 'MAF_bin': f"{maf_lo}-{maf_hi}"})
    print(f"  chr{chrom} MAF [{maf_lo},{maf_hi}]: {best['SNP']} (MAF={best['MAF']:.4f})")

focal_df = pd.DataFrame(focal)
focal_df.to_csv(f"{dir}/focal_snps_info.csv", index=False)

with open(f"{dir}/focal_snps.txt", 'w') as f:
    for sid in focal_df['SNP']:
        f.write(sid + '\n')

print(f"\n  Selected {len(focal_df)} focal SNPs -> {dir}/focal_snps.txt")
PYSELECT

echo "=== Done ==="

