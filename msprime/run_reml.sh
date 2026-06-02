#!/bin/bash
# run_reml.sh — Run SLEMM REML for all replicates in a simulation directory
#
# Usage:
#   bash run_reml.sh <simdir> [<num_threads>]
#
# Example:
#   bash run_reml.sh ./plink_out/Ne20/h2_0.3__beta2_0.01/N_5000 10
#
# Expects in <simdir>/:
#   sub_<N>.bed/bim/fam       — Genotypes
#   sim.pheno.csv             — IID, pheno_1, ..., pheno_<nrep>
#   focal_snp_genotypes.csv   — Focal SNP covariate file
#   snp_info.csv              — SNP IDs for GRM (single column, header "SNP")
#
# Runs SLEMM --reml for each pheno_<r> trait, with focal SNPs as covariates.

set -euo pipefail

simdir=$1
nthreads=${2:-1}

# --- Auto-detect N from bed file ---
bedfile=$(ls ${simdir}/sub_*.bed 2>/dev/null | head -1)
if [ -z "${bedfile}" ]; then
    echo "ERROR: No sub_*.bed found in ${simdir}"
    exit 1
fi
bprefix="${bedfile%.bed}"

# --- Auto-detect number of replicates from pheno CSV header ---
nrep=$(head -1 ${simdir}/sim.pheno.csv | tr ',' '\n' | grep -c 'pheno_')

# --- Extract h2 from directory name and set max_herit ---
h2=$(echo "${simdir}" | grep -oP 'h2_\K[0-9.]+')
if [ -z "${h2}" ]; then
    echo "ERROR: Cannot parse h2 from directory name: ${simdir}"
    exit 1
fi
max_h2=$(python3 -c "print(min(${h2} + 0.1, 0.99))")

# --- Generate snp_info.csv if not present ---
if [ ! -f "${simdir}/snp_info.csv" ]; then
    echo "Generating snp_info.csv..."
    echo "SNP" > ${simdir}/snp_info.csv
    awk '{print $2}' ${bprefix}.bim >> ${simdir}/snp_info.csv
fi

echo "=== run_reml: ${simdir} ==="
echo "  bfile: ${bprefix}"
echo "  replicates: ${nrep}"
echo "  threads: ${nthreads}"
echo "  h2=${h2}, max_herit=${max_h2}"

for rep in $(seq 1 ${nrep}); do
    mkdir -p ${simdir}/rep_${rep}

    slemm --reml --lrt \
        --bfile ${bprefix} \
        --phenotype_file ${simdir}/sim.pheno.csv \
        --trait pheno_${rep} \
        --covariate_file ${simdir}/focal_snp_genotypes.csv \
        --covariate_names all \
        --snp_info_file ${simdir}/snp_info.csv \
        --max_herit ${max_h2} \
        --out ${simdir}/rep_${rep}/reml \
        --num_threads ${nthreads} \
        2>/dev/null

    if [ $((rep % 20)) -eq 0 ] || [ ${rep} -eq 1 ] || [ ${rep} -eq ${nrep} ]; then
        echo "  Rep ${rep}/${nrep} done"
    fi
done

echo "=== run_reml complete ==="

