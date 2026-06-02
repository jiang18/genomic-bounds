#!/bin/bash
# run_all_scenarios.sh — Run phenotype simulation + REML + collect for all parameter combinations.
#
# Usage:
#   bash run_all_scenarios.sh <Ne> <nrep> <ncores>
#
# Example:
#   bash run_all_scenarios.sh Ne20 100 28
#
# Prerequisites (already done):
#   - msprime simulation -> sim_out/<Ne>/
#   - npy_to_plink.py -> plink_out/<Ne>/
#   - select_focal_snps.sh -> plink_out/<Ne>/focal_snps_info.csv
#   - compute_eigenvalues_by_N.py -> eigen_out/<Ne>/
#   - project_focal_snps.py -> eigen_out/<Ne>/focal_proj_N*.npz

set -e

if [ $# -lt 3 ]; then
    echo "Usage: bash run_all_scenarios.sh <Ne> <nrep> <ncores>"
    echo "Example: bash run_all_scenarios.sh Ne20 100 28"
    exit 1
fi

NE=$1
NREP=$2
NCORES=$3

PLINK_DIR=./plink_out/${NE}
EIGEN_DIR=./eigen_out/${NE}

# Parameter grid
N_VALUES="5000 10000 20000 50000"
H2_VALUES="0.1 0.3 0.5"
BETA2_VALUES="0.001 0.01"

# Count total scenarios
n_total=0
for N in $N_VALUES; do
for h2 in $H2_VALUES; do
for beta2 in $BETA2_VALUES; do
    n_total=$((n_total + 1))
done; done; done
echo "=== run_all_scenarios: ${NE}, ${NREP} reps, ${NCORES} cores ==="
echo "  N values:    ${N_VALUES}"
echo "  h2 values:   ${H2_VALUES}"
echo "  beta2 values: ${BETA2_VALUES}"
echo "  Total scenarios: ${n_total}"
echo ""

# Track progress
n_done=0

for N in $N_VALUES; do
for h2 in $H2_VALUES; do
for beta2 in $BETA2_VALUES; do
    n_done=$((n_done + 1))
    SIMDIR=${PLINK_DIR}/h2_${h2}__beta2_${beta2}/N_${N}

    echo "============================================================"
    echo "[${n_done}/${n_total}] N=${N} h2=${h2} beta2=${beta2}"
    echo "============================================================"

    # Step 1: Simulate phenotypes
    if [ ! -f "${SIMDIR}/sim.pheno.txt" ]; then
        echo "  -> sim_step1.sh"
        bash sim_step1.sh ${N} ${h2} ${beta2} ${NREP} ${PLINK_DIR} ${EIGEN_DIR}
    else
        echo "  -> sim_step1.sh: SKIP (sim.pheno.txt exists)"
    fi

    # Step 2: Run REML
    if [ ! -f "${SIMDIR}/rep_${NREP}/reml.reml.blue.csv" ]; then
        echo "  -> run_reml.sh"
        bash run_reml.sh ${SIMDIR} ${NCORES}
    else
        echo "  -> run_reml.sh: SKIP (all reps done)"
    fi

    # Step 3: Collect results
    # h2 for GRM = h2_total - 5*beta2 (focal SNPs are fixed effects, not in GRM)
    H2_GRM=$(awk "BEGIN {printf \"%.6f\", ${h2} - 5*${beta2}}")
    
    echo "  -> collect_results.py (h2_grm=${H2_GRM})"
    python3 collect_results.py ${SIMDIR} ${EIGEN_DIR} ${H2_GRM} ${beta2}

    echo ""
done; done; done

echo "=== All ${n_total} scenarios complete ==="

