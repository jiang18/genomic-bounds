#!/usr/bin/env python3
"""
simulate_genotypes.py

Simulate genotypes for theoretical GWAS validation using msprime.

Key design choices:
  - Two-epoch demography: recent Ne (target) + large ancestral Ne
    Ensures variants are old enough to be embedded in LD structure.
  - hudson coalescent: Ne is a drift-rate parameter, not a census size.
    Correct for arbitrary sample sizes.
  - BinaryMutationModel: all sites biallelic by construction; no post-hoc filter.
  - MAF > 0.05 filter: mimics chip ascertainment (common variants only).
  - Random thinning to target M: avoids artificial SNP spacing regularity.
  - Returns numpy array directly: no PLINK dependency needed for eigenanalysis.

Usage:
  python3 simulate_genotypes.py <Ne> <N_diploid> <M_total> <outdir>

Example:
  python3 simulate_genotypes.py 200 10000 50000 sim_out/Ne200_N10k
"""

import os
import sys
import numpy as np

try:
    import msprime
except ImportError:
    print("ERROR: pip install msprime")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Cattle chromosome lengths (ARS-UCD1.2 autosomes, 29 chromosomes)
# ---------------------------------------------------------------------------
CHR_LENGTHS_MB = [
    158.5, 137.1, 121.4, 120.7, 121.2,
    119.5, 112.6, 113.4, 105.7, 104.3,
    107.3,  87.2,  84.2,  83.5,  85.3,
     81.7,  75.2,  66.0,  63.7,  72.0,
     69.9,  61.4,  52.5,  62.7,  42.9,
     51.7,  45.4,  46.3,  51.5,
]
CHR_LENGTHS_BP = [int(round(x * 1e6)) for x in CHR_LENGTHS_MB]
TOTAL_LENGTH_BP = sum(CHR_LENGTHS_BP)
N_CHR = len(CHR_LENGTHS_BP)


def make_demography(Ne, Ne_ancestral=10_000, t_transition=None):
    """
    Two-epoch demographic model.

      0 ... t_transition generations ago : Ne          (recent, governs chip LD)
      t_transition ... inf               : Ne_ancestral (ensures old variants)

    Transition time defaults to 10 * Ne generations, which is long enough
    that most lineages have not yet coalesced in the recent epoch, so the
    LD structure at chip scale is determined by Ne, not Ne_ancestral.
    """
    if t_transition is None:
        t_transition = 10 * Ne
    dem = msprime.Demography()
    dem.add_population(name="A", initial_size=Ne)
    dem.add_population_parameters_change(
        time=t_transition,
        initial_size=Ne_ancestral,
        population="A",
    )
    return dem


def simulate_chromosome(chrom, L_bp, Ne, N_diploid, target_snps,
                         maf_min, recomb_rate, mutation_rate, seed):
    """
    Simulate one chromosome and return a standardised genotype matrix.

    Returns
    -------
    G : np.ndarray, shape (N_diploid, n_snps_retained)
        Standardised genotype matrix (mean 0, std 1 per SNP).
    positions : np.ndarray, shape (n_snps_retained,)
        Physical positions (bp) of retained SNPs.
    mafs : np.ndarray, shape (n_snps_retained,)
        Minor allele frequencies of retained SNPs.
    """
    dem = make_demography(Ne)

    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(N_diploid, ploidy=2)],
        demography=dem,
        sequence_length=L_bp,
        recombination_rate=recomb_rate,
        model="hudson",
        random_seed=seed,
    )
    # BinaryMutationModel: all sites biallelic by construction
    ts = msprime.sim_mutations(
        ts,
        rate=mutation_rate,
        model=msprime.BinaryMutationModel(),
        random_seed=seed + 100_000,
    )

    # Stream variants: all sites biallelic, only MAF filter needed
    positions_all, dosages_all, mafs_all = [], [], []
    for var in ts.variants():
        g_hap = var.genotypes                        # shape (2*N_diploid,)
        freq = g_hap.mean()
        maf = min(freq, 1.0 - freq)
        if maf < maf_min:
            continue
        dosage = g_hap[0::2] + g_hap[1::2]          # diploid dosage, shape (N_diploid,)
        positions_all.append(var.site.position)
        dosages_all.append(dosage)
        mafs_all.append(maf)

    n_pass = len(positions_all)
    if n_pass < target_snps:
        raise RuntimeError(
            f"chr{chrom}: only {n_pass} SNPs pass MAF>{maf_min}; "
            f"need {target_snps}. Increase mutation_rate."
        )

    # Random thinning to target_snps
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(n_pass, size=target_snps, replace=False))

    dosages = np.stack([dosages_all[i] for i in idx], axis=1)   # (N_diploid, target_snps)
    positions = np.array([positions_all[i] for i in idx])
    mafs = np.array([mafs_all[i] for i in idx])

    # Standardise
    mu = dosages.mean(axis=0)
    sd = dosages.std(axis=0)
    sd[sd == 0] = 1.0
    G = (dosages - mu) / sd

    return G, positions, mafs


def simulate_genome(Ne, N_diploid, M_total,
                    maf_min=0.05,
                    recomb_rate=1e-8,
                    mutation_rate=5e-8,
                    base_seed=42):
    """
    Simulate all 29 cattle autosomes and concatenate.

    Returns
    -------
    G : np.ndarray, shape (N_diploid, M_total)
    chrom_ids : np.ndarray, shape (M_total,)   (1-based chromosome index)
    positions : np.ndarray, shape (M_total,)   (bp within chromosome)
    mafs : np.ndarray, shape (M_total,)
    """
    # Distribute SNPs proportional to chromosome length
    chr_targets = [int(round(M_total * L / TOTAL_LENGTH_BP)) for L in CHR_LENGTHS_BP]
    diff = M_total - sum(chr_targets)
    for i in range(abs(diff)):
        chr_targets[i % N_CHR] += 1 if diff > 0 else -1
    assert sum(chr_targets) == M_total

    G_list, chrom_id_list, pos_list, maf_list = [], [], [], []

    for c in range(N_CHR):
        chrom = c + 1
        L = CHR_LENGTHS_BP[c]
        target = chr_targets[c]
        print(f"  chr{chrom:2d}: {CHR_LENGTHS_MB[c]:.1f} Mb, target {target} SNPs", flush=True)

        G_chr, pos_chr, maf_chr = simulate_chromosome(
            chrom=chrom,
            L_bp=L,
            Ne=Ne,
            N_diploid=N_diploid,
            target_snps=target,
            maf_min=maf_min,
            recomb_rate=recomb_rate,
            mutation_rate=mutation_rate,
            seed=base_seed * 100 + chrom,
        )
        G_list.append(G_chr)
        chrom_id_list.append(np.full(target, chrom, dtype=np.int32))
        pos_list.append(pos_chr)
        maf_list.append(maf_chr)

    G = np.concatenate(G_list, axis=1)
    chrom_ids = np.concatenate(chrom_id_list)
    positions  = np.concatenate(pos_list)
    mafs       = np.concatenate(maf_list)
    return G, chrom_ids, positions, mafs


def main():
    if len(sys.argv) != 5:
        print(__doc__.strip())
        sys.exit(2)

    Ne        = int(sys.argv[1])
    N_diploid = int(sys.argv[2])
    M_total   = int(sys.argv[3])
    outdir    = sys.argv[4]
    os.makedirs(outdir, exist_ok=True)

    print(f"Ne={Ne}, N={N_diploid}, M={M_total}")

    G, chrom_ids, positions, mafs = simulate_genome(
        Ne=Ne,
        N_diploid=N_diploid,
        M_total=M_total,
    )

    # Save
    np.save(os.path.join(outdir, "G.npy"),         G)
    np.save(os.path.join(outdir, "chrom_ids.npy"), chrom_ids)
    np.save(os.path.join(outdir, "positions.npy"), positions)
    np.save(os.path.join(outdir, "mafs.npy"),      mafs)

    print(f"\nSaved to {outdir}/")
    print(f"  G shape : {G.shape}")
    print(f"  MAF : min={mafs.min():.3f}  median={np.median(mafs):.3f}  max={mafs.max():.3f}")


if __name__ == "__main__":
    main()

