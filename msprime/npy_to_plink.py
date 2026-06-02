#!/usr/bin/env python3
"""
Convert .npy genotype data (from msprime simulations) to PLINK bed/bim/fam format.

Usage:
    python npy_to_plink.py sim_out/Ne20 plink_out/Ne20

Input files (in input_dir):
    G.npy          - N x M genotype matrix (0/1/2 dosage)
    positions.npy  - M-length array of base-pair positions
    chrom_ids.npy  - M-length array of chromosome IDs (0-indexed or 1-indexed)
    mafs.npy       - M-length array of minor allele frequencies

Output files (in output_dir):
    genotypes.bed  - PLINK binary genotype file
    genotypes.bim  - variant information
    genotypes.fam  - sample information
"""

import argparse
import numpy as np
import os


def write_bed(G, filepath):
    """
    Write PLINK .bed file from N x M genotype matrix (0/1/2 dosages).
    Vectorized implementation for large N.

    PLINK .bed format (SNP-major mode):
      - 3-byte magic header: 0x6C 0x1B 0x01
      - For each SNP: ceil(N/4) bytes, 4 samples per byte (LSB first)
      - Encoding per 2-bit pair:
          00 = homozygous A1 (minor) -> dosage 2
          01 = missing
          10 = heterozygous           -> dosage 1
          11 = homozygous A2 (major)  -> dosage 0
    """
    N, M = G.shape
    bytes_per_snp = (N + 3) // 4

    # Map dosage -> 2-bit code
    # dosage 0 -> 0b11 = 3, dosage 1 -> 0b10 = 2, dosage 2 -> 0b00 = 0
    code_map = np.array([3, 2, 0, 1], dtype=np.uint8)  # index 3 = missing fallback

    # Pad N to multiple of 4
    n_pad = bytes_per_snp * 4
    G_pad = np.zeros((n_pad, M), dtype=np.int8)
    G_pad[:N, :] = G
    # Padded samples get dosage 0 -> code 3 (0b11), which is the PLINK convention
    # for padding bits (they are ignored by PLINK)

    # Map all genotypes to 2-bit codes
    codes = code_map[np.clip(G_pad, 0, 3).astype(np.uint8)]  # (n_pad, M)

    # Reshape to (bytes_per_snp, 4, M) and pack 4 codes per byte
    codes = codes.reshape(bytes_per_snp, 4, M)
    packed = (codes[:, 0, :]
              | (codes[:, 1, :] << 2)
              | (codes[:, 2, :] << 4)
              | (codes[:, 3, :] << 6))  # (bytes_per_snp, M), uint8

    with open(filepath, "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        # Write SNP-major: iterate over SNPs (columns)
        for j in range(M):
            f.write(packed[:, j].tobytes())
            if (j + 1) % 5000 == 0:
                print(f"    {j+1}/{M} SNPs written...")


def write_bim(positions, chrom_ids, mafs, filepath):
    """
    Write PLINK .bim file.
    Columns: chr, snp_id, genetic_dist(cM), bp_pos, allele1(minor), allele2(major)
    """
    M = len(positions)
    with open(filepath, "w") as f:
        for j in range(M):
            chrom = int(chrom_ids[j])
            bp = int(positions[j])
            snp_id = f"snp_{chrom}_{bp}"
            # genetic distance: 0 placeholder (not needed for GWAS)
            # alleles: A = minor, B = major (arbitrary but consistent)
            f.write(f"{chrom}\t{snp_id}\t0\t{bp}\tA\tB\n")


def write_fam(N, filepath):
    """
    Write PLINK .fam file.
    Columns: FID, IID, father, mother, sex, phenotype
    """
    with open(filepath, "w") as f:
        for i in range(N):
            fid = f"FAM{i+1}"
            iid = f"IND{i+1}"
            f.write(f"{fid}\t{iid}\t0\t0\t0\t-9\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert .npy genotype data to PLINK bed/bim/fam"
    )
    parser.add_argument("input_dir", help="Directory with G.npy, positions.npy, etc.")
    parser.add_argument("output_dir", help="Output directory for PLINK files")
    parser.add_argument(
        "--prefix", default="genotypes", help="Output file prefix (default: genotypes)"
    )
    parser.add_argument(
        "--chrom-offset", type=int, default=0,
        help="Add this to chrom_ids values (default: 0, i.e. already 1-indexed)"
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load data
    print(f"Loading data from {args.input_dir}...")
    G = np.load(os.path.join(args.input_dir, "G.npy"))
    positions = np.load(os.path.join(args.input_dir, "positions.npy"))
    chrom_ids = np.load(os.path.join(args.input_dir, "chrom_ids.npy"))
    mafs = np.load(os.path.join(args.input_dir, "mafs.npy"))

    N, M = G.shape
    print(f"  N={N} samples, M={M} SNPs")
    print(f"  Chromosomes: {np.unique(chrom_ids)}")
    print(f"  MAF range: [{mafs.min():.4f}, {mafs.max():.4f}]")

    # Validate dimensions
    assert positions.shape[0] == M, f"positions length {positions.shape[0]} != M={M}"
    assert chrom_ids.shape[0] == M, f"chrom_ids length {chrom_ids.shape[0]} != M={M}"
    assert mafs.shape[0] == M, f"mafs length {mafs.shape[0]} != M={M}"

    # Detect if G is standardized (values outside [0,2]) and unstandardize
    if G.min() < -0.5 or G.max() > 2.5:
        print("  G appears standardized — converting back to 0/1/2 dosages...")
        # z = (g - 2p) / sqrt(2p(1-p))  =>  z*sd = g - 2p
        # So z*sd takes values {-2p, 1-2p, 2-2p} for dosages {0, 1, 2}.
        # With N=100k and MAF>=0.05, HWE guarantees all 3 genotypes present.
        # Recover p per-SNP: min(z*sd) = -2p  =>  p = -min/2.
        p = mafs.copy()
        sd = np.sqrt(2.0 * p * (1.0 - p))
        # Multiply in-place if G is float; otherwise make one copy
        if np.issubdtype(G.dtype, np.floating):
            G *= sd[np.newaxis, :]
            G_shifted = G  # just a rename, no copy
        else:
            G_shifted = G.astype(np.float64)
            G_shifted *= sd[np.newaxis, :]
            del G

        p_recovered = np.empty(M, dtype=np.float64)
        n_not3 = 0
        for j in range(M):
            col = G_shifted[:, j]
            vals = np.unique(np.round(col, 6))
            if len(vals) != 3:
                n_not3 += 1
                # Fallback: use MAF (should not happen under HWE with large N)
                p_recovered[j] = mafs[j]
            else:
                p_recovered[j] = -vals[0] / 2.0

        if n_not3 > 0:
            print(f"  WARNING: {n_not3}/{M} SNPs do not have 3 distinct genotypes"
                  f" — unexpected under HWE with large N and MAF>=0.05")
        else:
            print(f"  All {M} SNPs have 3 distinct genotypes (HWE sanity check passed)")

        # Unstandardize in-place to save memory:
        # G_shifted currently holds z*sd = g - 2p; add 2p, round, clip
        G_shifted += 2.0 * p_recovered[np.newaxis, :]
        np.rint(G_shifted, out=G_shifted)
        np.clip(G_shifted, 0, 2, out=G_shifted)
        G = G_shifted.astype(np.int8)  # int8 is 8x smaller than float64
        del G_shifted

        # Sanity check
        maf_check = np.minimum(G.mean(axis=0) / 2.0, 1.0 - G.mean(axis=0) / 2.0)
        maf_err = np.abs(maf_check - mafs).max()
        print(f"  Max MAF reconstruction error: {maf_err:.6f}")
        if maf_err > 0.02:
            print("  WARNING: large MAF reconstruction error — check input data")
        n_flipped = np.sum(np.abs(p_recovered - mafs) > 0.01)
        print(f"  SNPs with p != MAF (allele freq > 0.5): {n_flipped}/{M}")
    else:
        G = np.rint(G).astype(np.int8)

    # Adjust chromosome indexing
    chrom_ids_adj = chrom_ids + args.chrom_offset

    prefix = os.path.join(args.output_dir, args.prefix)

    # Write .fam
    fam_path = prefix + ".fam"
    print(f"Writing {fam_path}...")
    write_fam(N, fam_path)

    # Write .bim
    bim_path = prefix + ".bim"
    print(f"Writing {bim_path}...")
    write_bim(positions, chrom_ids_adj, mafs, bim_path)

    # Write .bed
    bed_path = prefix + ".bed"
    print(f"Writing {bed_path} (this may take a moment)...")
    write_bed(G, bed_path)

    print(f"Done. PLINK files written to {args.output_dir}/{args.prefix}.*")
    print(f"  Verify with: plink --bfile {prefix} --freq --out {prefix}_check")


if __name__ == "__main__":
    main()

