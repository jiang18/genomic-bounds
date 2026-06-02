#!/usr/bin/env python3
"""Compute exact per-SNP LOCO GRAMMAR-Gamma coefficients from PLINK BED sets."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np


BED_MAGIC = b"\x6c\x1b\x01"


def build_decode_table() -> np.ndarray:
    table = np.empty((256, 4), dtype=np.float64)
    mapping = {
        0b00: 0.0,
        0b01: np.nan,
        0b10: 1.0,
        0b11: 2.0,
    }
    for byte in range(256):
        for shift in range(4):
            code = (byte >> (2 * shift)) & 0b11
            table[byte, shift] = mapping[code]
    return table


DECODE_TABLE = build_decode_table()


@dataclass
class BimRecord:
    chrom: str
    snp: str
    cm: str
    bp: str
    a1: str
    a2: str


class BedReader:
    def __init__(self, prefix: Path):
        self.prefix = prefix
        self.bed_path = prefix.with_suffix(".bed")
        self.bim_path = prefix.with_suffix(".bim")
        self.fam_path = prefix.with_suffix(".fam")
        self.n_samples = sum(1 for _ in self.fam_path.open())
        with self.bim_path.open() as fh:
            self.bim = [BimRecord(*line.split()[:6]) for line in fh]
        self.n_variants = len(self.bim)
        self.bytes_per_variant = (self.n_samples + 3) // 4
        with self.bed_path.open("rb") as fh:
            magic = fh.read(3)
        if magic != BED_MAGIC:
            raise ValueError(f"{self.bed_path} is not a PLINK SNP-major BED file")

    def read_chunk(self, start: int, count: int) -> np.ndarray:
        if count <= 0:
            return np.empty((self.n_samples, 0), dtype=np.float64)
        with self.bed_path.open("rb") as fh:
            fh.seek(3 + start * self.bytes_per_variant)
            raw = fh.read(count * self.bytes_per_variant)
        arr = np.frombuffer(raw, dtype=np.uint8).reshape(count, self.bytes_per_variant)
        decoded = DECODE_TABLE[arr].reshape(count, self.bytes_per_variant * 4)[:, : self.n_samples]
        return decoded.T.copy()


def standardize_chunk(geno: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    means = np.nanmean(geno, axis=0)
    p = means / 2.0
    scale = np.sqrt(2.0 * p * (1.0 - p))
    centered = geno - means
    mask = np.isnan(centered)
    if mask.any():
        centered[mask] = 0.0
    centered /= scale
    return centered, means, scale


def build_crossproducts(reader: BedReader, chunk_size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = reader.n_samples
    m = reader.n_variants
    total = np.zeros((n, n), dtype=np.float64)
    means = np.empty(m, dtype=np.float64)
    scales = np.empty(m, dtype=np.float64)
    for start in range(0, m, chunk_size):
        count = min(chunk_size, m - start)
        geno = reader.read_chunk(start, count)
        z, chunk_means, chunk_scales = standardize_chunk(geno)
        means[start : start + count] = chunk_means
        scales[start : start + count] = chunk_scales
        total += z @ z.T
    return total, means, scales


def chrom_to_indices(reader: BedReader) -> dict[str, list[int]]:
    out: dict[str, list[int]] = defaultdict(list)
    for idx, rec in enumerate(reader.bim):
        out[rec.chrom].append(idx)
    return dict(out)


def read_standardized_indices(
    reader: BedReader,
    indices: list[int],
    means: np.ndarray,
    scales: np.ndarray,
    chunk_size: int,
) -> np.ndarray:
    z = np.empty((reader.n_samples, len(indices)), dtype=np.float64)
    pos = 0
    for block_start in range(0, len(indices), chunk_size):
        block_indices = indices[block_start : block_start + chunk_size]
        first = block_indices[0]
        if block_indices[-1] - first + 1 != len(block_indices):
            raise ValueError("Chromosome SNPs are not contiguous in BIM order")
        geno = reader.read_chunk(first, len(block_indices))
        block = geno - means[first : first + len(block_indices)]
        mask = np.isnan(block)
        if mask.any():
            block[mask] = 0.0
        block /= scales[first : first + len(block_indices)]
        z[:, pos : pos + len(block_indices)] = block
        pos += len(block_indices)
    return z


def compute_loco_coefficients(
    reader: BedReader,
    total_crossprod: np.ndarray,
    means: np.ndarray,
    scales: np.ndarray,
    h2_values: list[float],
    chunk_size: int,
) -> tuple[dict[float, np.ndarray], dict[str, dict[str, float]]]:
    n = reader.n_samples
    m = reader.n_variants
    coeffs = {h2: np.empty(m, dtype=np.float64) for h2 in h2_values}
    chrom_meta: dict[str, dict[str, float]] = {}
    chrom_indices = chrom_to_indices(reader)
    for chrom, indices in sorted(chrom_indices.items(), key=lambda kv: int(kv[0])):
        z_chr = read_standardized_indices(reader, indices, means, scales, chunk_size)
        chr_crossprod = z_chr @ z_chr.T
        m_other = m - len(indices)
        if m_other <= 0:
            raise ValueError(f"Chromosome {chrom} leaves no SNPs for LOCO GRM")
        g_other = (total_crossprod - chr_crossprod) / m_other
        eigenvalues, eigenvectors = np.linalg.eigh(g_other)
        y = eigenvectors.T @ z_chr
        chrom_meta[chrom] = {
            "n_snps_chr": len(indices),
            "n_snps_other": m_other,
            "grm_eig_min": float(eigenvalues.min()),
            "grm_eig_max": float(eigenvalues.max()),
        }
        y_sq = y * y
        for h2 in h2_values:
            denom = h2 * eigenvalues + (1.0 - h2)
            coeffs[h2][indices] = np.sum(y_sq / denom[:, None], axis=0) / n
    return coeffs, chrom_meta


def write_results(
    reader: BedReader,
    coeffs: dict[float, np.ndarray],
    output_tsv: Path,
    meta_json: Path,
    h2_values: list[float],
    chrom_meta: dict[str, dict[str, float]],
) -> None:
    with output_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        header = ["chrom", "snp", "cm", "bp", "a1", "a2"] + [f"c_loco_h2_{h2:g}" for h2 in h2_values]
        writer.writerow(header)
        for idx, rec in enumerate(reader.bim):
            row = [rec.chrom, rec.snp, rec.cm, rec.bp, rec.a1, rec.a2]
            row.extend(coeffs[h2][idx] for h2 in h2_values)
            writer.writerow(row)
    summary = {
        "n_samples": reader.n_samples,
        "n_snps": reader.n_variants,
        "n_chromosomes": len(chrom_meta),
        "h2_values": h2_values,
        "mean_coefficients": {str(h2): float(np.mean(coeffs[h2])) for h2 in h2_values},
        "sd_coefficients": {str(h2): float(np.std(coeffs[h2])) for h2 in h2_values},
        "min_coefficients": {str(h2): float(np.min(coeffs[h2])) for h2 in h2_values},
        "max_coefficients": {str(h2): float(np.max(coeffs[h2])) for h2 in h2_values},
        "chromosomes": chrom_meta,
    }
    meta_json.write_text(json.dumps(summary, indent=2))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bfile", required=True, type=Path, help="PLINK prefix without extension")
    parser.add_argument("--output-tsv", required=True, type=Path)
    parser.add_argument("--meta-json", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--h2", nargs="+", type=float, default=[0.1, 0.3, 0.5])
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    reader = BedReader(args.bfile)
    total_crossprod, means, scales = build_crossproducts(reader, args.chunk_size)
    coeffs, chrom_meta = compute_loco_coefficients(
        reader=reader,
        total_crossprod=total_crossprod,
        means=means,
        scales=scales,
        h2_values=args.h2,
        chunk_size=args.chunk_size,
    )
    write_results(reader, coeffs, args.output_tsv, args.meta_json, args.h2, chrom_meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
