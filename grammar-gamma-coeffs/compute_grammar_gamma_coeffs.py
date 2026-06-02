#!/usr/bin/env python3
"""Compute exact per-SNP GRAMMAR-Gamma coefficients from PLINK BED sets."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve


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
        """Return genotype dosage matrix of shape (n_samples, count)."""
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


def build_grm(reader: BedReader, chunk_size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = reader.n_samples
    m = reader.n_variants
    grm = np.zeros((n, n), dtype=np.float64)
    means = np.empty(m, dtype=np.float64)
    scales = np.empty(m, dtype=np.float64)
    for start in range(0, m, chunk_size):
        count = min(chunk_size, m - start)
        geno = reader.read_chunk(start, count)
        z, chunk_means, chunk_scales = standardize_chunk(geno)
        means[start : start + count] = chunk_means
        scales[start : start + count] = chunk_scales
        grm += z @ z.T
    grm /= m
    return grm, means, scales


def compute_coefficients(
    reader: BedReader,
    means: np.ndarray,
    scales: np.ndarray,
    grm: np.ndarray,
    h2_values: list[float],
    chunk_size: int,
) -> dict[float, np.ndarray]:
    n = reader.n_samples
    identity = np.eye(n, dtype=np.float64)
    factors = {}
    for h2 in h2_values:
        v = grm * h2 + identity * (1.0 - h2)
        factors[h2] = cho_factor(v, lower=True, check_finite=False)
    coeffs = {h2: np.empty(reader.n_variants, dtype=np.float64) for h2 in h2_values}
    for start in range(0, reader.n_variants, chunk_size):
        count = min(chunk_size, reader.n_variants - start)
        geno = reader.read_chunk(start, count)
        z = geno - means[start : start + count]
        mask = np.isnan(z)
        if mask.any():
            z[mask] = 0.0
        z /= scales[start : start + count]
        for h2 in h2_values:
            y = cho_solve(factors[h2], z, check_finite=False)
            coeffs[h2][start : start + count] = np.sum(z * y, axis=0) / n
    return coeffs


def write_results(
    reader: BedReader,
    coeffs: dict[float, np.ndarray],
    output_tsv: Path,
    meta_json: Path,
    h2_values: list[float],
) -> None:
    with output_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        header = ["chrom", "snp", "cm", "bp", "a1", "a2"] + [f"c_h2_{h2:g}" for h2 in h2_values]
        writer.writerow(header)
        for idx, rec in enumerate(reader.bim):
            row = [rec.chrom, rec.snp, rec.cm, rec.bp, rec.a1, rec.a2]
            row.extend(coeffs[h2][idx] for h2 in h2_values)
            writer.writerow(row)
    summary = {
        "n_samples": reader.n_samples,
        "n_snps": reader.n_variants,
        "h2_values": h2_values,
        "mean_coefficients": {str(h2): float(np.mean(coeffs[h2])) for h2 in h2_values},
        "sd_coefficients": {str(h2): float(np.std(coeffs[h2])) for h2 in h2_values},
        "min_coefficients": {str(h2): float(np.min(coeffs[h2])) for h2 in h2_values},
        "max_coefficients": {str(h2): float(np.max(coeffs[h2])) for h2 in h2_values},
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
    grm, means, scales = build_grm(reader, args.chunk_size)
    coeffs = compute_coefficients(reader, means, scales, grm, args.h2, args.chunk_size)
    write_results(reader, coeffs, args.output_tsv, args.meta_json, args.h2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
