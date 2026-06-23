# Scripts for Tables and Figures

## File overview

| File | Purpose |
|---|---|
| `style.py` | Shared rcParams, colours, constants, helper functions. Imported by all other scripts. |
| `compute_table1_eigenspectrum.py` | Table 1 (long form) and Tables S6/S7: eigenvalue spectral parameters at three Me definitions |
| `compute_table2_thresholds.py` | Table 2 (transposed, species as columns): detection and fine-mapping thresholds |
| `compute_tableS1_cl_summary.py` | Table S1: per-SNP c_l summary statistics (full-GRM and LOCO) for Fig 1 |
| `compute_tableS2_per_snp_ncp.py` | Table S2: 360-row per-SNP NCP for Fig 2, including MAF |
| `compute_tableS3_sigmoid_sum.py` | Table S3: S(N), S_A(N), S_B(N) values for Fig 3 |
| `compute_tableS4_cbar.py` | Table S4: c-bar full-GRM and LOCO + bound for Fig 4 |
| `compute_tableS5_group_b.py` | Table S5: Group B properties (S_A, S_B, S_B/S_A, U_B, Jensen gap) for Figs S1 and S2 |
| `fig1_cl_distribution.py` | Figure 1: 2 x 3 grid of per-SNP c_l KDEs (full-GRM and LOCO) |
| `fig2_ncp_scatter.py` | Figure 2: theoretical vs empirical NCP scatter (360 points) |
| `fig3_sigmoid_sum.py` | Figure 3: S(N) saturation at Me (3 panels per Ne) |
| `fig4_cbar_divergence.py` | Figure 4: full-GRM vs LOCO c-bar divergence (3 panels per Ne) |
| `fig5_gblup_reliability.py` | Figure 5: empirical vs predicted in-sample GBLUP reliability (36 scenarios) |
| `fig_supp_sb_sa_ratio.py` | Supplementary Figure S1: S_B/S_A ratio (3 panels per Ne) |
| `fig_supp_group_b.py` | Supplementary Figure S2: S_B vs U_B (3 panels per Ne) |

## Dependencies

```bash
pip install numpy pandas matplotlib scipy
```

No SciencePlots dependency --- style is implemented directly in `style.py`.
Fonts are embedded as TrueType (Type 42) in all PDF/PS output, as required
for journal submission. Arial is used on Windows; Liberation Sans or DejaVu
Sans are used as fallbacks on Linux.

## Directory structure expected

`figures-tables/` is a sibling of `msprime/` and `grammar-gamma-coeffs/`
in the repository root. The scripts read from those sibling pipelines
via relative paths:

```
<repo root>/
├── msprime/
│   ├── cbar_Ne20.csv  cbar_Ne50.csv  cbar_Ne100.csv
│   ├── eigen_out/
│   │   └── Ne20/  Ne50/  Ne100/
│   │           eigenvalues_N{N}.npy
│   └── plink_out/
│       └── Ne{ne}/
│           ├── focal_snps_info.csv
│           └── h2_{h2}__beta2_{beta2}/N_{N}/summary.csv
├── grammar-gamma-coeffs/
│   └── analysis_outputs/grammar_gamma/results/
│       ├── chinese_holstein_coeffs.tsv      chinese_holstein_loco_coeffs.tsv
│       ├── karacabey_merino_coeffs.tsv      karacabey_merino_loco_coeffs.tsv
│       └── german_holstein_coeffs.tsv
└── figures-tables/           (this folder)
    ├── style.py, compute_*.py, fig*.py
    ├── tables/               (CSV outputs)
    └── figures/              (PDF outputs)
```

## Figure and table numbering in the manuscript

| Display item | Section | Script | Output |
|---|---|---|---|
| Figure 1: per-SNP c_l KDEs (full + LOCO) | 4.1 | `fig1_cl_distribution.py` | `figures/fig1_cl_distribution.pdf` |
| Figure 2: theoretical vs empirical NCP | 4.2 | `fig2_ncp_scatter.py` | `figures/fig2_ncp_scatter.pdf` |
| Figure 3: S(N) saturation | 4.3 | `fig3_sigmoid_sum.py` | `figures/fig3_sigmoid_sum.pdf` |
| Figure 4: c-bar divergence | 4.3 | `fig4_cbar_divergence.py` | `figures/fig4_cbar_divergence.pdf` |
| Figure 5: GBLUP reliability | 4.4 | `fig5_gblup_reliability.py` | `figures/fig5_gblup_reliability.pdf` |
| Table 1: eigenspectrum parameters | 4.3 | `compute_table1_eigenspectrum.py` | `tables/table1_eigenspectrum.csv` |
| Table 2: detection and fine-mapping | 4.5 | `compute_table2_thresholds.py` | `tables/table2_thresholds.csv` |
| Supp. Figure S1: S_B/S_A ratio | 4.3 | `fig_supp_sb_sa_ratio.py` | `figures/fig_supp_sb_sa_ratio.pdf` |
| Supp. Figure S2: S_B vs U_B | 4.3 | `fig_supp_group_b.py` | `figures/fig_supp_group_b.pdf` |
| Supp. Table S1: c_l summary | 4.1 | `compute_tableS1_cl_summary.py` | `tables/tableS1_cl_summary.csv` |
| Supp. Table S2: per-SNP NCP | 4.2 | `compute_tableS2_per_snp_ncp.py` | `tables/tableS2_per_snp_ncp.csv` |
| Supp. Table S3: S(N) values | 4.3 | `compute_tableS3_sigmoid_sum.py` | `tables/tableS3_sigmoid_sum.csv` |
| Supp. Table S4: c-bar values | 4.3 | `compute_tableS4_cbar.py` | `tables/tableS4_cbar.csv` |
| Supp. Table S5: Group B properties (S_A, S_B, U_B, Jensen gap) | 4.3 | `compute_tableS5_group_b.py` | `tables/tableS5_group_b.csv` |
| Supp. Table S6: EIG98 spectrum | 4.3 | `compute_table1_eigenspectrum.py` | `tables/tableS6_eig98.csv` |
| Supp. Table S7: EIG99 spectrum | 4.3 | `compute_table1_eigenspectrum.py` | `tables/tableS7_eig99.csv` |
| Supp. Table S8: GBLUP reliability | 4.4 | ≈`../msprime/gblup_reliability.csv` | `tables/tableS8_gblup_reliability.csv` |

## Usage

### Tables

All commands are run from inside `figures-tables/`. The `--eigen_dir`,
`--plink_dir`, `--cbar_dir`, and `--coeff_dir` arguments point upward
into the sibling `msprime/` and `grammar-gamma-coeffs/` folders.

```bash
# Table 1 (long form) + Tables S6, S7
python compute_table1_eigenspectrum.py --eigen_dir ../msprime/eigen_out

# Table 2 (transposed)
python compute_table2_thresholds.py

# Supplementary tables
python compute_tableS1_cl_summary.py
python compute_tableS2_per_snp_ncp.py --plink_dir ../msprime/plink_out
python compute_tableS3_sigmoid_sum.py --eigen_dir ../msprime/eigen_out
python compute_tableS4_cbar.py
python compute_tableS5_group_b.py --eigen_dir ../msprime/eigen_out
```

### Figures

```bash
# Main figures
python fig1_cl_distribution.py --out figures/fig1_cl_distribution.pdf
python fig2_ncp_scatter.py --plink_dir ../msprime/plink_out --out figures/fig2_ncp_scatter.pdf
python fig3_sigmoid_sum.py --eigen_dir ../msprime/eigen_out --out figures/fig3_sigmoid_sum.pdf
python fig4_cbar_divergence.py --out figures/fig4_cbar_divergence.pdf
python fig5_gblup_reliability.py --csv ../msprime/gblup_reliability.csv --out figures/fig5_gblup_reliability.pdf

# Supplementary figures
python fig_supp_sb_sa_ratio.py --eigen_dir ../msprime/eigen_out --out figures/fig_supp_sb_sa_ratio.pdf
python fig_supp_group_b.py --eigen_dir ../msprime/eigen_out --out figures/fig_supp_group_b.pdf
```

### Raw data

`tables/raw/` holds the per-replicate simulation outputs underlying the
NCP analysis (`empirical_ncp.csv`, `ncp_validation_summary.csv`).  These
are kept for reproducibility but are not publication tables.

## Species parameters (from Supplementary Methods S15)

| Species | Ne | L (Morgans) | Me = 4NeL |
|---|---|---|---|
| Cattle | 100 | 25 | 10,000 |
| Pig | 50 | 20 | 4,000 |
| Chicken | 50 | 30 | 6,000 |
| Human | 10,000 | 35 | 1,400,000 |

These are the authoritative values used in Table 2 and throughout the
manuscript. The msprime simulations use L = 25 Morgans (cattle genome)
for all Ne values.

## Key design decisions

### h2 coverage
All figures and tables cover h2 = 0.1, 0.3, 0.5 simultaneously. The
real-data GRAMMAR-Gamma analysis also includes h2 = 0.7 and 0.9, but
these additional values are discussed in the text rather than shown in
every figure or table.

### N used in eigenvalue tables
Table 1 (and S7/S8) compute eigenvalue spectral parameters at the
*largest* available N for each Ne (100,000 for Ne=20 and 50; 200,000
for Ne=100) so that the eigenvalue spectrum is well-resolved.

### Genome length L
Default L = 25 Morgans (bovine genome). To change, edit `L_MORGANS` in
`style.py`.

### Figure output format
PDF is preferred for submission (vector graphics, fonts embedded as
TrueType). For quick inspection use `.png`. The `--out` argument accepts
any matplotlib-supported extension.
