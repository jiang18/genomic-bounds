# gwas-bounds

Code and archived data for the manuscript
**"Genomic Dimensionality Bounds Mixed-Model Association Power and
Fine-Mapping Resolution"** (Jiang, in submission).

## Overview

This repository accompanies a theoretical study showing that
mixed-model GWAS power saturation and fine-mapping resolution limits in
livestock are governed by the effective genomic dimensionality $M_e = 4 N_e L$
(where $N_e$ is the effective population size and $L$ the genome length
in Morgans, so $M_e$ is low when $N_e$ is small, as in livestock). The
framework gives these limits in closed form, including a per-SNP
detection floor on the proportion of phenotypic variance explained,
$q_{\min} \approx 30 \cdot h^2/M_e$, and connects them to a well-established theoretical formula for genomic prediction accuracy, explaining why
prediction is comparatively easy in the very populations where SNP-level
mapping is hard. It is validated through coalescent simulations
(msprime) and per-SNP GRAMMAR-Gamma coefficients from three publicly
available livestock chip-genotype datasets.

## Repository layout

```
gwas-bounds/
├── README.md                       (this file)
├── LICENSE                         MIT license
├── requirements.txt                Python dependencies
│
├── msprime/                        Coalescent simulation + phenotype simulation
│   ├── simulate_genotypes.py       msprime coalescent simulation (two-epoch demography)
│   ├── compute_eigenvalues_by_N.py LD-matrix eigenvalue spectra at multiple N
│   ├── npy_to_plink.py             Convert .npy genotype dumps to PLINK BED
│   ├── select_focal_snps.sh        Pre-select focal SNPs (5 per Ne) by MAF target
│   ├── project_focal_snps.py       Project focal SNPs onto LD eigenvectors
│   ├── run_all_scenarios.sh        Run all (h2, beta2, N) phenotype-sim scenarios
│   ├── run_reml.sh                 REML + score-chi-square per replicate (uses SLEMM)
│   ├── sim_step1.sh                Helper for stage 1 of the simulation pipeline
│   ├── compute_cbar_loco.py        Average GRAMMAR-Gamma coefficients (full-GRM and LOCO)
│   ├── collect_results.py          Aggregate per-replicate outputs into summary tables
│   │
│   ├── cbar_Ne20.csv, cbar_Ne50.csv, cbar_Ne100.csv
│   │                               Average GRAMMAR-Gamma tables (full-GRM and LOCO) per Ne
│   ├── eigen_out/                  Archived LD-matrix eigenvalue arrays per (Ne, N)
│   │   └── Ne{20,50,100}/eigenvalues_N{N}.npy
│   └── plink_out/                  Archived phenotype-simulation outputs per scenario
│       └── Ne{20,50,100}/
│           ├── focal_snps_info.csv
│           └── h2_{0.1,0.3,0.5}__beta2_{0.001,0.01}/N_{5000..100000}/
│               ├── summary.csv           per-replicate NCP and PVE summary
│               ├── per_replicate_chi2.csv
│               └── blue_bias_check.csv
│
├── grammar-gamma-coeffs/           Real-data per-SNP GRAMMAR-Gamma analysis
│   ├── GRAMMAR_Gamma_Analysis_Summary.md
│   ├── compute_grammar_gamma_coeffs.py     Full-GRM per-SNP c_l
│   ├── compute_loco_grammar_gamma_coeffs.py LOCO per-SNP c_LOCO
│   ├── analysis_outputs/grammar_gamma/results/
│   │   ├── chinese_holstein_coeffs.tsv      chinese_holstein_loco_coeffs.tsv
│   │   ├── karacabey_merino_coeffs.tsv      karacabey_merino_loco_coeffs.tsv
│   │   ├── german_holstein_coeffs.tsv
│   │   ├── *_meta.json                      Run-metadata for each analysis
│   │   ├── summary.tsv                      Per-dataset c_l summary
│   │   └── loco_vs_full_summary.tsv         Side-by-side full-GRM vs LOCO summary
│   └── datasets/livestock_public/
│       └── README.md               Pointers to public sources (Figshare DOIs);
│                                   the raw genotype files themselves are NOT
│                                   archived here. See this README for download
│                                   links to the Chinese Holstein, German Holstein,
│                                   and Karacabey Merino chip panels.
│
└── figures-tables/                 Scripts producing all published figures and tables
    ├── README.md                   Detailed file-by-file documentation
    ├── style.py                    Shared rcParams, palettes, and helper functions
    ├── compute_table1_eigenspectrum.py    Table 1 + Tables S6 and S7
    ├── compute_table2_thresholds.py       Table 2
    ├── compute_tableS1_cl_summary.py      Table S1
    ├── compute_tableS2_per_snp_ncp.py     Table S2
    ├── compute_tableS3_sigmoid_sum.py     Table S3
    ├── compute_tableS4_cbar.py            Table S4
    ├── compute_tableS5_group_b.py         Table S5
    ├── fig1_cl_distribution.py            Figure 1
    ├── fig2_ncp_scatter.py                Figure 2
    ├── fig3_sigmoid_sum.py                Figure 3
    ├── fig4_cbar_divergence.py            Figure 4
    ├── fig_supp_sb_sa_ratio.py            Figure S1
    ├── fig_supp_group_b.py                Figure S2
    ├── tables/                            Published CSV tables (regenerable)
    └── figures/                           Published PDF figures (regenerable)
```

## What is archived vs. what is regenerable

| Archived in this repo | Not archived (regenerable) |
|---|---|
| • LD-matrix eigenvalue arrays (`msprime/eigen_out/`)<br>• Per-replicate score $\chi^2$ statistics (`msprime/plink_out/.../per_replicate_chi2.csv`)<br>• Per-population focal-SNP information (`msprime/plink_out/.../focal_snps_info.csv`)<br>• Average GRAMMAR-Gamma tables for full-GRM and LOCO (`msprime/cbar_Ne*.csv`)<br>• Per-SNP GRAMMAR-Gamma coefficient tables for the three livestock chip panels (`grammar-gamma-coeffs/analysis_outputs/grammar_gamma/results/`)<br>• Published manuscript figures (PDF) and tables (CSV) | • Raw simulated genotype matrices (multi-GB; can be regenerated from documented msprime random seeds)<br>• Raw chip-genotype data for the three real-data panels (download from the Figshare DOIs listed in `grammar-gamma-coeffs/datasets/livestock_public/README.md`)<br>• QC intermediates produced by the GRAMMAR-Gamma analysis pipeline |

## Dependencies

- Python ≥ 3.10
- Python packages listed in `requirements.txt` (`numpy`, `pandas`, `matplotlib`, `scipy`, `msprime`)
- PLINK v1.9 or v2 (for genotype-format conversion in the simulation pipeline)
- SLEMM (for REML + score-chi-square in phenotype simulations);
  available at <https://github.com/jiang18/slemm>

Install Python dependencies with:
```bash
pip install -r requirements.txt
```

## Reproducing the published figures and tables

If you only want to regenerate the figures and tables from the
archived intermediate data (without re-running the simulations), run:

```bash
cd figures-tables

# Tables (run any subset)
python compute_table1_eigenspectrum.py --eigen_dir ../msprime/eigen_out
python compute_table2_thresholds.py
python compute_tableS1_cl_summary.py
python compute_tableS2_per_snp_ncp.py --plink_dir ../msprime/plink_out
python compute_tableS3_sigmoid_sum.py --eigen_dir ../msprime/eigen_out
python compute_tableS4_cbar.py
python compute_tableS5_group_b.py --eigen_dir ../msprime/eigen_out

# Figures (run any subset)
python fig1_cl_distribution.py        --out figures/fig1_cl_distribution.pdf
python fig2_ncp_scatter.py            --plink_dir ../msprime/plink_out --out figures/fig2_ncp_scatter.pdf
python fig3_sigmoid_sum.py            --eigen_dir ../msprime/eigen_out --out figures/fig3_sigmoid_sum.pdf
python fig4_cbar_divergence.py                                          --out figures/fig4_cbar_divergence.pdf
python fig_supp_sb_sa_ratio.py        --eigen_dir ../msprime/eigen_out --out figures/fig_supp_sb_sa_ratio.pdf
python fig_supp_group_b.py            --eigen_dir ../msprime/eigen_out --out figures/fig_supp_group_b.pdf
```

See `figures-tables/README.md` for per-script details.

## Reproducing the simulations and the real-data analysis from scratch

### msprime simulations (per Ne)

```bash
cd msprime

# Ne = 50 example; replace 50 with 20 or 100 for other populations
python simulate_genotypes.py 50 100000 50000 sim_out/Ne50
python compute_eigenvalues_by_N.py sim_out/Ne50 eigen_out/Ne50 \
       5000 10000 20000 50000 100000
python npy_to_plink.py sim_out/Ne50 plink_out/Ne50
bash select_focal_snps.sh plink_out/Ne50
python project_focal_snps.py sim_out/Ne50 eigen_out/Ne50 plink_out/Ne50 \
       5000 10000 20000 50000
bash run_all_scenarios.sh Ne50 100 28        # 100 replicates, 28 threads
```

The intermediate `sim_out/` directory holds the raw simulated
genotype matrices and is not archived in this repo (multi-GB per
population); it is regenerated by `simulate_genotypes.py` from the
documented random seed.

### Real-data per-SNP GRAMMAR-Gamma coefficients

```bash
cd grammar-gamma-coeffs

# Download genotype data per links in
# datasets/livestock_public/README.md, place under
# datasets/livestock_public/, then:
python compute_grammar_gamma_coeffs.py        # Full-GRM c_l for all three panels
python compute_loco_grammar_gamma_coeffs.py   # LOCO c_LOCO for the two panels with maps
```

See `grammar-gamma-coeffs/GRAMMAR_Gamma_Analysis_Summary.md` for
analysis details.

## Citing this repository

If you use this code or the archived data in your work, please cite
the accompanying paper (citation will be added on publication).

## License

MIT — see `LICENSE`. Archived intermediate data are likewise free to
reuse with attribution. Raw genotype data for the three livestock
panels remain under the licenses of their respective original
deposits (see `grammar-gamma-coeffs/datasets/livestock_public/README.md`).
