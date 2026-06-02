# GRAMMAR-Gamma Coefficient Analysis Summary

## Overview

This note summarizes the computation of GRAMMAR-Gamma coefficients across three livestock genotype datasets and documents the QC procedures used before coefficient calculation.

The full-GRM per-SNP coefficient was computed as

$$
 c_l(h^2) = \frac{x_l^T \left(Gh^2 + I(1-h^2)\right)^{-1} x_l}{N},
$$

for the working-model heritability values

- `h2 = 0.1`
- `h2 = 0.3`
- `h2 = 0.5`
- `h2 = 0.7`
- `h2 = 0.9`

The three datasets were:

- `Cattle_geno`: German Holstein bulls, used for full-GRM only
- `genotypedata.7z`: Chinese Holstein cows, used for full-GRM and LOCO
- `figshare_sheep_karacabey_merino_50K.vcf`: Karacabey Merino sheep, used for full-GRM and LOCO

Relevant local files:

- Full-GRM script: `compute_grammar_gamma_coeffs.py`
- LOCO script: `compute_loco_grammar_gamma_coeffs.py`
- Dataset summary: `datasets/livestock_public/README.md`
- Full-GRM summary table: `analysis_outputs/grammar_gamma/results/summary.tsv`
- LOCO vs full summary table: `analysis_outputs/grammar_gamma/results/loco_vs_full_summary.tsv`

## Datasets

| Dataset | Species / breed | Source format | Planned use |
| --- | --- | --- | --- |
| `Cattle_geno` | German Holstein bulls | PLINK BED/BIM/FAM | Full-GRM only |
| `genotypedata.7z` | Chinese Holstein cows | PED/MAP in archive | Full-GRM and LOCO |
| `figshare_sheep_karacabey_merino_50K.vcf` | Karacabey Merino sheep | VCF | Full-GRM and LOCO |

Notes:

- `Cattle_geno` has a synthetic `.bim` map and was therefore used only for full-GRM analyses.
- Chinese Holstein and Karacabey Merino have usable chromosome assignments and were analyzed under both full-GRM and LOCO.

## QC Procedure

### Common variant filters

The same core QC thresholds were applied across datasets:

- `MAF >= 0.01`
- `HWE p >= 1e-6`

### Filters not applied

The following were intentionally not applied:

- no SNP missingness filter such as `--geno`
- no sample missingness filter such as `--mind`
- no LD pruning
- no relatedness filtering
- no control-only HWE for the Chinese Holstein case-control dataset

### Missing genotype handling in coefficient computation

Missing genotypes were handled inside the coefficient scripts as follows:

1. PLINK missing calls were decoded as `NaN`.
2. SNP means were computed from observed genotypes only.
3. Missing centered genotypes were set to `0`, i.e. mean-imputed after centering.
4. SNPs were standardized by `sqrt(2p(1-p))`.

Thus, coefficient computation used mean-imputed standardized genotypes rather than dropping SNPs or samples for low missingness.

## Dataset-Specific QC

### 1. German Holstein (`Cattle_geno`)

QC log:

- `analysis_outputs/grammar_gamma/qc/german_holstein_qc.log`

PLINK command options recorded in the log:

- `--bfile datasets/livestock_public/Cattle_geno`
- `--maf 0.01`
- `--hwe 1e-6`
- `--make-bed`

Initial data:

- `5024` animals
- `42551` SNPs

QC results:

- total genotyping rate: `1.000000`
- removed by HWE: `157`
- removed by MAF: `178`
- final retained: `5024` animals, `42216` SNPs

Special note:

- no autosome-only filter was applied because the map is synthetic and this dataset was used only for full-GRM analysis

### 2. Chinese Holstein (`genotypedata.7z`)

Initial QC log:

- `analysis_outputs/grammar_gamma/qc/chinese_holstein_qc.log`

PLINK command options recorded in the log:

- `--file /tmp/chinese_holstein_raw/genotype`
- `--cow`
- `--maf 0.01`
- `--hwe 1e-6`
- `--make-bed`

Initial data:

- `2510` cattle
- `44046` SNPs

QC results before autosome restriction:

- total genotyping rate: `0.997582`
- removed by HWE: `16`
- removed by MAF: `16`
- retained after basic QC: `2510` cattle, `44014` SNPs

Case-control note:

- the dataset contains `364` cases and `2146` controls
- PLINK emitted an HWE warning because case/control counts were imbalanced
- global HWE filtering was retained intentionally; controls-only HWE was not used

Autosome restriction log:

- `analysis_outputs/grammar_gamma/qc/chinese_holstein_auto_qc.log`

Autosome filter:

- `--autosome --cow`

Autosome-only final set:

- `2510` cattle
- `42775` SNPs

These autosome-only data were used for both full-GRM and LOCO analyses.

### 3. Karacabey Merino (`figshare_sheep_karacabey_merino_50K.vcf`)

Initial QC log:

- `analysis_outputs/grammar_gamma/qc/karacabey_merino_qc.log`

PLINK command options recorded in the log:

- `--vcf datasets/livestock_public/figshare_sheep_karacabey_merino_50K.vcf`
- `--double-id`
- `--allow-extra-chr`
- `--maf 0.01`
- `--hwe 1e-6`
- `--make-bed`

Initial data:

- `734` individuals
- `35886` SNPs

QC results before autosome restriction:

- total genotyping rate: `0.995616`
- removed by HWE: `0`
- removed by MAF: `0`
- retained after basic QC: `734` individuals, `35886` SNPs

Additional notes:

- PLINK wrote a `.nosex` file because sex was ambiguous in the VCF-derived FAM file
- PLINK warned about nonmissing nonmale Y chromosome genotypes

Autosome restriction log:

- `analysis_outputs/grammar_gamma/qc/karacabey_merino_auto_qc.log`

Autosome filter:

- `--chr 1-26`

Autosome-only final set:

- `734` individuals
- `35774` SNPs

These autosome-only data were used for both full-GRM and LOCO analyses.

## Final Datasets Used for Coefficient Computation

### Full-GRM

| Dataset | Samples | SNPs | Autosomes only |
| --- | ---: | ---: | --- |
| German Holstein | 5024 | 42216 | No |
| Chinese Holstein | 2510 | 42775 | Yes |
| Karacabey Merino | 734 | 35774 | Yes |

### LOCO

| Dataset | Samples | SNPs | Chromosomes |
| --- | ---: | ---: | ---: |
| Chinese Holstein | 2510 | 42775 | 29 |
| Karacabey Merino | 734 | 35774 | 26 |

## Full-GRM Coefficient Computation

The full-GRM coefficient analysis was implemented in `compute_grammar_gamma_coeffs.py`.

Procedure:

1. Read PLINK BED data directly.
2. Decode genotypes into dosage values `0/1/2`, with missing calls as `NaN`.
3. Standardize each SNP using observed allele frequency.
4. Construct the exact GRM:

   $$
   G = \frac{ZZ^T}{M}.
   $$

5. For each working-model heritability value `h2`, form:

   $$
   V = Gh^2 + I(1-h^2).
   $$

6. Compute the per-SNP coefficient:

   $$
   c_l = \frac{x_l^T V^{-1} x_l}{N}.
   $$

The script uses a Cholesky factorization of `V` and processes SNPs in chunks.

Full-GRM output files:

- `analysis_outputs/grammar_gamma/results/german_holstein_coeffs.tsv`
- `analysis_outputs/grammar_gamma/results/chinese_holstein_coeffs.tsv`
- `analysis_outputs/grammar_gamma/results/karacabey_merino_coeffs.tsv`

Metadata files:

- `analysis_outputs/grammar_gamma/results/german_holstein_meta.json`
- `analysis_outputs/grammar_gamma/results/chinese_holstein_meta.json`
- `analysis_outputs/grammar_gamma/results/karacabey_merino_meta.json`

## LOCO Coefficient Computation

LOCO coefficients were computed only for the two datasets with usable chromosome maps, using `compute_loco_grammar_gamma_coeffs.py`.

Procedure:

1. Start from the autosome-filtered PLINK data.
2. Build the total standardized cross-product matrix.
3. For each chromosome:
   - identify all SNPs on that chromosome
   - remove those SNPs from the GRM construction
   - construct the leave-one-chromosome-out GRM from the remaining SNPs
4. For each working-model heritability value `h2`, compute chromosome-specific LOCO coefficients for SNPs on the held-out chromosome.

LOCO output files:

- `analysis_outputs/grammar_gamma/results/chinese_holstein_loco_coeffs.tsv`
- `analysis_outputs/grammar_gamma/results/karacabey_merino_loco_coeffs.tsv`

LOCO metadata files:

- `analysis_outputs/grammar_gamma/results/chinese_holstein_loco_meta.json`
- `analysis_outputs/grammar_gamma/results/karacabey_merino_loco_meta.json`

Comparison table:

- `analysis_outputs/grammar_gamma/results/loco_vs_full_summary.tsv`

## Full-GRM Summary Statistics

The consolidated full-GRM summary is in:

- `analysis_outputs/grammar_gamma/results/summary.tsv`

This file reports, for each dataset and each `h2` value:

- post-QC sample size
- post-QC SNP count
- whether autosomes-only filtering was used
- mean coefficient
- standard deviation
- coefficient of variation
- minimum coefficient
- maximum coefficient
- a short note about dataset-specific handling

### Mean coefficients

German Holstein:

- `h2=0.1`: `0.7388`
- `h2=0.3`: `0.6297`
- `h2=0.5`: `0.6089`
- `h2=0.7`: `0.6402`
- `h2=0.9`: `0.7680`

Chinese Holstein:

- `h2=0.1`: `0.7729`
- `h2=0.3`: `0.7069`
- `h2=0.5`: `0.7064`
- `h2=0.7`: `0.7488`
- `h2=0.9`: `0.8635`

Karacabey Merino:

- `h2=0.1`: `0.8634`
- `h2=0.3`: `0.8089`
- `h2=0.5`: `0.8122`
- `h2=0.7`: `0.8494`
- `h2=0.9`: `0.9297`

## LOCO vs Full-GRM Summary

The file `analysis_outputs/grammar_gamma/results/loco_vs_full_summary.tsv` summarizes the full-GRM and LOCO results for Chinese Holstein and Karacabey Merino.

### Chinese Holstein mean coefficients

Full-GRM:

- `0.7729`, `0.7069`, `0.7064`, `0.7488`, `0.8635`

LOCO:

- `0.8331`, `0.8554`, `0.9523`, `1.1376`, `1.5428`

### Karacabey Merino mean coefficients

Full-GRM:

- `0.8634`, `0.8089`, `0.8122`, `0.8494`, `0.9297`

LOCO:

- `0.8873`, `0.8713`, `0.9173`, `1.0119`, `1.1854`

## Practical QC Takeaways

The QC strategy was intentionally simple and reproducible.

Strengths:

- consistent MAF and HWE thresholds across datasets
- autosome restriction applied where chromosome maps were usable
- no aggressive pruning or sample removal that would distort the empirical coefficient distribution

Caveats:

- missingness filtering was not applied, though missingness was low overall
- Chinese Holstein used global HWE despite case-control imbalance
- `Cattle_geno` could not be used for chromosome-based LOCO because the map is synthetic

In short, the final coefficient summaries reflect realistic chip data after light, transparent QC rather than heavily curated marker sets.
