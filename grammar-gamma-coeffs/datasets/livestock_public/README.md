# Livestock Genotype Datasets

This folder currently contains three genotype datasets for validating GRAMMAR-Gamma under
full-GRM and LOCO mixed models.

## Summary

| Dataset | Files | Species / breed | Platform | Samples | SNPs | Reference | Planned use |
| --- | --- | --- | --- | ---: | ---: | --- | --- |
| `Cattle_geno` | `Cattle_geno.bed/.bim/.fam` | Cattle, German Holstein bulls | Illumina Bovine SNP50 BeadChip | 5,024 | 42,551 | Zhang et al. 2015, DOI `10.1534/g3.114.016261` | Full-GRM only |
| `genotypedata.7z` | `genotype.ped`, `genotype.map` in archive | Cattle, Chinese Holstein cows | Illumina BovineSNP50 BeadChip | 2,510 | ~44,045 | Huang et al. 2019, DOI `10.1093/jas/skz031`; Figshare DOI `10.6084/m9.figshare.5353498.v1` | Full-GRM and LOCO |
| `figshare_sheep_karacabey_merino_50K.vcf` | VCF | Sheep, Karacabey Merino | Illumina OVINE 50K BeadChip | 734 | 35,886 | Figshare DOI `10.6084/m9.figshare.29184098.v1`; related paper DOI `10.1186/s12917-025-04957-9` | Full-GRM and LOCO |

## Notes for analysis

- `Cattle_geno` has usable genotypes but a fake `.bim` map: all SNPs are assigned to
  chromosome `1` with synthetic evenly spaced positions. It should therefore be used only
  for full-GRM analyses that do not require credible chromosome boundaries.
- `genotypedata.7z` is the Chinese Holstein SNP50 dataset added to this folder. The archive
  contains `genotype.ped` and `genotype.map`, which are suitable for both full-GRM and
  LOCO analyses after standard PLINK conversion.
- `figshare_sheep_karacabey_merino_50K.vcf` is a single-breed sheep dataset with a usable
  marker map and can be analyzed under both full-GRM and LOCO.

## Source links

- Zhang et al. 2015: `https://pmc.ncbi.nlm.nih.gov/articles/PMC4390577/`
- Huang et al. 2019: `https://pmc.ncbi.nlm.nih.gov/articles/PMC6396242/`
- Chinese Holstein Figshare record: `https://figshare.com/articles/dataset/genotype_data_of_2510_dairy_cows/5353498`
- Karacabey Merino Figshare record: `https://figshare.com/articles/dataset/Karacabey_merino_50K_genotype_data/29184098`
