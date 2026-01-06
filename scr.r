#SNPs Data Selection

"
SNPs Data Selection
Standard GWAS QC:
  - Genotyping rate > 0.90
  - Remove SNPs with > 10% missing data
  - Heterozygosity: remove samples > mean(proportion of heterozygous sites)
    ± 3 * sd(proportion of heterozygous sites)
  - Exclude individuals with sex discordance
  - INFO score > 0.3–0.8
  - MAF > 0.005–0.001
  - LD score
  - HWE P > 1e-6
"

#Target Data
"
* More than 100 individuals
* Genome build
* No sample overlap with SNPs GWAS samples

Standard GWAS For samples:
  -  Individuals genotype rate >0.99
  -  Remove sample with >10% missing data
  -  Heterozygosity. Remove samples > mean(proportion of heterozygous sites ) ± 3 * sd(proportion of heterozygous sites)
  -  Exclude individuals with sex discordance
  -  Remove duplicates and relatedness samples. PI_HAT > 0.8 For SNPs:
  -  Imputation information score (INFO score). > 0.3 – 0.8
  -  Minor Allele Frequency (MAF). > 0.01 - 0.005
  -  Linkage-disequilibrium (LD) score
  -  Hardy-Weinberg Equilibrium (HWE) p > 1x10-6
"
