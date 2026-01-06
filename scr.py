#SNPs Data Selection

"""
SNPs Data Selection
Standard GWAS QC:
  - Genotyping rate > 0.99
  - Remove SNPs with > 5% missing data
  - Heterozygosity: remove samples > mean(proportion of heterozygous sites)
    ± 3 * sd(proportion of heterozygous sites)
  - Exclude individuals with sex discordance
  - INFO score > 0.3–0.8
  - MAF > 0.005–0.001
  - LD score
  - HWE P > 1e-6
"""
