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

# load data
CRC_SNP <- read_excel("data/resources/1_CRC_SNPs_list_205_Nature_2023.xlsx")

# The first 5 rows and first 7 columns are shown
CRC_SNP[c(1:5),c(1:7)]

# .txt file generation
SNP_file_CRC <- CRC_SNP[,c(1,2)]
head(SNP_file_CRC)

# two avoid scientific annotation
SNP_file_CRC$`POS (GRCh37)` <- format(SNP_file_CRC$`POS (GRCh37)`, scientific=FALSE)

# create folder data/1_snps_extraction/ and save as .txt
dir.create("data/1_snps_extraction/", recursive=TRUE)
write.table(SNP_file_CRC, "data/1_snps_extraction/1_SNP_file_CRC.txt", quote = F, col.names = F, row.names = F, qmethod = "double", sep = "\t")



