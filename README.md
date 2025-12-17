
# TArray

TArray is a small collection of shell pipelines for joint quality control (QC) and harmonisation of **whole-genome sequencing (WGS)** and **SNP-array** data, developed in the context of building a Turkish population–specific SNP array and downstream PRS / GWAS analyses.

The scripts implement standard GWAS/WGS QC steps (sample- and variant-level) based on published recommendations and QC table (missingness, MAF, HWE, Ti/Tv, depth, singletons, ambiguous SNPs, etc.).
<img width="1083" height="433" alt="sd drawio (1)" src="https://github.com/user-attachments/assets/ade1a4f9-af36-4b02-9b42-e9ce818be114" />

---

## Contents

- `wgs_qc.sh`  
  End-to-end QC for a multi-sample WGS VCF:
  - sample-level QC using `bcftools stats` (Ti/Tv, depth, singletons, total variants, missingness),
  - variant-level filtering on QUAL, depth/genotype quality, call rate,
  - conversion to PLINK2 format for downstream GWAS/PRS.

- `snps_qc.sh`  
  QC pipeline for SNP-array data in PLINK2 format:
  - sample and SNP missingness filters,
  - heterozygosity / sex checks,
  - Hardy–Weinberg filtering,
  - MAF filters for common variants,
  - removal of strand-ambiguous SNPs to facilitate merging with WGS.

- `count_matched_snps.sh`  
  Helper script to compare SNP sets between WGS and SNP array:
  - counts overlapping SNPs,
  - lists SNPs present only in WGS or only in the array,
  - useful for array design, imputation planning, and platform concordance checks.

---

## Requirements

- `bcftools` (for VCF filtering and stats)  
- `plink2` (for GWAS/PRS-ready formats and SNP-level QC)  
- `bash` and standard Unix tools (`grep`, `awk`, etc.)  

---

## Quick start

1. **Clone the repository**

   ```bash
   git clone https://github.com/isky20/TArray.git
   cd TArray
