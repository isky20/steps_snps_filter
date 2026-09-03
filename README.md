
# TArray

TArray is a small collection of shell pipelines for joint quality control (QC) and harmonisation of **whole-genome sequencing (WGS)** and **SNP-array** data, developed in the context of building a Turkish population–specific SNP array and downstream PRS / GWAS analyses.

The scripts implement standard GWAS/WGS QC steps (sample- and variant-level) based on published recommendations and QC table (missingness, MAF, HWE, Ti/Tv, depth, singletons, ambiguous SNPs, etc.).
<img width="1083" height="433" alt="sd drawio (2)" src="https://github.com/user-attachments/assets/fcbd39fe-aed0-4b35-b00b-11d435a18a42" />

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
---
# steps:

1- Fix target array size (1.2M) and genome build (GRCh38).

2- Collect inputs: Paragenomic backbone SNP list + Turkish WGS VCF + “common disorders” PRS scoring files (PGS Catalog).

3- Standardise everything to the same build and variant ID format (prefer chr:pos:ref:alt) and normalise alleles (uppercase, left-normalise, biallelic SNPs).

4- QC Turkish WGS to create a clean “truth set” (PASS, biallelic SNPs) and compute Turkish AF/MAF.

5- Annotate the backbone SNPs with Turkish AF/MAF and presence/quality in WGS.

6- Define your v1 “common disorders” package (e.g., 30–50 diseases) and for each disease select 1–2 best PRS models from PGS Catalog.

7- Standardise all PRS files (build, IDs, effect allele, weight) into one consistent scoring format.

8- Calculate PRS–backbone overlap (which PRS SNPs are already on the backbone vs missing).

7- Set an array “slot budget” by bucket (60–80% backbone tags, 10–25% PRS add-ons, 5–20% Turkish-enriched/clinical).

9- Keep the Paragenomic SNP set as the core backbone, but remove obvious redundancies in Turkish LD (free slots without losing coverage).

10- Add missing PRS variants that are feasible to genotype; if a PRS SNP is problematic, choose a Turkish LD proxy (r² ≥ 0.8).

11- Add Turkish-enriched variants from WGS (common/moderate frequency in Turkish, underrepresented globally, relevant to target diseases).

12- Run manufacturer probe-design screening; drop failed SNPs and replace them with the best Turkish proxies.

13- Freeze “candidate array v1” with labels per SNP (Backbone / PRS / Turkish-enriched / Clinical).

14- Validate
- Validate by WGS simulation: split WGS into reference + test, mask test to only candidate array SNPs, impute, compare imputed vs true WGS (dosage r²/INFO by MAF bins).

- Validate PRS performance: compute PRS from true WGS and from imputed simulated-array, then check correlation/error per disorder.

- Iterate: adjust tags/add-ons until imputation and PRS metrics meet your targets.

15- Finalise “array v1” deliverables: final SNP manifest + replacement log + validation report + file/ID/build documentation.
