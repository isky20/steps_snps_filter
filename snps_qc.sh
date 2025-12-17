#!/usr/bin/env bash
set -euo pipefail

###############################
# 0. (Optional) Adjust header #
###############################
# Example: rename a header field (e.g. change HWP tag name in the VCF header).
sed 's/H.W.p-Value/HWpValue/g' input.vcf > input_renamed.vcf

bgzip -@ 8 -c allrun1to19.vcf > allrun1to19.vcf.gz
tabix -p vcf allrun1to19.vcf

#########################################
# 1. Add extra INFO fields to VCF header
#########################################
#cat > extra_info_header.hdr <<EOF
##INFO=<ID=chr_id,Number=1,Type=String,Description="Chromosome ID from array annotation">
##INFO=<ID=start,Number=1,Type=String,Description="Start position from array annotation">
##INFO=<ID=stop,Number=1,Type=String,Description="Stop position from array annotation">
##INFO=<ID=strand,Number=1,Type=String,Description="Strand information">
##INFO=<ID=dbsnp_loctype,Number=1,Type=String,Description="dbSNP location type">
##INFO=<ID=in_hapmap,Number=1,Type=String,Description="Present in HapMap">
##INFO=<ID=strand_vs_dbsnp,Number=1,Type=String,Description="Strand compared to dbSNP">
##INFO=<ID=probe_count,Number=1,Type=String,Description="Number of probes">
##INFO=<ID=cytoband,Number=1,Type=String,Description="Cytoband">
##INFO=<ID=chrx_par,Number=1,Type=String,Description="X chromosome pseudoautosomal region flag">
##INFO=<ID=flank,Number=1,Type=String,Description="Flanking sequence">
##INFO=<ID=allele_a,Number=1,Type=String,Description="Allele A from array">
##INFO=<ID=allele_b,Number=1,Type=String,Description="Allele B from array">
##INFO=<ID=ref_allele,Number=1,Type=String,Description="Reference allele from annotation">
##INFO=<ID=alt_allele,Number=1,Type=String,Description="Alternative allele from annotation">
##INFO=<ID=associated_gene,Number=1,Type=String,Description="Associated gene symbol">
##INFO=<ID=genetic_map,Number=1,Type=String,Description="Genetic map position">
##INFO=<ID=microsatellite,Number=1,Type=String,Description="Microsatellite flag/info">
##INFO=<ID=heterozygous_allele_frequencies,Number=.,Type=String,Description="Heterozygous allele frequencies">
##INFO=<ID=allele_frequency_count,Number=.,Type=String,Description="Allele frequency counts">
##INFO=<ID=allele_frequencies,Number=.,Type=String,Description="Allele frequencies">
##INFO=<ID=minor_allele,Number=1,Type=String,Description="Minor allele">
##INFO=<ID=minor_allele_frequency,Number=1,Type=String,Description="Minor allele frequency">
##INFO=<ID=omim,Number=.,Type=String,Description="OMIM annotation">
##INFO=<ID=biomedical,Number=.,Type=String,Description="Biomedical annotation">
##INFO=<ID=annotation_notes,Number=.,Type=String,Description="Annotation notes">
##INFO=<ID=ordered_alleles,Number=.,Type=String,Description="Ordered alleles">
##INFO=<ID=allele_count,Number=.,Type=String,Description="Allele count">
#EOF
###########################################
bcftools annotate \
  -h extra_info_header.hdr \
  -Oz -o vcf_with_info_header.vcf.gz \
  allrun1to19_chr1_22.vcf.gz

tabix vcf_with_info_header.vcf.gz

###########################################
# 2. Strip INFO for PLINK (keep GT only) #
###########################################
bcftools annotate \
  -x INFO \
  -Oz -o vcf_noinfo.vcf.gz  vcf_with_info_header.vcf.gz                          

tabix vcf_noinfo.vcf.gz             
###########################################
###########################################
##########START############################
############################HERE###########
###########################################

sed 's/H.W.p-Value/HWpValue/g' input.vcf > allrun1to19_v1.vcf

bgzip -@ 8 -c allrun1to19_v1.vcf > allrun1to19_v1.vcf.gz
tabix -p vcf allrun1to19_v1.vcf.gz

singularity exec ~/isky20/02.software/sifs/bcftools.sif bcftools view \
-r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
-Ov -o allrun1to19_v2.vcf   allrun1to19_v1.vcf.gz


singularity exec ~/isky20/02.software/sifs/bcftools.sif   bcftools annotate -x INFO -Ov -o allrun1to19_v3.vcf allrun1to19_v2.vcf --force


bgzip -@ 8 -c allrun1to19_v3.vcf > allrun1to19_v3.vcf.gz
tabix -p vcf allrun1to19_v3.vcf.gz

##############################################################
# 3. Convert SNP-array VCF to PLINK2 pgen (raw starting point)
##############################################################
singularity exec plink_combo.sif plink2 \
  --vcf vcf_with_info_header.vcf.gz \
  --snps-only just-acgt \
  --max-alleles 2 \
  --make-pgen \
  --out snps_array_qc0 \
  --allow-extra-chr \
  --threads 32 \
  --split-par hg38

#######################
# 🟦 SAMPLE-LEVEL QC  #
#######################

##############
# 4. Autosome
##############
singularity exec ~/isky20/02.software/sifs/plink2.sif plink2 \
  --pfile snps_array_qc0 \
  --autosome \
  --make-pgen \
  --out snps_array_step1 \
  --allow-extra-chr

####################################################
# 5. Remove samples with high missingness (> 5%)  #
####################################################
singularity exec ~/isky20/02.software/sifs/plink2.sif plink2 \
  --pfile snps_array_step1 \
  --mind 0.05 \
  --make-pgen \
  --out snps_array_step2_mind

######################################
# 6. Heterozygosity outliers (PLINK) #
######################################
singularity exec ~/isky20/02.software/sifs/plink2.sif plink2 \
  --pfile snps_array_step2_mind \
  --make-bed \
  --out snps_array_step2_mind_bed


# Create het.het
singularity exec ~/isky20/02.software/sifs/plink2.sif plink \
  --bfile snps_array_step2_mind_bed \
  --het \
  --out het

###########################################
# 6b. Heterozygosity outliers (R script)  #
###########################################
Rscript - << 'EOF_R'
het <- read.table("het.het", header=TRUE, stringsAsFactors=FALSE)

# mean/SD of F
meanF <- mean(het$F, na.rm=TRUE)
sdF   <- sd(het$F,   na.rm=TRUE)

# outliers
out <- het[abs(het$F - meanF) > 3*sdF, c("FID","IID")]
out <- unique(out)

# write (no header)
if (nrow(out) > 0) {
  write.table(out, "het_outliers.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
} else {
  # safest: write a dummy non-existing ID so --remove doesn't crash on empty file
  writeLines("0 __NONE__", "het_outliers.txt")
}
EOF_R


############################################
# 6c. Remove heterozygosity outliers      #
############################################
singularity exec ~/isky20/02.software/sifs/plink2.sif plink2 \
  --pfile snps_array_step2_mind \
  --remove het_outliers.txt \
  --make-pgen \
  --out snps_array_step3_nohet

############################################
# 7. Remove related individuals (>0.0884) #
############################################
singularity exec ~/isky20/02.software/sifs/plink2.sif plink2 \
  --pfile snps_array_step3_nohet \
  --king-cutoff 0.0884 \
  --make-pgen \
  --out snps_array_step4_unrelated

####################
# 🟧 SNP-LEVEL QC  #
####################

#####################################################################
# 8. Filter SNPs: missingness <10%, MAF>5%, HWE p > 1e−6
#####################################################################
singularity exec ~/isky20/02.software/sifs/plink2.sif plink2 \
  --pfile snps_array_step4_unrelated \
  --geno 0.10 \
  --maf 0.01 \
  --hwe 1e-6 \
  --make-pgen \
  --out snps_array_step6_snpqc

#############################################################
# 9. Remove strand-ambiguous SNPs (A/T, T/A, C/G, G/C)
#############################################################
awk 'NR>1 && (($4=="A" && $5=="T") || \
              ($4=="T" && $5=="A") || \
              ($4=="C" && $5=="G") || \
              ($4=="G" && $5=="C")) {print $3}' \
    snps_array_step6_snpqc.pvar > ambiguous_snps.txt

singularity exec ~/isky20/02.software/sifs/plink2.sif plink2 \
  --pfile snps_array_step6_snpqc \
  --exclude ambiguous_snps.txt \
  --make-pgen \
  --out snps_array_step7_noambig

##################################
# 10. Count remaining SNPs       #
##################################
singularity exec ~/isky20/02.software/sifs/plink2.sif plink2 \
  --pfile snps_array_step7_noambig \
  --write-snplist \
  --out snps_array

wc -l snps_array.snplist
