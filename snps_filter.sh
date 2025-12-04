###############################
# 0. (Optional) Adjust header #
###############################
# Example: rename a header field (e.g. change HWP tag name in the VCF header).
# This does a global search/replace in the whole VCF file.
sed 's/OLDWORD/NEWWORD/g' input.vcf > input_renamed.vcf


#########################################
# 1. Add extra INFO fields to VCF header
#########################################
# Create a small header file describing additional INFO annotations
# that you will later attach to the VCF (e.g. from array annotation).
cat > extra_info_header.hdr <<EOF
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
EOF

# Attach the new INFO header lines to the original VCF
bcftools annotate \
  -h extra_info_header.hdr \
  -Oz -o vcf_with_info_header.vcf.gz \
  input_renamed.vcf

# Index the compressed VCF (bgzip + tabix index)
tabix vcf_with_info_header.vcf.gz


###########################################
# 2. Strip INFO for PLINK (keep GT only)
###########################################
# Remove all INFO fields to make a smaller, “GT-only” VCF for PLINK.
bcftools annotate \
  -x INFO \
  -Oz -o vcf_noinfo.vcf.gz \
  vcf_with_info_header.vcf.gz

# Index again after modification
tabix vcf_noinfo.vcf.gz


##############################################################
# 3. Convert SNP-array VCF to PLINK2 pgen (raw starting point)
##############################################################
# -–snps-only just-acgt : keep only SNPs with A/C/G/T alleles
# --max-alleles 2       : bi-allelic variants only
# --allow-extra-chr     : allow non-standard chromosome names
# --split-par hg38      : handle pseudoautosomal regions on X/Y (hg38)
singularity exec plink_combo.sif plink2 \
  --vcf vcf_noinfo.vcf.gz \
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
# Keep autosomal variants only (chr1–22) for downstream sample QC.
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_qc0 \
  --autosome \
  --make-pgen \
  --out snps_array_step1 \
  --allow-extra-chr


####################################################
# 5. Remove samples with high missingness (> 5%)
####################################################
# --mind 0.05 removes samples with >5% missing genotypes.
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step1 \
  --mind 0.05 \
  --make-pgen \
  --out snps_array_step2_mind


######################################
# 6. Remove heterozygosity outliers
######################################
# --het computes per-sample inbreeding coefficient F.
# Outliers will be identified in R (see R code below).
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step2_mind \
  --het \
  --out het

# After running the R script, het_outliers.txt will contain FID IID of outlier samples.
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step2_mind \
  --remove het_outliers.txt \
  --make-pgen \
  --out snps_array_step3_nohet


############################################
# 7. Remove related individuals (>0.0884)
############################################
# --king-cutoff 0.0884 removes one sample from each pair above the cutoff
# (roughly closer than 3rd-degree relatives).
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step3_nohet \
  --king-cutoff 0.0884 \
  --make-pgen \
  --out snps_array_step4_unrelated


####################
# 🟧 SNP-LEVEL QC  #
####################

#####################################################################
# 8. Filter SNPs by missingness <10%, MAF >5%, HWE p > 1e−6
#####################################################################
# --geno 0.10 : remove SNPs with >10% missingness
# --maf 0.05  : keep SNPs with minor allele frequency ≥5%
# --hwe 1e-6  : remove SNPs deviating strongly from HWE
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step4_unrelated \
  --geno 0.10 \
  --maf 0.05 \
  --hwe 1e-6 \
  --make-pgen \
  --out snps_array_step6_snpqc


#############################################################
# 9. Remove strand-ambiguous SNPs (A/T, T/A, C/G, G/C)
#############################################################
# Extract ambiguous SNP IDs from the .pvar file and save to ambiguous_snps.txt
awk 'NR>1 && (($4=="A" && $5=="T") || \
              ($4=="T" && $5=="A") || \
              ($4=="C" && $5=="G") || \
              ($4=="G" && $5=="C")) {print $3}' \
    snps_array_step6_snpqc.pvar > ambiguous_snps.txt

# Exclude ambiguous SNPs and create final SNP-QCed dataset
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step6_snpqc \
  --exclude ambiguous_snps.txt \
  --make-pgen \
  --out snps_array_step7_noambig


##################################
# 10. Count remaining SNPs
##################################
# Write all variant IDs to a list, then count how many SNPs passed QC.
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step7_noambig \
  --write-snplist \
  --out snps_array_step7_noambig

wc -l snps_array_step7_noambig.snplist
# (lines = number of SNPs that survived all QC steps)
