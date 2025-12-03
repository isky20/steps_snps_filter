#Adjust the header info HWp

sed 's/OLDWORD/NEWWORD/g' input.vcf > output.vcf

#Add extra ID to info header in vcf

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

bcftools annotate -h extra_info_header.hdr  -Oz -o vcf_with_info_header.vcf.gz input.vcf
tadix vcf_with_info_header.vcf.gz

#Create VCF available for plink 

bcftools annotate \
  -x INFO \
  -Oz -o vcf_noinfo.vcf.gz \
  vcf_with_info_header.vcf.gz

tabix  vcf_noinfo.vcf.gz

#Convert your SNP-array VCF to PLINK format (autosomal, bi-allelic SNPs)
singularity exec plink_combo.sif plink2
--vcf vcf_noinfo.vcf.gz
--snps-only just-acgt
--max-alleles 2
--make-pgen
--out snps_array_qc0
--allow-extra-chr
--threads 32
--split-par hg38


🟦 SAMPLE-LEVEL QC
1️⃣ autosome   
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_qc0   --autosome   --make-pgen   --out snps_array_step1 --allow-extra-chr
  
2️⃣ Remove samples with high missingness (> 5%)
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step1   --mind 0.05  --make-pgen   --out snps_array_step2_mind

3️⃣ Remove heterozygosity outliers
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step2_mind \
  --het \
  --out het
# Detect outliers in R  (|F − mean(F)| > 3SD), produce
#===============================================================
# Read het.het (do NOT treat '#' as comment)
het <- read.table("het.het",
                  header = TRUE,
                  comment.char = "",
                  stringsAsFactors = FALSE)

# Fix column name: X.IID -> IID
names(het)[1] <- "IID"

# Check
print(colnames(het))
# Should be: "IID" "O.HOM." "E.HOM." "OBS_CT" "F"

# Mean and SD of F
meanF <- mean(het$F, na.rm = TRUE)
sdF   <- sd(het$F,   na.rm = TRUE)

# Flag outliers: |F - mean(F)| > 3 SD
het$outlier <- abs(het$F - meanF) > 3 * sdF

# Subset outliers, keep as data.frame (drop = FALSE!)
out <- het[het$outlier, c("IID"), drop = FALSE]

if (nrow(out) > 0) {
  out2 <- data.frame(FID = out$IID, IID = out$IID)
  write.table(out2, "het_outliers.txt",
              col.names = FALSE,
              row.names  = FALSE,
              quote      = FALSE)
} else {
  # no outliers -> create empty file
  file.create("het_outliers.txt")
}
#===============================================================

singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step2_mind \
  --remove het_outliers.txt \
  --make-pgen \
  --out snps_array_step3_nohet

4️⃣ Remove related individuals (>0.0884)
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step3_nohet \
  --king-cutoff 0.0884 \
  --make-pgen \
  --out snps_array_step4_unrelated

🟧 SNP-LEVEL QC
6️⃣ Filter SNPs: missingness <10%, MAF>5%, HWE>1e−6
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step4_unrelated \
  --geno 0.10 \
  --maf 0.05 \
  --hwe 1e-6 \
  --make-pgen \
  --out snps_array_step6_snpqc
7️⃣ Remove ambiguous SNPs (A/T, T/A, C/G, G/C)
awk 'NR>1 && (($4=="A" && $5=="T") ||
  ($4=="T" && $5=="A") ||
  ($4=="C" && $5=="G") ||
  ($4=="G" && $5=="C")) {print $3}' snps_array_step6_snpqc.pvar > ambiguous_snps.txt

singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step6_snpqc \
  --exclude ambiguous_snps.txt \
  --make-pgen \
  --out snps_array_step7_noambig

#count
singularity exec plink_combo.sif plink2 \
  --pfile snps_array_step7_noambig \
  --write-snplist \
  --out snps_array_step7_noambig

wc -l snps_array_step7_noambig.snplist
