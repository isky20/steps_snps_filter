# with plink2 
plink2 \
  --pfile snps_array_step7_noambig \
  --extract wgs_ps_step1.snplist \
  --write-snplist \
  --out snps_array_wgs_overlap

wc -l snps_array_wgs_overlap.snplist

# Array-only SNPs = all array SNPs minus overlap
plink2 \
  --pfile snps_array_step7_noambig \
  --exclude snps_array_wgs_overlap.snplist \
  --write-snplist \
  --out snps_array_only_not_in_wgs

  wc -l snps_array_only_not_in_wgs


# All SNPs from WGS
plink2 \
  --pfile wgs_ps_step1 \
  --write-snplist \
  --out wgs_all

wc -l  wgs_all



# WGS-only SNPs = all WGS SNPs minus overlap
plink2 \
  --pfile wgs_ps_step1 \
  --exclude snps_array_wgs_overlap.snplist \
  --write-snplist \
  --out wgs_only_not_in_array

wc -l  wgs_only_not_in_array
