# Easy way with plink2 
plink2 \
  --pfile snps_array_step7_noambig \
  --extract wgs_ps_step1.snplist \
  --write-snplist \
  --out snps_array_wgs_overlap

wc -l snps_array_wgs_overlap.snplist

