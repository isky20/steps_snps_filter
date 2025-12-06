sort snps_array.snplist > snps_array.sorted
sort WGS_snps.snplist   > WGS_snps.sorted


# 1) Overlap: present in BOTH snps_array AND WGS
comm -12 snps_array.sorted WGS_snps.sorted > snps_array_wgs_overlap.snplist

# 2) snps_array only: in array, NOT in WGS
comm -23 snps_array.sorted WGS_snps.sorted > snps_array_only.snplist

# 3) WGS only: in WGS, NOT in array
comm -13 snps_array.sorted WGS_snps.sorted > WGS_only.snplist



wc -l snps_array_wgs_overlap.snplist
wc -l snps_array_only.snplist
wc -l WGS_only.snplist
