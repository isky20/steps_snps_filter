# Easy way with plink2 
plink2 \
  --pfile snps_array_step7_noambig \
  --extract wgs_ps_step1.snplist \
  --write-snplist \
  --out snps_array_wgs_overlap

wc -l snps_array_wgs_overlap.snplist

#===================by AF=======================
# SNP array AF
plink2 \
  --pfile snps_array_step7_noambig \
  --extract shared_snps.list \
  --freq \
  --out snps_array_shared

# WGS AF
plink2 \
  --pfile wgs_ps_step1 \
  --extract shared_snps.list \
  --freq \
  --out wgs_shared

#================================================
import numpy as np

def load_afreq(path):
    data = {}
    with open(path) as f:
        header = next(f)
        for line in f:
            if not line.strip():
                continue
            c = line.strip().split()
            snp = c[1]
            ref = c[2]
            alt = c[3]
            af  = float(c[4])   # ALT allele frequency
            data[snp] = (ref, alt, af)
    return data

array_af = load_afreq("snps_array_shared.afreq")
wgs_af   = load_afreq("wgs_shared.afreq")

shared_snps = sorted(set(array_af.keys()) & set(wgs_af.keys()))

freq_diffs = []
mismatch_snps = []
allele_mismatch_snps = []

# threshold for "bad" allele frequency mismatch (you can change 0.20 to 0.10 etc.)
THR = 0.20

for snp in shared_snps:
    ref_a, alt_a, af_a = array_af[snp]
    ref_w, alt_w, af_w = wgs_af[snp]

    # Same REF/ALT
    if (ref_a == ref_w) and (alt_a == alt_w):
        diff = abs(af_a - af_w)
        freq_diffs.append(diff)
        if diff > THR:
            mismatch_snps.append((snp, ref_a, alt_a, af_a, af_w, diff))

    # Swapped REF/ALT (e.g. array ALT == WGS REF)
    elif (ref_a == alt_w) and (alt_a == ref_w):
        # WGS ALT freq corresponds to (1 - WGS ALT) if alleles flipped
        diff = abs(af_a - (1.0 - af_w))
        freq_diffs.append(diff)
        if diff > THR:
            mismatch_snps.append((snp, ref_a, alt_a, af_a, 1.0 - af_w, diff))

    else:
        # Alleles don’t match at all (likely bad mapping / different build)
        allele_mismatch_snps.append((snp, ref_a, alt_a, ref_w, alt_w))

# Write lists
with open("freq_mismatch_snps.txt", "w") as out:
    out.write("SNP\tREF\tALT\tAF_array\tAF_wgs_aligned\tabs_diff\n")
    for snp, ref, alt, af_a, af_w_aligned, d in mismatch_snps:
        out.write(f"{snp}\t{ref}\t{alt}\t{af_a:.4f}\t{af_w_aligned:.4f}\t{d:.4f}\n")

with open("allele_mismatch_snps.txt", "w") as out:
    out.write("SNP\tREF_array\tALT_array\tREF_wgs\tALT_wgs\n")
    for snp, ra, aa, rw, aw in allele_mismatch_snps:
        out.write(f"{snp}\t{ra}\t{aa}\t{rw}\t{aw}\n")

# Quick summary
if freq_diffs:
    freq_diffs = np.array(freq_diffs)
    print("Number of shared SNPs compared:", len(freq_diffs))
    print("Mean |AF_array - AF_wgs|:", float(freq_diffs.mean()))
    print("Median |AF_array - AF_wgs|:", float(np.median(freq_diffs)))
    print("Number with |ΔAF| >", THR, ":", len(mismatch_snps))
else:
    print("No shared SNPs to compare.")
#================================================

cut -f1 freq_mismatch_snps.txt | tail -n +2 > mismatched_freq_snps.txt

plink2 \
  --pfile snps_array_step7_noambig \
  --exclude mismatched_freq_snps.txt \
  --make-pgen \
  --out snps_array_step8_freqok


wc -l mismatched_freq_snps.txt

plink2 \
  --pfile snps_array_step8_freqok \
  --write-snplist \
  --out snps_array_step8_freqok

wc -l snps_array_step8_freqok.snplist
