#!/usr/bin/env bash
set -euo pipefail

########################################
# Config
########################################

# Raw WGS multi-sample VCF (msVCF)
RAW_VCF="wgs_ms_raw.vcf.gz"

# Prefix for output files
PREFIX="wgs"

echo "[INFO] Raw WGS VCF: ${RAW_VCF}"
echo "[INFO] Prefix      : ${PREFIX}"

########################################
# STEP 1 – Sample-level QC with bcftools
########################################
echo "[STEP 1] Running bcftools stats..."

bcftools norm -m-any -f GRCh38.fa \
  -Oz -o "${PREFIX}_samples_qc_norm.vcf.gz" \
  "${RAW_VCF}"

tabix -p vcf "${PREFIX}_samples_qc_norm.vcf.gz"


# Per-sample stats from raw WGS VCF
bcftools stats -s - "${PREFIX}_samples_qc_norm.vcf.gz" > "${PREFIX}.stats"

echo "[STEP 1] Extracting per-sample QC metrics (singletons, Ti/Tv, depth, total variants)..."

# Extract:
#  sample, nSingletons, TiTv, meanDP, totalVar
grep '^PSC' "${PREFIX}.stats" \
  | awk 'BEGIN{
           OFS="\t";
           print "sample","nSingletons","TiTv","meanDP","totalVar";
         }
         {
           # Columns from PSC line:
           # [3]sample [4]nRefHom [5]nNonRefHom [6]nHets
           # [7]nTransitions [8]nTransversions [9]nIndels
           # [10]average depth [11]nSingletons ...

           sample      = $3;
           nNonRefHom  = $5 + 0;
           nHets       = $6 + 0;
           nTs         = $7 + 0;
           nTv         = $8 + 0;
           avgDP       = $10 + 0;
           nSingletons = $11 + 0;

           # Total non-reference variants
           totalVar = nNonRefHom + nHets;

           # Ti/Tv: NA if there are no transversions
           if (nTv > 0) {
             titv = nTs / nTv;
           } else {
             titv = "NA";
           }

           print sample, nSingletons, titv, avgDP, totalVar;
         }' > "${PREFIX}_sample_qc_metrics.txt"


########################################
# STEP 2 – Python MAD filter (5×MAD)
########################################
echo "[STEP 2] Running Python MAD-based outlier detection..."

python3 - << 'EOF_PY'
import numpy as np

metrics_file = "wgs_sample_qc_metrics.txt"  # produced above
prefix       = "wgs"

samples = []
values  = []

with open(metrics_file) as f:
    header = next(f)  # skip header
    for line in f:
        if not line.strip():
            continue
        flds = line.strip().split()
        # Expect: sample nSingletons TiTv meanDP totalVar
        samples.append(flds[0])
        values.append([
            float(flds[1]),  # nSingletons
            float(flds[2]),  # TiTv
            float(flds[3]),  # meanDP
            float(flds[4])   # totalVar
        ])

vals = np.array(values)  # shape: (N, 4)

# Median per metric
med = np.median(vals, axis=0)

# Absolute deviation from median
dev = np.abs(vals - med)

# Median absolute deviation (MAD) per metric
mad = np.median(dev, axis=0)

# Threshold = 5 × MAD
thr = 5 * mad

# Mark as outlier if any metric deviates > 5×MAD
outlier_mask = (dev > thr).any(axis=1)

bad_samples  = [s for s, o in zip(samples, outlier_mask) if o]
good_samples = [s for s, o in zip(samples, outlier_mask) if not o]

print(f"[PYTHON] N samples: {len(samples)}")
print(f"[PYTHON] N bad   : {len(bad_samples)}")
print(f"[PYTHON] N good  : {len(good_samples)}")

with open(f"{prefix}_bad_samples.txt", "w") as fb:
    if bad_samples:
        fb.write("\n".join(bad_samples) + "\n")

with open(f"{prefix}_good_samples.txt", "w") as fg:
    if good_samples:
        fg.write("\n".join(good_samples) + "\n")
EOF_PY

########################################
# STEP 3 – Keep only good samples in msVCF
########################################
echo "[STEP 3] Subsetting VCF to good samples..."

bcftools view \
  -S "${PREFIX}_good_samples.txt" \
  -Oz -o "${PREFIX}_samples_qc.vcf.gz" \
  "${PREFIX}_samples_qc_norm.vcf.gz"

tabix -p vcf "${PREFIX}_samples_qc.vcf.gz"

########################################
# STEP 4 – Genotype & site-level QC
########################################
echo "[STEP 4] Applying genotype & variant filters (DP,GQ,QUAL,F_MISSING)..."

bcftools view -Ou "${PREFIX}_samples_qc.vcf.gz" \
  | bcftools +setGT -Ou -- -t q -n . -i 'FMT/GQ<20' \
  | bcftools +fill-tags -Ou -- -t F_MISSING \
  | bcftools view -Oz -o "${PREFIX}_var_qc.vcf.gz" \
      -i 'QUAL>=30 && F_MISSING<0.20'

tabix -p vcf "${PREFIX}_var_qc.vcf.gz"

########################################
# STEP 5 – Convert to PLINK2 pgen + SNP QC
########################################
echo "[STEP 5] Converting to PLINK2 pgen and applying SNP-level filters..."

plink2 \
  --vcf "${PREFIX}_var_qc.vcf.gz" \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno 0.20 \
  --maf 0.05 \
  --hwe 5e-5 midp \
  --make-pgen \
  --out "${PREFIX}_ps_step1"

echo "[DONE] WGS QC finished. Final PLINK2 dataset: ${PREFIX}_ps_step1"

############################################
# STEP 6 – Remove related individuals (>0.0884)
############################################
echo "[STEP 6] Removing related individuals (KING cutoff 0.0884)..."

plink2 \
  --pfile "${PREFIX}_ps_step1" \
  --king-cutoff 0.0884 \
  --make-pgen \
  --out "${PREFIX}_ps_step2_unrelated"

echo "[DONE] Final unrelated WGS dataset: ${PREFIX}_ps_step2_unrelated"


plink2 \
  --pfile "${PREFIX}_ps_step2_unrelated" \
  --write-snplist \
  --out WGS_snps

wc -l WGS_snps.snplist
