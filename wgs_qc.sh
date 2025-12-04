#Run sample-level QC with bcftools stats on that msVCF
bcftools stats -s - wgs_ms_raw.vcf.gz > wgs.stats
#extract total variant count with awk
grep '^PSC' wgs.stats \
  | awk 'BEGIN{
           OFS="\t";
           print "sample","nSingletons","TiTv","meanDP","totalVar";
         }
         {
           sample      = $3;
           nRefHom     = $4;
           nNonRefHom  = $5;
           nHets       = $6;
           nTs         = $7;
           nTv         = $8;
           nIndels     = $9;
           avgDP       = $10;
           nSingletons = $11;
           titv = (nTv > 0 ? nTs / nTv : 0);
           totalVar = nNonRefHom + nHets ;
           print sample, nSingletons, titv, avgDP, totalVar;
         }' > wgs_sample_qc_metrics.txt
#------extract good sample-----------
import numpy as np
# Input / output
metrics_file = "wgs_sample_qc_metrics.txt"   # columns: sample  nSingletons  TiTv  meanDP  totalVar
prefix       = "wgs"
# Read table (skip header)
samples = []
values  = []
with open(metrics_file) as f:
    header = next(f)
    for line in f:
        if not line.strip():
            continue
        flds = line.strip().split()
        samples.append(flds[0])
        values.append([float(flds[1]), float(flds[2]), float(flds[3]), float(flds[4])])
vals = np.array(values)   # shape (N, 4) = [nSingletons, TiTv, meanDP, totalVar]
# medians per column
med = np.median(vals, axis=0)
# absolute deviations from median
dev = np.abs(vals - med)
# MAD per column
mad = np.median(dev, axis=0)
# thresholds = 5 * MAD
thr = 5 * mad
# boolean mask: outlier if any metric > 5 * MAD from median
outlier_mask = (dev > thr).any(axis=1)
bad_samples  = [s for s, o in zip(samples, outlier_mask) if o]
good_samples = [s for s, o in zip(samples, outlier_mask) if not o]
with open(f"{prefix}_bad_samples.txt", "w") as fb:
    fb.write("\n".join(bad_samples) + "\n")
with open(f"{prefix}_good_samples.txt", "w") as fg:
    fg.write("\n".join(good_samples) + "\n")
#Filter out bad samples from the msVCF
bcftools view \
  -S wgs_good_samples.txt \
  -Oz -o wgs_samples_qc.vcf.gz \
  wgs_ms_raw.vcf.gz
tabix -p vcf wgs_samples_qc.vcf.gz
bcftools view -Ou wgs_samples_qc.vcf.gz \
  | bcftools +setGT -Ou -- -t q -n . -i 'FMT/DP<8 || FMT/GQ<20' \
  | bcftools +fill-tags -Ou -- -t F_MISSING \
  | bcftools view -Oz -o wgs_var_qc.vcf.gz \
      -i 'QUAL>=30 && F_MISSING<0.20 && INFO/gnomAD_FILTER="PASS"'
plink2 \
  --vcf wgs_var_qc.vcf.gz \
  --snps-only just-acgt \
  --max-alleles 2 \
  --geno 0.20 \
  --maf 0.05 \
  --hwe 5e-5 midp \
  --make-pgen \
  --out wgs_ps_step1
