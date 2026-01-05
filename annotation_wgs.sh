#1) Find Nirvana folder (v4.3)
ls -lah /opt/dragen/4.3*/share/nirvana/

#2) Create the data folder
mkdir -p /data/nirvana

#3) Download Nirvana annotation bundle (GRCh38)
/opt/dragen/4.3*/share/nirvana/Downloader --ga GRCh38 --out /data/nirvana

#4) Prepare your input VCF (bgzip + tabix index)
bgzip -c cohort.vcf > cohort.vcf.gz
tabix -p vcf cohort.vcf.gz

#5) Run Nirvana (standalone, GRCh38)
/opt/dragen/4.3*/share/nirvana/Nirvana \
  -c /data/nirvana/Cache/GRCh38/Both \
  --sd /data/nirvana/SupplementaryAnnotation/GRCh38 \
  -r /data/nirvana/References/Homo_sapiens.GRCh38.Nirvana.dat \
  -i cohort.vcf.gz \
  -o cohort.nirvana

