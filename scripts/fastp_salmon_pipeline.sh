#!/bin/bash
set -euo pipefail

THREADS=8
SRRS=(SRR5219989 SRR5219985 SRR5220016 SRR5219984)

# Build Salmon index once (if not already present)
if [ ! -d "salmon_index" ]; then
  salmon index -t gencode.v49.transcripts.fa.gz -i salmon_index -p $THREADS
fi

for SRR in "${SRRS[@]}"; do
  echo "=== Processing $SRR ==="

  # Download SRA
  prefetch "$SRR" --progress --transport https

  # Convert to FASTQ (paired-end)
  fasterq-dump "$SRR.sra" --split-files --threads $THREADS

  # Compress FASTQ (optional but recommended)
  gzip -f "${SRR}_1.fastq"
  gzip -f "${SRR}_2.fastq"

  # Trimming + QC
  fastp \
    -i "${SRR}_1.fastq.gz" \
    -I "${SRR}_2.fastq.gz" \
    -o "Trim${SRR}_1.fastq.gz" \
    -O "Trim${SRR}_2.fastq.gz" \
    --detect_adapter_for_pe \
    --thread $THREADS \
    --html "${SRR}_fastp.html" \
    --json "${SRR}_fastp.json" \
    --report_title "Fastp QC report for ${SRR}"

  # Salmon quantification
  salmon quant -i salmon_index -l A \
    -1 "Trim${SRR}_1.fastq.gz" \
    -2 "Trim${SRR}_2.fastq.gz" \
    -p $THREADS \
    -o "salmon_${SRR}"

done

echo "All samples processed."