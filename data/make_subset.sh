#!/bin/bash
# How the small teaching data in data/ were made (instructor use only - participants do NOT run this).
#
# Source: three real RNA-seq runs of LIVER from three Danish Landrace boars (one per boar-taint group),
#   GEO GSE113316 / ENA PRJNA451072 (the same animals as the Day 2 R tutorials):
#     high   = animal 7107, SRR7041911 (GSM2461212, "Sample 7107 - Liver - high boar taint")
#     medium = animal 2111, SRR7041924 (GSM2461225, "Sample 2111 - Liver - medium boar taint")
#     low    = animal 1105, SRR7041939 (GSM2461240, "Sample 1105 - Liver - low boar taint")
#   Illumina HiSeq 2500, 2 x 100 bp, 15 to 20 million read pairs (1.4 to 1.7 GB) per run.
#
# What this script does (needs: curl, hisat2, hisat2-build, samtools, python3):
#   1. streams ONLY the first 8,000,000 read pairs of each run from ENA (about 350 MB per mate, not the full run)
#   2. downloads pig chromosome 12 (Ensembl release 79, Sscrofa10.2) and builds a HISAT2 index of it
#   3. aligns the 8 M pairs and keeps the pairs that map properly to chromosome 12
#   4. writes the first 100,000 of those pairs, UNCHANGED (same read names, bases and qualities), to data/reads/
#   5. cuts the chromosome-12 part out of the Ensembl release 79 GTF
# Run:  bash make_subset.sh  [output folder, default ./subset_build]

set -eu
OUT=${1:-subset_build}
N_HEAD=8000000      # read pairs streamed from the start of every run
N_KEEP=100000       # read pairs kept per sample
CHR=12
ENS=https://ftp.ensembl.org/pub/release-79
ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR704
mkdir -p "$OUT/heads" "$OUT/ref" "$OUT/aln" "$OUT/final/reads" "$OUT/final/ref"
cd "$OUT"

# 1. first N_HEAD pairs of each run (zcat | head stops the download early)
for s in high_7107:SRR7041911:001 medium_2111:SRR7041924:004 low_1105:SRR7041939:009; do
  IFS=: read name run dir <<< "$s"
  for m in 1 2; do
    [ -s heads/${name}_R$m.fastq ] || curl -sS "$ENA/$dir/$run/${run}_$m.fastq.gz" | zcat 2>/dev/null | head -n $((N_HEAD*4)) > heads/${name}_R$m.fastq || true
  done
done

# 2. reference: chromosome 12 + annotation, HISAT2 index
[ -s ref/chr$CHR.fa.gz ] || curl -sS -o ref/chr$CHR.fa.gz $ENS/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa10.2.dna.chromosome.$CHR.fa.gz
[ -s ref/ens79.gtf.gz ] || curl -sS -o ref/ens79.gtf.gz $ENS/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.79.gtf.gz
[ -s ref/chr$CHR.fa ] || zcat ref/chr$CHR.fa.gz | sed 's/^>\([^ ]*\) .*/>\1/' > ref/chr$CHR.fa
[ -s ref/chr$CHR.1.ht2 ] || hisat2-build -p 8 ref/chr$CHR.fa ref/chr$CHR > ref/build.log 2>&1

# 3. align and list the properly paired read names (SRR number + read number, same in R1 and R2)
for name in high_7107 medium_2111 low_1105; do
  hisat2 -p 8 -x ref/chr$CHR -1 heads/${name}_R1.fastq -2 heads/${name}_R2.fastq 2> aln/$name.hisat2.log | samtools view -b -f 2 -F 256 - | samtools view - | cut -f1 | sort -u > aln/$name.names
  echo "$name: $(wc -l < aln/$name.names) read pairs map properly to chromosome $CHR"
done

# 4. keep the first N_KEEP of those pairs (in the original file order), unchanged
for name in high_7107 medium_2111 low_1105; do
  python3 - "$name" "$N_KEEP" <<'EOF'
import sys, gzip
name, n_keep = sys.argv[1], int(sys.argv[2])
keep = set(l.strip() for l in open(f"aln/{name}.names"))
r1, r2 = open(f"heads/{name}_R1.fastq"), open(f"heads/{name}_R2.fastq")
o1 = gzip.open(f"final/reads/{name}_R1.fastq.gz", "wt", compresslevel=9)
o2 = gzip.open(f"final/reads/{name}_R2.fastq.gz", "wt", compresslevel=9)
kept = 0
while kept < n_keep:
    a = [r1.readline() for _ in range(4)]
    b = [r2.readline() for _ in range(4)]
    if not a[0]:
        break
    key = a[0].split()[0][1:]            # '@SRR7041911.1 HWI-...' -> 'SRR7041911.1'
    if key in keep:
        o1.writelines(a); o2.writelines(b); kept += 1
print(f"{name}: kept {kept} pairs")
EOF
done

# 5. annotation of chromosome 12 only (header lines + all lines of chromosome 12)
zcat ref/ens79.gtf.gz | awk -F'\t' -v c=$CHR '/^#/ || $1==c' | gzip -9 > final/ref/chr$CHR.gtf.gz
cp ref/chr$CHR.fa.gz final/ref/chr$CHR.fa.gz
ls -l final/reads final/ref
