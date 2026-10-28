#!/bin/bash
# How the small teaching data in data/ were made (instructor use only - participants do NOT run this).
#
# Source: two real RNA-seq runs of ONE Holstein-Friesian calf (138 days old), GEO GSE289870 / ENA PRJNA1225222:
#     liver  = SRR32383427 (Liver_AI_Female_2,  GSM8800637)
#     muscle = SRR32383411 (Muscle_AI_Female_2, GSM8800653)
#   Illumina NovaSeq 6000, TruSeq stranded mRNA, 2 x 101 bp, about 23 million read pairs (2.3 GB) per run.
#
# What this script does (needs: curl, hisat2, hisat2-build, samtools, python3):
#   1. streams ONLY the first 4,000,000 read pairs of each run from ENA (about 200 MB per mate, not the full 2.3 GB)
#   2. downloads cattle chromosome 19 (Ensembl release 114, ARS-UCD2.0) and builds a HISAT2 index of it
#   3. aligns the 4 M pairs and keeps the pairs that map properly to chromosome 19
#   4. writes the first 100,000 of those pairs, UNCHANGED (same read names, bases and qualities), to data/reads/
#   5. cuts the chromosome-19 part out of the Ensembl GTF
# Run:  bash make_subset.sh  [output folder, default ./subset_build]

set -eu
OUT=${1:-subset_build}
N_HEAD=4000000      # read pairs streamed from the start of every run
N_KEEP=100000       # read pairs kept per sample
ENS=https://ftp.ensembl.org/pub/release-114
ENA=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR323
mkdir -p "$OUT/heads" "$OUT/ref" "$OUT/aln" "$OUT/final/reads" "$OUT/final/ref"
cd "$OUT"

# 1. first N_HEAD pairs of each run (zcat | head stops the download early)
for s in liver:SRR32383427:027 muscle:SRR32383411:011; do
  IFS=: read name run dir <<< "$s"
  for m in 1 2; do
    [ -s heads/${name}_R$m.fastq ] || curl -sS "$ENA/$dir/$run/${run}_$m.fastq.gz" | zcat 2>/dev/null | head -n $((N_HEAD*4)) > heads/${name}_R$m.fastq || true
  done
done

# 2. reference: chromosome 19 + annotation, HISAT2 index
[ -s ref/chr19.fa.gz ] || curl -sS -o ref/chr19.fa.gz $ENS/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD2.0.dna.primary_assembly.19.fa.gz
[ -s ref/ens114.gtf.gz ] || curl -sS -o ref/ens114.gtf.gz $ENS/gtf/bos_taurus/Bos_taurus.ARS-UCD2.0.114.gtf.gz
[ -s ref/chr19.fa ] || zcat ref/chr19.fa.gz | sed 's/^>\([^ ]*\) .*/>\1/' > ref/chr19.fa
[ -s ref/chr19.1.ht2 ] || hisat2-build -p 8 ref/chr19.fa ref/chr19 > ref/build.log 2>&1

# 3. align and list the properly paired read names (SRR number + read number, same in R1 and R2)
for name in liver muscle; do
  hisat2 -p 8 -x ref/chr19 -1 heads/${name}_R1.fastq -2 heads/${name}_R2.fastq 2> aln/$name.hisat2.log | samtools view -b -f 2 -F 256 - | samtools view - | cut -f1 | sort -u > aln/$name.names
  echo "$name: $(wc -l < aln/$name.names) read pairs map properly to chromosome 19"
done

# 4. keep the first N_KEEP of those pairs (in the original file order), unchanged
for name in liver muscle; do
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
    key = a[0].split()[0][1:]            # '@SRR32383427.1 A00125:...' -> 'SRR32383427.1'
    if key in keep:
        o1.writelines(a); o2.writelines(b); kept += 1
print(f"{name}: kept {kept} pairs")
EOF
done

# 5. annotation of chromosome 19 only (header lines + all lines of chromosome 19)
zcat ref/ens114.gtf.gz | awk -F'\t' '/^#/ || $1=="19"' | gzip -9 > final/ref/chr19.gtf.gz
cp ref/chr19.fa.gz final/ref/chr19.fa.gz
ls -l final/reads final/ref
