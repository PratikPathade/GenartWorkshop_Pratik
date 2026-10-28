#!/bin/bash
# GenART Day 2 - RNA-seq: from reads (FASTQ) to a gene count matrix. The shell cells of the notebook, in order.
# Run:  git clone https://github.com/PratikPathade/GenartWorkshop_Pratik.git && cd GenartWorkshop_Pratik && bash Day2_RNAseq_FASTQ_to_counts_shell_tutorial.sh
set -e
export LAB=$HOME/genart-rnaseq REPO=https://github.com/PratikPathade/GenartWorkshop_Pratik.git
export PATH=$LAB/env/bin:$LAB/bin:$PATH MAMBA_ROOT_PREFIX=$LAB/mamba
export THREADS=$(nproc)
trap 'case "$BASH_COMMAND" in "for "*|"[ "*) ;; *) printf "\n\033[1;32m$ %s\033[0m\n" "$BASH_COMMAND";; esac' DEBUG
ROOT=$PWD; mkdir -p RNAseq_calf; cd RNAseq_calf


# ====================================================================================================
# GenART Day 2 - RNA-seq: from reads (FASTQ) to a gene count matrix, all in shell commands


# ====================================================================================================
# 0. Setup

rm -rf data
cp -r $ROOT/data data
ls -lh data/reads data/ref

mkdir -p $LAB/bin
curl -fsSL -o $LAB/bin/micromamba https://github.com/mamba-org/micromamba-releases/releases/download/2.9.0-0/micromamba-linux-64
chmod +x $LAB/bin/micromamba
micromamba --version

[ -x $LAB/env/bin/hisat2 ] || micromamba create -y -p $LAB/env -c conda-forge -c bioconda python=3.11 fastqc=0.13.0 multiqc=1.35 fastp=1.4.0 hisat2=2.2.3 samtools=1.24 subread=2.1.1

fastqc --version
multiqc --version
fastp --version
hisat2 --version | head -n 1
samtools --version | head -n 1
featureCounts -v 2>&1 | grep featureCounts

mkdir -p ref index results/qc_raw results/trim results/qc_trim results/align results/counts
ls


# ====================================================================================================
# PART A - Look at the files

zcat data/reads/liver_R1.fastq.gz | head -n 8

zcat data/reads/liver_R2.fastq.gz | head -n 4

ls -lh data/reads
for f in data/reads/*.fastq.gz; do echo "$f   $(zcat $f | wc -l | awk '{print $1/4}') reads"; done

zcat data/reads/liver_R1.fastq.gz | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c

zcat data/reads/liver_R1.fastq.gz | head -n 40000 | awk 'NR%4==0' | fold -w1 | sort | uniq -c | sort -k1,1nr
for c in F : , '#'; do echo "symbol $c  = Q$(( $(printf '%d' "'$c") - 33 ))"; done


# ====================================================================================================
# Step 2 of 16: FASTA - the reference genome

gunzip -c data/ref/chr19.fa.gz > ref/chr19.fa
ls -lh ref
grep '>' ref/chr19.fa

samtools faidx ref/chr19.fa
cat ref/chr19.fa.fai
samtools faidx ref/chr19.fa 19:20000001-20000090


# ====================================================================================================
# Step 3 of 16: GTF - the gene annotation

zcat data/ref/chr19.gtf.gz | head -n 5
zcat data/ref/chr19.gtf.gz | grep -v '^#' | head -n 3 | cut -c1-260

zcat data/ref/chr19.gtf.gz | grep -v '^#' | cut -f3 | sort | uniq -c | sort -k1,1nr
zcat data/ref/chr19.gtf.gz | awk -F'\t' '$3=="gene"' | grep -o 'gene_biotype "[^"]*"' | sort | uniq -c | sort -k1,1nr

zcat data/ref/chr19.gtf.gz | grep 'gene_name "VTN"' | awk -F'\t' '$3=="gene" || $3=="exon"' | cut -f1,3,4,5,7 | head -n 8


# ====================================================================================================
# PART B - From reads to a gene count matrix

fastqc -t $THREADS --extract --outdir results/qc_raw data/reads/*.fastq.gz
ls results/qc_raw


# ====================================================================================================
# Step 5 of 16: Read the FastQC result

ls results/qc_raw/liver_R1_fastqc
cat results/qc_raw/{liver_R1,liver_R2,muscle_R1,muscle_R2}_fastqc/summary.txt | awk -F'\t' '{m[$2]=m[$2]"   "$1; if(!($2 in o)){o[$2]=++n; name[n]=$2}} END{printf "%-34s %s\n","module","liver_R1 liver_R2 muscle_R1 muscle_R2"; for(i=1;i<=n;i++) printf "%-34s %s\n", name[i], m[name[i]]}'

sed -n '1,12p' results/qc_raw/liver_R1_fastqc/fastqc_data.txt

sed -n '/>>Per base sequence quality/,/>>END_MODULE/p' results/qc_raw/liver_R1_fastqc/fastqc_data.txt | sed -n '1,14p'


# ====================================================================================================
# Step 6 of 16: MultiQC - all reports in one

multiqc --no-ansi results/qc_raw -o results/multiqc_raw -n multiqc_raw --force
ls results/multiqc_raw results/multiqc_raw/multiqc_raw_data

column -t -s $'\t' results/multiqc_raw/multiqc_raw_data/multiqc_general_stats.txt


# ====================================================================================================
# Step 7 of 16: fastp - trim the reads

for s in liver muscle; do
  echo "=================== sample: $s"
  fastp -i data/reads/${s}_R1.fastq.gz -I data/reads/${s}_R2.fastq.gz \
        -o results/trim/${s}_R1.fastq.gz -O results/trim/${s}_R2.fastq.gz \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --cut_tail --cut_tail_mean_quality 20 --length_required 36 \
        --thread $THREADS --json results/trim/${s}.fastp.json --html results/trim/${s}.fastp.html --report_title $s
done

ls -lh results/trim


# ====================================================================================================
# Step 8 of 16: FastQC and MultiQC again - did trimming help?

fastqc -t $THREADS --extract --outdir results/qc_trim results/trim/*.fastq.gz
multiqc --no-ansi results/qc_raw results/qc_trim results/trim -o results/multiqc_trim -n multiqc_raw_and_trimmed --force

cat results/qc_raw/liver_R1_fastqc/summary.txt results/qc_trim/liver_R1_fastqc/summary.txt results/qc_raw/muscle_R1_fastqc/summary.txt results/qc_trim/muscle_R1_fastqc/summary.txt | awk -F'\t' '{m[$2]=m[$2]"   "$1; if(!($2 in o)){o[$2]=++n; name[n]=$2}} END{printf "%-34s %s\n","module","liver_raw liver_trim muscle_raw muscle_trim"; for(i=1;i<=n;i++) printf "%-34s %s\n", name[i], m[name[i]]}'


# ====================================================================================================
# Step 9 of 16: HISAT2 - build the genome index

hisat2-build -p $THREADS ref/chr19.fa index/chr19 > index/hisat2-build.log 2>&1
head -n 8 index/hisat2-build.log
tail -n 6 index/hisat2-build.log
ls -lh index


# ====================================================================================================
# Step 10 of 16: HISAT2 - align the reads

for s in liver muscle; do
  echo "=================== sample: $s"
  hisat2 -p $THREADS -x index/chr19 -1 results/trim/${s}_R1.fastq.gz -2 results/trim/${s}_R2.fastq.gz 2> results/align/${s}.hisat2.log | samtools sort -@ 2 -o results/align/${s}.bam -
  samtools index results/align/${s}.bam
done
ls -lh results/align

cat results/align/liver.hisat2.log
cat results/align/muscle.hisat2.log


# ====================================================================================================
# Step 11 of 16: Look inside the BAM file

samtools view -H results/align/liver.bam
samtools view results/align/liver.bam | head -n 3 | cut -f1-9,12-

samtools flags 99 147

samtools view -c results/align/liver.bam
samtools view results/align/liver.bam | awk '$6 ~ /N/' | wc -l
samtools view results/align/liver.bam | awk '$5==60 && $9>0 && $9<5000 && $6 ~ /^[0-9]+M[0-9][0-9][0-9][0-9]?N[0-9]+M$/' | head -n 4 | cut -f1-9


# ====================================================================================================
# Step 12 of 16: Alignment QC

for s in liver muscle; do
  samtools flagstat results/align/$s.bam > results/align/$s.flagstat.txt
  samtools stats results/align/$s.bam > results/align/$s.stats.txt
done
cat results/align/liver.flagstat.txt
cat results/align/muscle.flagstat.txt

grep ^SN results/align/liver.stats.txt | cut -f2,3 | sed -n '1,8p;12,13p;15,17p;20,22p'


# ====================================================================================================
# Step 13 of 16: Which strand? (library strandedness)

for s in 0 1 2; do
  featureCounts -T $THREADS -p --countReadPairs -s $s -t exon -g gene_id -a data/ref/chr19.gtf.gz -o results/counts/test_s$s.txt results/align/liver.bam > /dev/null 2>&1
  echo "-s $s"; grep -E '^(Assigned|Unassigned_NoFeatures|Unassigned_Ambiguity)' results/counts/test_s$s.txt.summary
done


# ====================================================================================================
# Step 14 of 16: featureCounts - count read pairs per gene

featureCounts -T $THREADS -p --countReadPairs -s 2 -t exon -g gene_id --extraAttributes gene_name -a data/ref/chr19.gtf.gz -o results/counts/counts.txt results/align/liver.bam results/align/muscle.bam

cat results/counts/counts.txt.summary

multiqc --no-ansi results/qc_raw results/qc_trim results/trim results/align results/counts/counts.txt.summary -o results/multiqc_final -n multiqc_final --force


# ====================================================================================================
# Step 15 of 16: Make the gene count matrix

head -n 5 results/counts/counts.txt | cut -f1,2,6-

grep -v '^#' results/counts/counts.txt | cut -f1,8- | sed '1s#results/align/##g; 1s#\.bam##g' > results/gene_count_matrix.tsv
grep -v '^#' results/counts/counts.txt | awk -F'\t' -v OFS='\t' 'NR==1{print "gene_id","gene_name","length","liver","muscle"; next} {print $1, ($7=="NA"?$1:$7), $6, $8, $9}' > results/gene_count_matrix_annotated.tsv
head -n 6 results/gene_count_matrix.tsv
head -n 6 results/gene_count_matrix_annotated.tsv


# ====================================================================================================
# Step 16 of 16: Look at the matrix

echo "genes (rows): $(( $(wc -l < results/gene_count_matrix.tsv) - 1 ))"
awk -F'\t' 'NR>1 {l+=$2; m+=$3; if($2>0) gl++; if($3>0) gm++; if($2>0 && $3>0) both++; if($2+$3>=10) ten++} END{print "counted read pairs  liver:", l, " muscle:", m; print "genes with counts   liver:", gl, " muscle:", gm, " both:", both; print "genes with >= 10 counts in total:", ten}' results/gene_count_matrix.tsv

echo "--- liver"
sort -t$'\t' -k4,4nr results/gene_count_matrix_annotated.tsv | head -n 8 | column -t -s $'\t'
echo "--- muscle"
sort -t$'\t' -k5,5nr results/gene_count_matrix_annotated.tsv | head -n 8 | column -t -s $'\t'

head -n 1 results/gene_count_matrix_annotated.tsv
grep -w -E 'VTN|APOH|ASGR1|G6PC1|MYH1|MYH2|ENO3|SCN4A|CACNG1' results/gene_count_matrix_annotated.tsv | sort -t$'\t' -k2,2 | column -t -s $'\t'

awk -F'\t' -v OFS='\t' 'NR==FNR{ if(FNR>1){l+=$4; m+=$5} next } FNR==1{print; next} {printf "%s\t%s\t%s\t%.1f\t%.1f\n", $1,$2,$3,$4*1e6/l,$5*1e6/m}' results/gene_count_matrix_annotated.tsv results/gene_count_matrix_annotated.tsv > results/gene_cpm_annotated.tsv
head -n 5 results/gene_cpm_annotated.tsv | column -t -s $'\t'


# ====================================================================================================
# PART C - Figures


# ====================================================================================================
# PART D - Wrap-up

ls -lh results
du -sh results/*  | sort -h


# ====================================================================================================
# Questions
