# GenART Day 2 - RNA-seq: from reads (FASTQ) to a gene count matrix (shell commands)
### Boar taint in pigs: one liver sample from each group (High, Medium, Low)

Hands-on tutorial for the GenART workshop (TANUVAS, Chennai): real RNA-seq reads of the **liver of three Danish Landrace boars**
(one with **High**, one with **Medium**, one with **Low** boar taint; the same animals as in the Day 2 R tutorials) go through
FastQC, MultiQC, fastp, HISAT2, samtools and featureCounts to a **gene count matrix**. Every step is a plain shell command; each command is printed and its full output is shown.

| Part | What |
|---|---|
| **A** | Look at the files: FASTQ reads, FASTA reference, GTF annotation |
| **B** | 16 steps: FastQC - MultiQC - fastp - FastQC again - HISAT2 index - HISAT2 alignment - BAM files - alignment QC - library strandedness - featureCounts - gene count matrix - look at the matrix |
| **C** | 1 figure: what happened to the reads (HISAT2 alignment and featureCounts assignment) |

## Run it in Google Colab (nothing to install)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pratikpathade/GenartWorkshop_Pratik/blob/main/Day2_RNAseq_boar_taint_FASTQ_to_counts_shell_tutorial.ipynb)

1. Click the badge, then **Connect** (top right).
2. Run the cells from top to bottom (play button, or *Runtime > Run all*). Tested end to end in 4-7 minutes on 2 CPU cores, including the tool install; a free Colab session can be slower.

The notebook `Day2_RNAseq_boar_taint_FASTQ_to_counts_shell_tutorial_with_outputs.ipynb` is the same notebook with all outputs saved from a test run (what you should see) - handy for reading without running or for comparing your own output.

## Run it as a script (Linux / Colab terminal)

```bash
git clone https://github.com/PratikPathade/GenartWorkshop_Pratik.git
cd GenartWorkshop_Pratik
bash Day2_RNAseq_boar_taint_FASTQ_to_counts_shell_tutorial.sh
```

Everything is made in `./RNAseq_pig`; the tools are installed to `~/genart-rnaseq` (micromamba: FastQC 0.13.0, MultiQC 1.35, fastp 1.4.0, HISAT2 2.2.3, samtools 1.24, Subread/featureCounts 2.1.1).

## Data (`data/`, about 40 MB)

| Files | What |
|---|---|
| `data/reads/high_7107, medium_2111, low_1105` (R1 and R2, FASTQ) | Real liver reads of three Danish Landrace boars (Illumina HiSeq 2500, 2 x 100 bases): animal 7107 (High boar taint), 2111 (Medium), 1105 (Low); 100,000 read pairs each, only pairs from **chromosome 12** |
| `data/ref/chr12.fa.gz`, `chr12.gtf.gz` | Pig chromosome 12 (Sscrofa10.2) and its gene annotation, both from Ensembl release 79 (the release of the Day 2 R tutorials) |
| `data/make_subset.sh` | How the small reads were cut from the original 1.4-1.7 GB runs (instructor use) |

Sources and what was changed: `data/SOURCES.md`.

With **one animal per group** this tutorial teaches the pipeline, not the biology of boar taint. Differential expression needs replicates: the Day 2 R tutorials continue with the counts of all nine boars (three per group).

**Please cite** GEO series GSE113316 / ENA project PRJNA451072 for the reads (Drag et al. 2017, *Scientific Reports*, doi:10.1038/s41598-017-11928-0) and Ensembl release 79 for the reference.
