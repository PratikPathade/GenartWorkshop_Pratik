# GenART Day 2 - RNA-seq: from reads (FASTQ) to a gene count matrix (shell commands)

Hands-on tutorial for the GenART workshop (TANUVAS, Chennai): real RNA-seq reads of a **calf liver** and a **calf muscle** sample go through
FastQC, MultiQC, fastp, HISAT2, samtools and featureCounts to a **gene count matrix**. Every step is a plain shell command; each command is printed and its full output is shown.

| Part | What |
|---|---|
| **A** | Look at the files: FASTQ reads, FASTA reference, GTF annotation |
| **B** | 16 steps: FastQC - MultiQC - fastp - FastQC again - HISAT2 index - HISAT2 alignment - BAM files - alignment QC - library strandedness - featureCounts - gene count matrix - look at the matrix |
| **C** | 5 figures (read quality, alignment and counting, liver vs muscle, top genes, counts along chromosome 19) |

## Run it in Google Colab (nothing to install)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pratikpathade/GenartWorkshop_Pratik/blob/main/Day2_RNAseq_FASTQ_to_counts_shell_tutorial.ipynb)

1. Click the badge, then **Connect** (top right).
2. Run the cells from top to bottom (play button, or *Runtime > Run all*). Setup (installing the tools) takes about 2-4 minutes, the analysis about 2-3 minutes. Tested end to end in 2.5 minutes on 2 CPU cores; a free Colab session can be slower.

The notebook `Day2_RNAseq_FASTQ_to_counts_shell_tutorial_with_outputs.ipynb` is the same notebook with all outputs saved from a test run (what you should see) - handy for reading without running or for comparing your own output.

## Run it as a script (Linux / Colab terminal)

```bash
git clone https://github.com/PratikPathade/GenartWorkshop_Pratik.git
cd GenartWorkshop_Pratik
bash Day2_RNAseq_FASTQ_to_counts_shell_tutorial.sh
```

Everything is made in `./RNAseq_calf`; the tools are installed to `~/genart-rnaseq` (micromamba: FastQC 0.13.0, MultiQC 1.35, fastp 1.4.0, HISAT2 2.2.3, samtools 1.24, Subread/featureCounts 2.1.1).

## Data (`data/`, about 37 MB)

| Files | What |
|---|---|
| `data/reads/liver_R1, liver_R2, muscle_R1, muscle_R2` (FASTQ) | Real reads of one Holstein-Friesian calf (138 days): liver and skeletal muscle, 100,000 read pairs each (2 x 101 bp, NovaSeq, reverse-stranded), only pairs from **chromosome 19** |
| `data/ref/chr19.fa.gz`, `chr19.gtf.gz` | Cattle chromosome 19 (ARS-UCD2.0) and its gene annotation, both from Ensembl release 114 |
| `data/make_subset.sh` | How the small reads were cut from the original 2.3 GB runs (instructor use) |

Sources and what was changed: `data/SOURCES.md`.

**Please cite** GEO series GSE289870 / ENA project PRJNA1225222 (University College Dublin) for the reads and Ensembl release 114 for the reference.
