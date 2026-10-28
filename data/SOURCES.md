# Data: sources and what was changed

## Reads (`data/reads/`) - real RNA-seq of one Holstein-Friesian calf

| Sample | Tissue | ENA run | GEO sample | Read pairs in the tutorial | File size (R1 + R2) |
|---|---|---|---|---|---|
| `liver`  | liver           | SRR32383427 | GSM8800637 (Liver_AI_Female_2)  | 100,000 | 4.9 + 5.0 MB |
| `muscle` | skeletal muscle | SRR32383411 | GSM8800653 (Muscle_AI_Female_2) | 100,000 | 3.6 + 3.6 MB |

**Source.** GEO series **GSE289870** = European Nucleotide Archive project **PRJNA1225222** (University College Dublin; dairy calves born after artificial insemination or in-vitro embryo production, liver, muscle and blood RNA-seq at about 4 months of age). The two samples are liver and muscle of the same female calf (AI group, 138 days old at biopsy). Library: polyA RNA, TruSeq stranded mRNA (**reverse-stranded**, `featureCounts -s 2`), Illumina NovaSeq 6000, paired-end 2 x 101 bp, about 23 million read pairs (2.3 GB) per run. Public data. The authors' own processing was FastQC + Trimmomatic + STAR (ARS-UCD1.3) + featureCounts.

**What was done to the reads** (script `data/make_subset.sh`, instructor use only):

1. Only the **first 4,000,000 read pairs** of each run were streamed from ENA (about 200 MB per mate instead of the full file).
2. These were aligned with HISAT2 to cattle **chromosome 19** (Ensembl release 114, ARS-UCD2.0). 206,696 liver pairs and 512,004 muscle pairs aligned properly.
3. The **first 100,000 of these pairs** (original file order) were written out: **both mates, exactly as the sequencer wrote them** (read names, bases, qualities). Nothing was trimmed, corrected or re-sampled.

**What this means for the data** (say it openly in class):

* Every read pair comes from chromosome 19 by construction, so the alignment rate in the tutorial (> 99 %) is much higher than in a real experiment (85-95 %), and the numbers of "genes detected" are those of one chromosome (1,771 genes), not of the genome.
* Because only the first part of each file was used, the reads come from the first tiles of one lane; they are not a random sample of the whole run, but the gene counts behave like a normal (small) RNA-seq sample.
* Muscle is dominated by a few very abundant genes (the myosin heavy chain MYH1 alone has more than half of the counts on chromosome 19) - a real property of skeletal muscle that is useful for talking about normalisation.
* One sample per tissue = no replicates: the matrix shows the biology but cannot be used for statistics.

## Reference (`data/ref/`)

| File | Source |
|---|---|
| `chr19.fa.gz` | Ensembl release 114, `Bos_taurus.ARS-UCD2.0.dna.primary_assembly.19.fa.gz` (unchanged), 63,449,741 bases |
| `chr19.gtf.gz` | Ensembl release 114, `Bos_taurus.ARS-UCD2.0.114.gtf.gz`, header lines and all lines of chromosome 19 (1,771 genes) |

Ensembl data are free to use (https://www.ensembl.org/info/about/legal/disclaimer.html).

## Please cite
* The GEO series GSE289870 / ENA project PRJNA1225222.
* Ensembl: Harrison PW et al. (2024) Ensembl 2024. Nucleic Acids Research 52:D891-D899.
* Tools: FastQC (Andrews), MultiQC (Ewels et al. 2016), fastp (Chen et al. 2018), HISAT2 (Kim et al. 2019), samtools (Danecek et al. 2021), featureCounts (Liao et al. 2014).
