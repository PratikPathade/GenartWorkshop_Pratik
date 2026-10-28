# Sources of the data in this folder

## Reads (`reads/`)

Liver RNA-seq of three Danish Landrace boars, one from each boar-taint group. They are three of the nine boars used in the Day 2 R tutorials (same animal numbers).

| File prefix | Group | Animal | GEO sample | ENA run | Reads in the original run |
|---|---|---|---|---|---|
| `high_7107` | High | 7107 | GSM2461212 | SRR7041911 | 18,363,245 pairs |
| `medium_2111` | Medium | 2111 | GSM2461225 | SRR7041924 | 19,949,217 pairs |
| `low_1105` | Low | 1105 | GSM2461240 | SRR7041939 | 15,034,351 pairs |

* Series: GEO **GSE113316** = ENA / BioProject **PRJNA451072** (SRA study SRP141481), "Differential expression and co-expression gene networks reveal candidate biomarkers of boar taint in non-castrated pigs (RNA-seq data set)".
* Total RNA from liver and testis of 48 Danish Landrace pigs with low, medium and high genetic merit of boar taint; Illumina HiSeq 2500, paired-end 2 x 100 bases.
* Contributors: Markus Drag, Ruta Skinkyte-Juskiene, Duy N. Do, Lisette J. Kogelman, Haja N. Kadarmideen.

**What was changed:** nothing in the reads themselves. Only the first 8,000,000 read pairs of each run were used; they were aligned to pig chromosome 12 and the first 100,000 pairs that align properly to it were kept, unchanged (same names, bases and qualities). `make_subset.sh` does exactly this.
Because the reads were chosen by their alignment to chromosome 12, the alignment rate in the tutorial is much higher than in a real experiment.

## Reference (`ref/`)

* `chr12.fa.gz`: *Sus scrofa* chromosome 12, assembly **Sscrofa10.2**, Ensembl release 79 (`Sus_scrofa.Sscrofa10.2.dna.chromosome.12.fa.gz`), unchanged.
* `chr12.gtf.gz`: the lines of chromosome 12 from the Ensembl release 79 gene annotation (`Sus_scrofa.Sscrofa10.2.79.gtf.gz`).
* Release 79 is the version that matches the gene IDs of the count tables in the Day 2 R tutorials.

## Please cite

* Drag M, Skinkyte-Juskiene R, Do DN, et al. (2017). Differential expression and co-expression gene networks reveal candidate biomarkers of boar taint in non-castrated pigs. *Scientific Reports*. doi:10.1038/s41598-017-11928-0 (PubMed 28939879)
* Drag M, Hansen MB, Kadarmideen HN, et al. (2018). Systems genomics study reveals expression quantitative trait loci, regulator genes and pathways associated with boar taint in pigs. *PLoS ONE*. doi:10.1371/journal.pone.0192673 (PubMed 29438444)
* GEO series GSE113316 / ENA project PRJNA451072 for the reads, and Ensembl release 79 for the reference.
