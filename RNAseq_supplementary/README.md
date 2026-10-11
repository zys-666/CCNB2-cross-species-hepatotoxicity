# RNA-seq Supplementary Data — Reviewer Response (BZ20250424728812)

RNA-seq results for the NTH (2-aminobenzothiazole) exposure experiment in crucian carp
(***Carassius auratus***), project BZ20250424728812, samples C-1 … L-3. Provided in response
to three reviewer comments:

- **①** Methodological detail: sequencing quality metrics, mapping rate, reference genome information
- **②** Library numbers per group, RNA pooling status, sequencing depth, mapping rate,
  reference genome/annotation version, gene-ID conversion strategy, complete DEG/enrichment results
- **③** Quality-control figures: library size, PCA/MDS, sample correlation, dispersion

---

## Species and reference genome — please note

| | |
|---|---|
| Sequenced organism (exposed to NTH) | ***Carassius auratus*** — crucian carp / goldfish |
| Reference genome used for alignment | ***Carassius gibelio*** — gibel carp, strain **Cgi1373**, assembly **carGib1.2-hapl.c**, NCBI **GCF_023724105.1** |

Reference genome details as delivered: 1,583,351,535 bp, 51 sequences, GC 37.60%,
chromosome-level assembly; ecotype "wild population from Czech Republic", female,
submitted by IGB Leibniz-Institute of Freshwater Ecology and Inland Fisheries.

*C. auratus* and *C. gibelio* are closely related congeneric species, so the large majority of
clean reads still align and gene models transfer largely one-to-one. However, **all gene IDs
(`LOC1279…`), symbols and functional annotations in this dataset refer to *C. gibelio* gene
models**, not to a *C. auratus* assembly. This is stated explicitly in the manuscript.

---

## Differential expression

Raw gene counts were generated with **featureCounts v2.0.6**, and differential expression was
tested with **edgeR (R 3.6.3)** directly on the counts using TMM normalization. Genes were
declared differentially expressed at **FDR < 0.05 and |log2FC| > 1**, with the log2 fold change
calculated from FPKM values.

| Comparison | DEGs (total / up / down) |
|---|---|
| H_vs_C (5 mg/L vs control) | 7,234 / 2,041 / 5,193 |
| L_vs_C (0.5 mg/L vs control) | 8,141 / 3,162 / 4,979 |
| L_vs_H (0.5 vs 5 mg/L) | 9,298 / 6,243 / 3,055 |

These are the DEG sets reported in the manuscript, and all downstream analyses
(PPI / MCODE, GO and KEGG enrichment, GSEA, Venn) were performed on them.

---

## Directory guide

### `00_metadata/`
| File | Contents |
|---|---|
| `table.json.js` | Authoritative project metadata: reference genome identity, software versions and parameters, sample grouping, DEG/enrichment summary counts |

### `01_QC/` — raw & clean read quality (reviewer ①)
| File | Contents |
|---|---|
| `rawdata.stat.main.xls` | Per-sample raw reads, total bases, N%, A/T/C/G%, error rate, **Q20%, Q30%**, GC% |
| `cleandata.stat.main.xls` | Same metrics after Trimmomatic QC |
| `per_base_quality.png` | Per-base sequence quality (all samples) |
| `per_base_sequence_content.png` | Per-base sequence content (all samples) |

### `02_Align/` — mapping statistics (reviewer ①②)
| File | Contents |
|---|---|
| `align_sum.txt` | Per-sample clean reads, mapped reads, **mapping rate (%)** |

### `03_RNAQuality/` — post-alignment assessment (reviewer ①)
| File | Contents |
|---|---|
| `rRNA_sum.txt` | Per-sample rRNA reads and **rRNA residual rate (%)** |
| `region.stat.txt` | Read distribution over genomic regions (exon / intron / TSS / TES / intergenic) |

### `04_ExpressionQC/` — expression-level quality control (reviewer ③)
| File | Contents |
|---|---|
| `box.png` | **Library size** — expression distribution boxplot |
| `pca.png` | **PCA** — principal component analysis (group-coloured) |
| `correlation.heatmap.png` | **Sample correlation** heatmap |
| `exp.density.png` | **Dispersion** — expression density curves |
| `genes.tre.png` | Sample clustering tree |
| `fpkm_exp_stat.txt` | Number of detected genes at FPKM > 0, 1, 10, 100 |
| `genes.fpkm.correlation.xls` | Full Pearson correlation matrix |
| `pca_Group_sites.fpkm.xls` | PCA sample coordinates |
| `pca_Group_importance.fpkm.xls` | PCA variance explained per component |

### `05_Expression/` — expression matrix & annotation maps (reviewer ②)
| File | Contents |
|---|---|
| `genes.fpkm.txt` | Genome-wide FPKM matrix |
| `genes.read.txt` | Genome-wide **raw count** matrix (featureCounts output) |
| `ref_genome.func.xls` | **GeneID → protein product** mapping |
| `ref_genome.GO.list` | **GeneID → GO term** mapping |
| `ref_genome.pathway.txt` | **GeneID → KEGG KO / pathway** mapping |

Gene-ID conversion strategy: quantification and all downstream analysis were carried out on the
original reference gene IDs. IDs were **not** replaced by gene symbols; unnamed genes retain their
NCBI `LOC` numbers and named genes keep their symbols. Gene functions were taken from the reference
annotation, supplemented by best-hit NR protein annotation where the reference lacked a description.

### `06_DEG/` — differential expression (reviewer ②)
| File | Contents |
|---|---|
| `DEG.stat.xls` | DEG counts per comparison (total / up / down) |
| `deg.fpkm.bar.png` | DEG count barplot (FPKM-based, as used in the manuscript) |
| `deg.fpkm.volcano.png` | Volcano plot (FPKM-based, as used in the manuscript) |
| `details/*.fpkm.DEG.xls` | **Full DEG tables** for the three comparisons (logFC, P value, FDR, product annotation) |

### `07_Enrichment/` — functional enrichment (reviewer ②)
| File | Contents |
|---|---|
| `GO/*.xls` | GO enrichment — all DEGs / up-regulated / down-regulated |
| `GO_level2/*.level2.txt` | GO level-2 category statistics |
| `KEGG/*.xls` | KEGG pathway enrichment — all / up / down |
| `KEGG_pathway/*.pathway.table.xls` | KEGG pathway annotation tables |

### `08_GSEA/`
Six GSEA archives (GO and KEGG for each of H_vs_C, L_vs_C, L_vs_H), each containing the ranked
gene list, gene sets (.gmt), enrichment plots, heat maps and result tables.

### `09_Venn/`
DEG Venn diagram and per-comparison overlap detail tables.

### `99_methods/`
| File | Contents |
|---|---|
| `真核有参转录英文分析方法描述.pdf` | Provider's generic English methods **template**. Its version numbers are boilerplate defaults and do **not** match this project: the template states Trimmomatic v0.36 and HTSeq, whereas the project-specific software report records Trimmomatic **v0.39** and **featureCounts v2.0.6**. Where the two disagree, the project-specific table in `00_metadata/table.json.js` is authoritative. |
| `GSEA分析及结果说明.pdf` | GSEA analysis and result description |

---

## Key parameters

Taken from the sequencing provider's **project-specific** analysis report (`00_metadata/table.json.js`):

| Step | Software | Version | Parameters |
|---|---|---|---|
| Read QC | Trimmomatic | 0.39 | `ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75` |
| Alignment | Hisat2 | 2.2.1 | `-p 16 --dta-cufflinks` |
| rRNA assessment | Bowtie2 | — | clean reads vs rRNA database |
| Post-alignment QC | Qualimap | 2.2.1 | default |
| Quantification | featureCounts (subread) | 2.0.6 | `-p -C -B -P -O -T 16 -Q 20` |
| Differential expression | edgeR | R 3.6.3 | \|logFC\| ≥ 1 & FDR ≤ 0.05 |
| GO enrichment | Goatools | 0.9.9 | default |
| KEGG enrichment | KOBAS | 3 | default |
| Alternative splicing | rMATS | 4.3.0 | default |
| SNP calling | BCFtools | 1.17 | `call -a FORMAT/GQ -c -v -V indels` |

Other characteristics:

- **Platform**: Illumina NovaSeq 6000, paired-end 150 bp × 2
- **Design**: 9 independent libraries, 3 biological replicates per group
  (C: C-1/C-2/C-3; H: H-1/H-2/H-3; L: L-1/L-2/L-3); **no RNA pooling**;
  comparisons H_vs_C, L_vs_C, L_vs_H
- **Sequencing depth**: 41.2–53.4 M clean reads per sample (6.12–7.93 Gb)
- **Mapping rate**: 73.20–82.08%
- **Enrichment thresholds**: GO — Goatools; KEGG — KOBAS; Bonferroni-corrected *P* < 0.05

## Not included in this repository
Raw genome sequence files and the full SNP annotation tables are excluded — they are large,
not required for the reviewer response, and exceed the per-file limit of this repository.
