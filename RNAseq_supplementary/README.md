# RNA-seq Supplementary Data — Reviewer Response (BZ20250424728812)

RNA-seq results for the crucian carp (**Carassius gibelio**) NTH experiment
(project BZ20250424728812, samples C-1…L-3), provided in response to three reviewer comments:

- **①** Methodological detail: sequencing quality metrics, mapping rate, reference genome information
- **②** Library numbers per group, RNA pooling status, sequencing depth, mapping rate,
  reference genome/annotation version, gene-ID conversion strategy, complete DEG/enrichment results
- **③** Quality-control figures: library size, PCA/MDS, sample correlation, dispersion

---

## Directory guide

### `00_metadata/`
| File | Contents |
|---|---|
| `table.json.js` | Authoritative project metadata: reference genome identity, all software versions & parameters, sample grouping, DEG/enrichment summary counts |

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
| `tpm_exp_stat.txt`, `fpkm_exp_stat.txt` | Number of detected genes at TPM/FPKM > 0, 1, 10, 100 |
| `genes.tpm.correlation.xls`, `genes.fpkm.correlation.xls` | Full Pearson correlation matrices |
| `pca_Group_pc1-2.tpm.pdf` | PCA scatter plot (TPM) |
| `pca_Group_sites.tpm.xls`, `pca_Group_sites.fpkm.xls` | PCA sample coordinates |
| `pca_Group_importance.tpm.xls`, `pca_Group_importance.fpkm.xls` | PCA variance explained per component |

### `05_Expression/` — expression matrices & annotation maps (reviewer ②)
| File | Contents |
|---|---|
| `genes.tpm.txt` | Genome-wide TPM matrix (43,903 genes × 9 samples) |
| `genes.fpkm.txt` | Genome-wide FPKM matrix |
| `genes.read.txt` | Genome-wide raw HTSeq count matrix |
| `ref_genome.func.xls` | **GeneID → protein product** mapping |
| `ref_genome.GO.list` | **GeneID → GO term** mapping |
| `ref_genome.pathway.txt` | **GeneID → KEGG KO / pathway** mapping |

### `06_DEG/` — differential expression (reviewer ②)
| File | Contents |
|---|---|
| `DEG.stat.xls` | DEG counts per comparison (TPM & FPKM; total / up / down) |
| `deg.tpm.bar.png` | DEG count barplot |
| `deg.tpm.volcano.png`, `deg.fpkm.volcano.png` | Volcano plots |
| `details/*.tpm.DEG.xls`, `details/*.fpkm.DEG.xls` | **Full DEG tables** (logFC, P value, FDR, product annotation) |
| `all_expression/*.all.exp.xls` | Statistics for all genes, not only DEGs |

### `07_Enrichment/` — functional enrichment (reviewer ②)
| File | Contents |
|---|---|
| `GO/*.xls` | GO enrichment — all DEGs / up-regulated / down-regulated |
| `GO_level2/*.level2.txt` | GO level-2 category statistics |
| `KEGG/*.xls` | KEGG pathway enrichment — all / up / down |
| `KEGG_pathway/*.pathway.table.xls` | KEGG pathway annotation tables |

### `08_GSEA/`
Six GSEA archives (GO and KEGG for each of H_vs_C, L_vs_C, L_vs_H), each containing
the ranked gene list, gene sets (.gmt), enrichment plots, heat maps and result tables.

### `09_Venn/`
DEG Venn diagram and per-comparison overlap detail tables.

### `99_methods/`
| File | Contents |
|---|---|
| `真核有参转录英文分析方法描述.pdf` | Vendor (Shanghai BIOZERON) English methods description |
| `GSEA分析及结果说明.pdf` | GSEA analysis and result description |

---

## Key parameters (summary)

- **Platform**: Illumina NovaSeq 6000, paired-end 150 bp × 2
- **Reference genome**: *Carassius gibelio* strain Cgi1373, assembly **carGib1.2-hapl.c**,
  NCBI **GCF_023724105.1**; 1,583,351,535 bp; 51 sequences; GC 37.60%; chromosome-level
- **QC**: Trimmomatic 0.39 (`ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75`)
- **Alignment**: Hisat2 2.2.1 (`-p 16 --dta-c`)
- **Post-alignment QC**: Qualimap v2.2.1
- **Quantification**: HTSeq (gene-level read counts) → FPKM / TPM
- **Differential expression**: edgeR, threshold **|logFC| > 1 and FDR < 0.05**
- **Enrichment**: GO — Goatools; KEGG — KOBAS; both with Bonferroni-corrected *P* < 0.05
- **Splicing**: rMATS (JCEC)
- **Design**: 9 independent libraries, 3 biological replicates per group
  (C: C-1/C-2/C-3; H: H-1/H-2/H-3; L: L-1/L-2/L-3); **no RNA pooling**;
  comparisons H_vs_C, L_vs_C, L_vs_H

## Not included in this repository
Raw genome sequence files and the full SNP annotation tables are excluded — they are large,
not required for the reviewer response, and exceed the per-file limit of this repository.
