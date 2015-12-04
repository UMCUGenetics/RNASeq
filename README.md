RNA-SEQ PIPELINE


0. Tasks
========

This pipeline performs the following tasks:
  - perform quality control on FastQ files (using FastQC)
  - align reads of each sample in a run against reference genome (using STAR) and add read groups (using Picard)
  - perform quality control on generated BAM files (using Picard)
  - count reads in features (using HTSeq-count)
  - normalize read counts (using DESeq)
  - calculate RPKMs (using edgeR)
  - perform DE analysis for standard designs (using DESeq2)


1. Running the pipeline
=======================

Run by typing: 
perl /hpc/cog_bioinf/common_scripts/RNA_seq_analysis/scripts/RNAseqAnalyse.pl -input [run directory] - output [output directory]

CAUTION: Paths to input and output directories must be absolute paths!

To see additional parameters, type:
perl /hpc/cog_bioinf/common_scripts/RNA_seq_analysis/scripts/RNAseqAnalyse.pl -h


2. Mapping
==========

RNA-seq reads are aligned to the reference genome using STAR, which was designed to specifically address many of the challenges of RNA-seq data mapping, and uses a novel strategy for spliced alignments.

3. Detecting fusion genes
=========================

The pipeline searches for fusion genes / chimeric junctions by default. The default minimum segment length is 1 and can be edited with option -chimSegmentMin.


4. Read counting
================

Counting sequencing reads in features (genes/transcripts/exons) is done with htseq-count using the union mode.
The raw read counts table can be found in the directory: <rundir>/read_counts/<run>_raw_counts.txt. This table contains the Ensembl (gene/transcript/exon) IDs (rows) and the raw read 
counts per library (columns).


5. Normalizing read counts
==========================

The raw read counts are normalized using the DESeq method included in the DESeq Bioconductor package and is based on the hypothesis that most genes are not DE. A DESeq scaling factor for
a given lane is computed as the median of the ratio, for each gene, of its read count over its geometric mean across all lanes. The underlying idea is that non-DE
genes should have similar read counts across samples, leading to a ratio of 1. Assuming most genes are not DE, the median of this ratio for the lane provides an
estimate of the correction factor that should be applied to all read counts of this lane to fulfill the hypothesis. By calling the estimateSizeFactors() and sizeFactors()
functions in the DESeq Bioconductor package, this factor is computed for each lane, and raw read counts are divided by the factor associated with their sequencing lane. 


6. Calculate RPKMs
==================

RPKM is a method of quantifying gene expression from RNA sequencing data by normalizing for total read length and the number of sequencing reads. RPKMs are calculated using the RPKM()
function included in the edgeR Bioconductor package.
RPKM = [# of mapped reads]/([length of transcript]/1000)/([total reads]/10^6)

CAUTION: Make sure the RPKM normalization is actually necessary in your analysis. If you are comparing gene expression among samples only, there really is no reason to normalize by length
as you will be dividing each gene among the samples by a constant (gene length). You only need to use RPKM when you are comparing transcript expression within one sample.


7. Differential expression analysis
===================================

Differential expression analysis can be done only for standard designs, such as:
   sample1	test
   sample2	control
   sample3	test
   sample4	test
   sample5	control

DE analysis is done using the DESeq2 Bioconductor package. It takes the merged raw read counts (from HTseq-count) as an input and creates the following output inside the /<run>/DEanalysis folder:
 - <run>_DEanalysis_all.txt : result table of DE analysis (ordered by p-value). The columns are:
   - Ensembl (gene/transcript/exon) ID
   - baseMean -> the average of the normalized count values, dividing by size factors, taken over all samples
   - log2FoldChange -> effect size estimate (tells us how much the gene's expression seems to have changed due to for example treatment in comparison to untreated samples
   - lfcSE -> standard error estimate for the log2 fold change estimate
   - stat
   - pvalue -> indicates the probability that a fold change as strong as the observed one, or oven stronger, would be seen under the situation described by the null hypothesis.
   - padj -> adjusted p-value using the Benjamini-Hochberg adjustment
   - gene_id
   - gene_name
 - <run>_MAplot.png: shows the log2 fold changes attributable to a given variable over the mean of normalized counts. Points are colored red if the adjusted p value is less than 0.1,
   points that fall out of the window are plotted as open triangles pointing either up or down.
 - <run>_sampletosample_distances.png: heatmap of the sample-to-sample distances.
 - <run>_PCAplot.png: principal component plot of the samples.


8. Additional tools
===================

Please contact Annelies Barendregt (F.A.S.Smouter@umcutrecht.nl) if you want to add additional tools/scripts/options.
