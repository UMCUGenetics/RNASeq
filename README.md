## RNA-SEQ PIPELINE
This pipeline performs the following tasks:
- perform quality control on FastQ files (using FastQC)
- align reads of each sample in a run against reference genome (using STAR) and add read groups (using Picard)
- perform quality control on generated BAM files (using Picard)
- count reads in features (using HTSeq-count)
- normalize read counts (using DESeq)
- calculate RPKMs (using edgeR)
- perform DE analysis for standard designs (using DESeq2)
- variant calling, filtering and annotation

## Download
Use git clone:
```bash
git clone git@github.com:UMCUGenetics/RNASeq.git
```

## Installation
[Download](#download) the RNAseq pipeline.
Make sure all [dependencies](#dependencies) are installed and the right paths are set in the pipeline (RNAseqAnalyse.pl) in the "Get options" section.

#### Genome files
Generate genome indexes files using the instructions in section [Generate genome indexes](#generate-genome-indexes). The genome indexes are saved to disk and need only be generated once for each genome/annotation combination.
Next to the files you had to collect to generate the genome indexes, you need:
- refFlat file (for using bamMetrics)
- Interval list (for using bamMetrics)
- Genesizes file (for calculating RPKMs)

## Usage
#### Run pipeline
```bash
perl RNAseqAnalyse.pl -input [/path/to/rundir] -outputDir [/path/to/outputdirname] -mail [email]
```
To see additional parameters, just type:
```bash
perl RNAseqAnalyse.pl
```

## Dependencies
#### Core tools
- Opengrid engine
- Perl 5
- Python 2.7
- R 3.2.2
- Java 1.7

#### Bio tools
- [STAR 2.4.2a](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.2a)
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Sambamba](http://lomereiter.github.io/sambamba/)
- [bamMetrics](https://github.com/CuppenResearch/bamMetrics)
- [GATK (genomeanalysis toolkit) >= 3.2-2](https://www.broadinstitute.org/gatk/)
- [Picard >= 1.119](http://broadinstitute.github.io/picard/) 
- [IGVtools](https://www.broadinstitute.org/igv/igvtools)
- [SnpEff / SnpSift](http://snpeff.sourceforge.net/)

#### Perl modules
- strict
- warnings
- POSIX
- Getopt::Long
- File::Basename
- File::Path
- Cwd
- List::MoreUtils
- Time::localtime

#### Python packages
- [HTSeq 0.6.1](https://pypi.python.org/pypi/HTSeq)

#### R packages
- DESeq 1.18.0
- DESeq2 1.6.3
- edgeR 3.8.6
- ggplot2
- gplots
- RColorBrewer
- Bioconductor annotations for:
    - Human: org.Hs.eg.db
    - Rat: org.Rn.eg.db
    - Mouse: org.Mm.eg.db
    - Zebrafish: org.Dr.eg.db
    - Dog: org.Cf.eg.db
    - Arabidopsis: org.At.tair.db

#### Databases
- dbNSFP 2.9

## Generate genome indexes
Create a directory where you want to store the indexes (e.g. /GENOMES/STAR/Homo_sapiens.GRCh37).

Collect the following files for your genome:
- Fasta file containing the genome reference sequences
- GTF file containing annotated transcripts

Run STAR:
```bash
STAR \
    --runMode genomeGenerate \
    --genomeDir /path/to/genomeDir \
    --genomeFastaFiles /path/to/genome/fasta.fa \
    --runThreadN 4 \
    --sjdbGTFfile /path/to/annotations.gtf
```

## Output description

#### Mapping
RNA-seq reads are aligned to the reference genome using STAR, which was designed to specifically address many of the challenges of RNA-seq data mapping, and uses a novel strategy for spliced alignments.

#### Detecting fusion genes
The pipeline searches for fusion genes / chimeric junctions by default. The default minimum segment length is 1 and can be edited with option -chimSegmentMin.

#### Read counting
Counting sequencing reads in features (genes/transcripts/exons) is done with htseq-count using the union mode.
The raw read counts table can be found in the directory: <rundir>/read_counts/<run>_raw_counts.txt. This table contains the Ensembl (gene/transcript/exon) IDs (rows) and the raw read 
counts per library (columns).

#### Normalizing read counts
The raw read counts are normalized using the DESeq method included in the DESeq Bioconductor package and is based on the hypothesis that most genes are not DE. A DESeq scaling factor for
a given lane is computed as the median of the ratio, for each gene, of its read count over its geometric mean across all lanes. The underlying idea is that non-DE
genes should have similar read counts across samples, leading to a ratio of 1. Assuming most genes are not DE, the median of this ratio for the lane provides an
estimate of the correction factor that should be applied to all read counts of this lane to fulfill the hypothesis. By calling the estimateSizeFactors() and sizeFactors()
functions in the DESeq Bioconductor package, this factor is computed for each lane, and raw read counts are divided by the factor associated with their sequencing lane. 

#### Calculate RPKMs
RPKM is a method of quantifying gene expression from RNA sequencing data by normalizing for total read length and the number of sequencing reads. RPKMs are calculated using the RPKM()
function included in the edgeR Bioconductor package.
RPKM = [# of mapped reads]/([length of transcript]/1000)/([total reads]/10^6)

CAUTION: Make sure the RPKM normalization is actually necessary in your analysis. If you are comparing gene expression among samples only, there really is no reason to normalize by length
as you will be dividing each gene among the samples by a constant (gene length). You only need to use RPKM when you are comparing transcript expression within one sample.

#### Differential expression analysis
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
- <run>_MAplot.png: shows the log2 fold changes attributable to a given variable over the mean of normalized counts. Points are colored red if the adjusted p value is less than 0.1, points that fall out of the window are plotted as open triangles pointing either up or down.
- <run>_sampletosample_distances.png: heatmap of the sample-to-sample distances.
- <run>_PCAplot.png: principal component plot of the samples.

#### Variant calling, filtering and annotation
For the variant calling and filtering, GATK's HaplotypeCaller and VariantFiltration tool is used.
SnpEff is used to add genetic variant annotation and effect predictions to the vcf.
In case of human genome data, also the vcf will be annotated with the dbNSFP database using SnpSift.

#### Additional tools
Please contact Sander Boymans (S.W.Boymans@umcutrecht.nl) if you want to add additional tools/scripts/options.
