# CAVD_COMP_PIPELINES
Computational pipelines for the analysis of RNASeq data in the context of the research project: "Identification of a peripheral blood gene signature
predicting aortic valve calcification".

Pipelines are distributed in four folders:
1. PREPROCESSING
2. ALIGN_AND_QUANTIFY
3. DEA_RNASeqEdgeR-4.0.4
4. CLUSTERING_AND_ANNOTATION

## PREPROCESSING
This folder contains a single shell script that calls cutadapt, to eliminate Illumina adapter remains, and FastQC,
to perform some quality tests on the reads, after preprocessing.

## ALIGN_AND_QUANTIFY
This folder contains a single shell script that calls RSEM to align preprocessed reads against a transcriptome reference and quantify gene expression levels.

## DEA_RNASeqEdgeR-4.0.4
This folder contains four R scripts that call EdgeR or limma, to normalize gene expression counts estimated by the RSEM pipeline,
and to perform differential expression analysis.

The output consists in the following collection of files:
* Two files named condA_vs_condB.txt and condA_vs_condB.xls for each of the contrasts being performed, describing raw counts and normalized expression values for all genes
in the samples representing the two conditions being compared, as well as fold change and log fold change values, and raw p_value and Benjamini-Hochberg (BH) adjusted p_values.
Genes are considered to be differentially expressed if the change in expression between the two conditions is associated to BH adjusted p_values < 0.05.
* A file named DEG_AllContrasts.xlsx with as many sheets as contrasts being performed, each of them containing the subset of analysed genes that were classified as
differentially expressed (with BH adjusted p_value < 0.05).
* A file with normalized counts for all genes with detectable expression, in all samples considered for analysis.
* Several types of diagnostic plots, contained within a folder called "Diagnostics".

## CLUSTERING_AND_ANNOTATION
This folder contains three R scripts and one Perl script that identify groups of genes with similar expression profiles, using the kmeans method.
In the context of the current research project, gene clusters were manually subject to functional enrichment analyses with the Canonical Pathway
analysis component of IPA (Qiagen).
After reformating IPA output files, the last script filters enrichment results and performs a concordance analysis to identify groups of clusters
enriched in similar collections of pathways.
Such groups of clusters were used to define meta-clusters, that were re-anotated with IPA.

