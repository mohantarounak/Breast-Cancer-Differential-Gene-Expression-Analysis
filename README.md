# Breast Cancer DGE Analysis
## Dataset: GSE42568
# Abstract
This study presents a comprehensive Differential Gene Expression (DGE) 
analysis of breast cancer versus normal breast tissue using publicly 
available microarray data (GSE42568) from the NCBI Gene Expression 
Omnibus (GEO) database. A total of 121 samples (17 normal and 104 tumor) 
were analyzed using the limma package in R. We identified 3,277 
significantly differentially expressed genes (adj.P.Val < 0.05, 
|logFC| > 1.5). Gene Ontology (GO) and KEGG pathway enrichment analyses 
revealed significant involvement of biological processes including blood 
circulation, extracellular matrix organization, and ERK1/ERK2 cascade, 
providing insights into molecular mechanisms underlying breast cancer 
progression.

# Introduction
Breast cancer is one of the most common malignancies worldwide and a 
leading cause of cancer-related mortality among women. Understanding 
the molecular mechanisms driving breast cancer development and 
progression is critical for identifying novel biomarkers and therapeutic 
targets. Gene expression profiling using microarray technology provides 
a powerful approach to study genome-wide transcriptional changes between 
cancerous and normal tissues.

Differential Gene Expression (DGE) analysis enables the identification 
of genes that are significantly upregulated or downregulated in disease 
conditions compared to normal controls. In this study, we performed DGE 
analysis using the GSE42568 microarray dataset, which contains expression 
profiles of 104 breast cancer samples and 17 normal breast tissue samples. 
Our goal was to identify key differentially expressed genes and biological 
pathways associated with breast cancer.

## Methods
#1. Dataset
   - Database: NCBI Gene Expression Omnibus (GEO)
   - Dataset ID: GSE42568
   - Platform: Affymetrix Human Genome U133 Plus 2.0 Array
   - Samples: 17 Normal + 104 Tumor breast tissue samples
## Tools & Packages : 
   - R version 4.6.1
   - GEOquery: Data retrieval
   - limma: DGE Analysis
   - pheatmap: Heatmap generation
   - EnhancedVolcano: Volcano plot
   - clusterProfiler: GO & KEGG enrichment
   - hgu133plus2.db: Probe ID annotation
   - org.Hs.eg.db: Gene ID conversion

## Analysis: Differential Gene Expression
 # Analysis Pipeline
   ## Step 1: Data Loading & Preprocessing
           - Downloaded GSE42568 series matrix file
           - Extracted expression matrix and sample metadata

   ## Step 2: DGE Analysis (limma)
           - Built design matrix (Normal vs Tumor)
           - Applied linear model fitting (lmFit)
           - Computed empirical Bayes statistics (eBayes)
           - Extracted DEGs using topTable (FDR adjusted p-value)

   ## Step 3: Filtering Significant DEGs
           - Criteria: adj.P.Val < 0.05 and |logFC| > 1.5
           - Total significant DEGs identified: 3,277

   ## Step 4: Visualization
           - Volcano plot: LogFC vs adjusted p-value
           - Heatmap: Top 50 DEGs expression pattern

   ## Step 5: Functional Enrichment Analysis
           - GO Biological Process enrichment
           - KEGG Pathway enrichment analysis
## Comparison: Tumor vs Normal Breast Tissue
