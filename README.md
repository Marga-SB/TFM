# TFM: ChIP-seq Analysis workflow

The present project includes the scripts employed for Chip-seq analysis

## 1.- Peak calling step: pipeline nf-core/chipseq
script_XRCC5.sh
script_TOP1-seq.sh
    
## 2.- Statistical analysis: R scripts
### 2.1.- Script to detect overlapping peaks across replicates 
R_peak_intersection.R

### 2.2.- Script for peak annotation and overlapping analyses on annotated genes. Graphical representation of annotated peaks relative to gene features and Gene Ontology Enrichment Analysis   
R_TSS_geneFeatures.R

### 2.3 Enhancers annotation of distal intergenic peaks
R_enhancer_annotation.R

    
