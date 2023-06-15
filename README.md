## TFM. ChIP-seq Analysis workflow research: Mapping DNA double-strand breaks induced by Topoisomerase I during transcription reveals a nom-random nature of these lesions alongside the genome  

The present project includes the scripts employed for the analysis of Chromatin Immunoprecipitation sequencing (ChIP-seq) data to decipher the process involved in genome instability due to TOP1-induced DSBs during transcription.

## 1.- Peak calling step: pipeline nf-core/chipseq
ChiP-seq datasets were processed employing the nf-core/chipseq pipeline (Ewels et al., 2020).

script_XRCC4.sh

script_TOP1-seq.sh
    
## 2.- Statistical analysis: R scripts
### 2.1.- Script to detect overlapping peaks across replicates 
R_peak_intersection.R

### 2.2.- Script for peak annotation and overlapping analyses on annotated genes. Graphical representation of annotated peaks relative to gene features and Gene Ontology Enrichment Analysis   
R_TSS_geneFeatures.R

### 2.3.- Script for enhancers annotation of distal intergenic peaks
R_enhancer_annotation.R

## Pipeline and Packages references

Carlson M (2019). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.8.2 

Ewels, P. A., Peltzer, A., Fillinger, S., Patel, H., Alneberg, J., Wilm, A., Garcia, M. U., Di Tommaso, P., & Nahnsen, S. (2020). The nf-core framework for community-curated bioinformatics pipelines. Nature biotechnology, 38(3), 276–278. https://doi.org/10.1038/s41587-020-0439-x 

Wang, Q., Li, M., Wu, T., Zhan, L., Li, L., Chen, M., Xie, W., Xie, Z., Hu, E., Xu, S., & Yu, G. (2022). Exploring Epigenomic Datasets by ChIPseeker. Current protocols, 2(10), e585. https://doi.org/10.1002/cpz1.585 

Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. Innovation (Cambridge (Mass.)), 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141 

Yu, G., & He, Q. Y. (2016). ReactomePA: an R/Bioconductor package for reactome pathway analysis and visualization. Molecular bioSystems, 12(2), 477–479. https://doi.org/10.1039/c5mb00663e 

Zhu, L. J., Gazin, C., Lawson, N. D., Pagès, H., Lin, S. M., Lapointe, D. S., & Green, M. R. (2010). ChIPpeakAnno: a Bioconductor package to annotate ChIP-seq and ChIP-chip data. BMC bioinformatics, 11, 237. https://doi.org/10.1186/1471-2105-11-237 
