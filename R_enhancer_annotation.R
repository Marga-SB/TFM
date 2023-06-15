###################################################
# Enhancers annotation of distal intergenic peaks #
###################################################

# Packages installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("GenomicRanges")

library(GenomicRanges)


## Consensus peaks across replicates variables: transformation into dataframe
peaks.XRCC4_CPT <- data.frame(XRCC4_CPT.peaks)
peaks.TOP1cc <- data.frame(TOP1cc.peaks)


# Download enhancers of RPE cells: http://www.enhanceratlas.org/data/download/enhancer/hs/RPE.bed
enhancers_RPE <- read.csv(file = 'enhancer_RPE.csv', sep = '\t',
                          header = T)


# Download enhancers of HCT116 cells : http://www.enhanceratlas.org/data/download/enhancer/hs/HCT116.bed
enhancers_HCT116 <- read.csv(file = 'enhancer_HCT116.csv', sep = '\t',
                             header = T)

# Build GRanges objects

enhancers_RPE.gr <- GRanges(seqnames = enhancers_RPE$sequence, 
                            ranges = IRanges(start = enhancers_RPE$start, 
                                             end = enhancers_RPE$end))

enhancers_HCT116.gr <- GRanges(seqnames = enhancers_HCT116$sequence, 
                               ranges = IRanges(start = enhancers_HCT116$start, 
                                                end = enhancers_HCT116$end))

XRCC4_CPT.gr <- GRanges(seqnames = peaks.XRCC4_CPT$seqnames, 
                        ranges = IRanges(start = peaks.XRCC4_CPT$start, end = peaks.XRCC4_CPT$end))

TOP1cc.gr <- GRanges(seqnames = peaks.TOP1cc$seqnames, 
                     ranges = IRanges(start = peaks.TOP1cc$start, end = peaks.TOP1cc$end))



## Overlaps of Distal intergenic peaks & Enhancers RPE cells
overlaps_enhancers_RPE <- intersect(enhancers_RPE.gr, XRCC4_CPT.gr)

# Determine number of enhancers annotated
length(overlaps_enhancers_RPE)

# Overlaps of Distal intergenic peaks & Enhancers HCT116 cells
overlaps_enhancers_HCT116 <- intersect(enhancers_HCT116.gr, TOP1cc.gr)

# Determine number of enhancers annotated
length(overlaps_enhancers_HCT116)


