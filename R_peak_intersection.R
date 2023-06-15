########################################################
# Script to detect overlapping peaks across replicates #
########################################################

# Packages installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("ChIPpeakAnno")
BiocManager::install("ChIPseeker")
install.packages("ggplot2")



# Load libraries
library(ChIPpeakAnno)
library(ChIPseeker)
library(VennDiagram)
library(ggplot2)


# LOAD DATA: source nf-core pipeline (MACS2 PEAKS ANNOTATION FILES)

## Consensus peaks across replicates: Absolute Path bed files
XRCC4_CPT.peaks <- readPeakFile(
  peakfile = "/mnt/beegfs/home/marsabbon/chip-seq/results_FDR_xrcc4_150/bwa/mergedLibrary/macs2/broadPeak/consensus/XRCC4/XRCC4.consensus_peaks.bed",
  header=FALSE)
TOP1cc.peaks <- readPeakFile(
  peakfile = "/mnt/beegfs/home/marsabbon/chip-seq/results_baranello_broad_peak/bwa/mergedLibrary/macs2/broadPeak/consensus/TOP1/TOP1.consensus_peaks.bed",
  header = FALSE)

## Replicates' loading
XRCC4_CPT_REP1.peaks <- readPeakFile(
  peakfile = "/mnt/beegfs/home/marsabbon/chip-seq/results_FDR_xrcc4_150/bwa/mergedLibrary/macs2/broadPeak/CPTplusV5_IP_REP1_peaks.broadPeak",
  header=FALSE)
XRCC4_CPT_REP2.peaks <- readPeakFile(
  peakfile = "/mnt/beegfs/home/marsabbon/chip-seq/results_FDR_xrcc4_150/bwa/mergedLibrary/macs2/broadPeak/CPTplusV5_IP_REP2_peaks.broadPeak",
  header=FALSE)
XRCC4_CPT_REP3.peaks <- readPeakFile(
  peakfile = "/mnt/beegfs/home/marsabbon/chip-seq/results_FDR_xrcc4_150/bwa/mergedLibrary/macs2/broadPeak/CPTplusV5_IP_REP3_peaks.broadPeak",
  header=FALSE)


TOP1cc_REP1.peaks <- readPeakFile(
  peakfile = "/mnt/beegfs/home/marsabbon/chip-seq/results_baranello_broad_peak/bwa/mergedLibrary/macs2/broadPeak/TOPcc_REP1_peaks.broadPeak",
  header=FALSE)
TOP1cc_REP2.peaks <- readPeakFile(
  peakfile = "/mnt/beegfs/home/marsabbon/chip-seq/results_baranello_broad_peak/bwa/mergedLibrary/macs2/broadPeak/TOPcc_REP2_peaks.broadPeak",
  header=FALSE)

#################################################
# 1.  Peaks overlaps analysis across replicates #
#################################################

##library(ChIPpeakAnno)

## findOverlapsOfPeaks function for peak intersection analysis

# 1.1 XRCC4 ChIP-seq overlapping peaks
ol_XRCC4 <-findOverlapsOfPeaks(XRCC4_CPT_REP1.peaks, 
                          XRCC4_CPT_REP2.peaks, 
                          XRCC4_CPT_REP3.peaks, 
                          maxgap=1000) # > 1 KB gap between peak-intersection its discarded

# Save the list of overlapping peaks
overlappingPeaks_XRCC4 <- ol_XRCC4$overlappingPeaks

# Graphical analysis of the set of overlapping peaks : overlapping features
names(overlappingPeaks_XRCC4)
pie1(table(overlappingPeaks_XRCC4[["XRCC4_CPT_REP1.peaks///XRCC4_CPT_REP2.peaks"]]$overlapFeature))
pie1(table(overlappingPeaks_XRCC4[["XRCC4_CPT_REP1.peaks///XRCC4_CPT_REP3.peaks"]]$overlapFeature))
pie1(table(overlappingPeaks_XRCC4[["XRCC4_CPT_REP2.peaks///XRCC4_CPT_REP3.peaks"]]$overlapFeature))


# Calculation of totalTest parameter: number of possible binding sites calculated as
# human genome length GRCh38: 3099734149 divided average length of peaks
peaks.XRCC4_CPT <- data.frame(XRCC4_CPT.peaks)
XRCC4_binding_sites <- 3099734149 / mean(peaks.XRCC4_CPT$width)


# Venn Diagram representation and statistical analysis of intersection significance
result_XRCC4 <-makeVennDiagram(ol_XRCC4, totalTest=XRCC4_binding_sites ,
                          category = c('REP1','REP2','REP3'),
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1.5,
                          col=c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.6), alpha('#21908dff',0.6), alpha('#fde725ff',0.6)),
                          cex = 1,
                          fontfamily = "sans",
                          cat.cex = 1,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.085),
                          cat.fontfamily = "sans",
                          rotation = 1)

# p-value Hypergeometric test
result_XRCC4$p.value

# 1.2 TOP1-seq overlapping peaks
ol_TOP1 <- findOverlapsOfPeaks(TOP1cc_REP1.peaks,
                               TOP1cc_REP2.peaks,
                               maxgap=1000)

# Save the list of overlapping peaks
overlappingPeaks_TOP1 <- ol_TOP1$overlappingPeaks

# Graphical analysis of the set of overlapping peaks
names(overlappingPeaks_TOP1)
pie1(table(overlappingPeaks_TOP1[["TOP1cc_REP1.peaks///TOP1cc_REP2.peaks"]]$overlapFeature))

# Calculation of totalTest parameter: number posible binding sites
peaks.TOP1cc <- data.frame(TOP1cc.peaks)
TOP1_binding_sites <- 3099734149/ mean(peaks.TOP1cc$width)

# Venn Diagram representation and statistical analysis of intersection significance
result_TOP1 <-makeVennDiagram(ol_TOP1, totalTest=TOP1_binding_sites,
                          category = c('REP1','REP2'),
                          height = 480 , 
                          width = 480 , 
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1.5,
                          col=c("#440154ff", '#21908dff'),
                          fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                          cex = 1,
                          fontfamily = "sans",
                          cat.cex = 1,
                          cat.pos = c(-180, 180),
                          cat.default.pos = "outer",
                          alpha = c(0.8,0.9))

# p-value Hypergeometric test
result_TOP1$p.value


# 1.3 Peaks intersection of XRCC4 and active TOP1(TOP1cc)
ol_TOP1vsXRCC4 <- findOverlapsOfPeaks(TOP1cc.peaks,
                                      XRCC4_CPT.peaks,
                                      maxgap=1000)

# Save the list of overlapping peaks
overlappingPeaks_TOP1vsXRCC4 <- ol_TOP1vsXRCC4$overlappingPeaks

# Graphical analysis of the set of overlapping peaks
names(overlappingPeaks_TOP1vsXRCC4)
pie1(table(overlappingPeaks_TOP1vsXRCC4[["TOP1cc.peaks///XRCC4_CPT.peaks"]]$overlapFeature))


# Venn Diagram representation and statistical analysis of intersection significance
result_TOP1vsXRCC4 <-makeVennDiagram(ol_TOP1vsXRCC4, totalTest=TOP1_binding_sites,
                              category = c('TOP1','XRCC4'),
                              main = 'TOP1vsXRCC4',
                              height = 480 , 
                              width = 480 , 
                              resolution = 300,
                              compression = "lzw",
                              lwd = 1.5,
                              col=c("#440154ff", '#21908dff'),
                              fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
                              cex = 1,
                              fontfamily = "sans",
                              cat.cex = 1,
                              cat.pos = c(-180, 180),
                              cat.default.pos = "outer",
                              alpha = c(0.8,0.9))

# p-value Hypergeometric test
result_TOP1vsXRCC4$p.value


