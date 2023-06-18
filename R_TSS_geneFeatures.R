###########################################################################
# Script for peak annotation and overlapping analyses on annotated genes. #
# Graphical representation of annotated peaks relative to gene features   #
# and Gene Ontology Enrichment Analysis                                   #
###########################################################################

# Packages installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", force = T)
BiocManager::install("ChIPpeakAnno")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("ggplot2")
BiocManager::install("GenomicRanges")

# Load libraries for peaks annotation
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

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

#####################
# Peaks annotation  #
#####################

# Human Genome GRCh38 annotated genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# Annotation of XRCC4 ChIP-seq replicates

XRCC4_CPT_REP1.peakAnno <- annotatePeak(peak = XRCC4_CPT_REP1.peaks, 
                                        tssRegion=c(-1000, 1000),
                                        TxDb=txdb)
XRCC4_CPT_REP2.peakAnno <- annotatePeak(peak = XRCC4_CPT_REP2.peaks, 
                                        tssRegion=c(-1000, 1000),
                                        TxDb=txdb)
XRCC4_CPT_REP3.peakAnno <- annotatePeak(peak = XRCC4_CPT_REP3.peaks, 
                                        tssRegion=c(-1000, 1000),
                                        TxDb=txdb)

# Annotation of TOP1-seq replicates

TOP1cc_REP1.peakAnno <- annotatePeak(peak = TOP1cc_REP1.peaks, 
                                     tssRegion=c(-1000, 1000),
                                     TxDb=txdb)
TOP1cc_REP2.peakAnno <- annotatePeak(peak = TOP1cc_REP2.peaks, 
                                     tssRegion=c(-1000, 1000),
                                     TxDb=txdb)
# Consensus peaks annotation

XRCC4_CPT.peakAnno <- annotatePeak(peak = XRCC4_CPT.peaks, 
                                   tssRegion=c(-1000, 1000),
                                   TxDb=txdb)

TOP1cc.peakAnno <- annotatePeak(peak = TOP1cc.peaks, 
                                tssRegion=c(-1000, 1000),
                                TxDb=txdb)

# Transform annotated objects to dataframe

## Consensus peaks
XRCC4_CPT.annotation <- as.data.frame(XRCC4_CPT.peakAnno)
TOP1cc.annotation <- as.data.frame(TOP1cc.peakAnno)

## Peaks replicates XRCC4 ChIP-seq
XRCC4_CPT_REP1.annotation <- as.data.frame(XRCC4_CPT_REP1.peakAnno)
XRCC4_CPT_REP2.annotation <- as.data.frame(XRCC4_CPT_REP2.peakAnno)
XRCC4_CPT_REP3.annotation <- as.data.frame(XRCC4_CPT_REP3.peakAnno)

## Peaks replicates TOP1-seq
TOP1cc_REP1.annotation <- as.data.frame(TOP1cc_REP1.peakAnno)
TOP1cc_REP2.annotation <- as.data.frame(TOP1cc_REP2.peakAnno)

## Target genes selection
target.genes_XRCC4_CPT <- unique(XRCC4_CPT.annotation$geneId)
target.genes_TOP1cc <- unique(TOP1cc.annotation$geneId)


#GRanges construction for peaks intersection analysis

## Load library for GRanges objects construction
library(GenomicRanges)

XRCC4_REP1_genes_target.gr <- GRanges(seqnames = XRCC4_CPT_REP1.annotation$seqnames, 
                                       ranges = IRanges(start = XRCC4_CPT_REP1.annotation$start, 
                                                        end = XRCC4_CPT_REP1.annotation$end))
XRCC4_REP2_genes_target.gr <- GRanges(seqnames = XRCC4_CPT_REP2.annotation$seqnames, 
                                       ranges = IRanges(start = XRCC4_CPT_REP2.annotation$start, 
                                                        end = XRCC4_CPT_REP2.annotation$end))
XRCC4_REP3_genes_target.gr <- GRanges(seqnames = XRCC4_CPT_REP3.annotation$seqnames, 
                                       ranges = IRanges(start = XRCC4_CPT_REP3.annotation$start, 
                                                        end = XRCC4_CPT_REP3.annotation$end))


TOP1cc_REP1_genes_target.gr <- GRanges(seqnames = TOP1cc_REP1.annotation$seqnames, 
                     ranges = IRanges(start = TOP1cc_REP1.annotation$start, 
                                      end = TOP1cc_REP1.annotation$end))
TOP1cc_REP2_genes_target.gr <- GRanges(seqnames = TOP1cc_REP2.annotation$seqnames, 
                                  ranges = IRanges(start = TOP1cc_REP2.annotation$start, 
                                                   end = TOP1cc_REP2.annotation$end))


# Use findOverlapsOfPeaks function (ChIPpeakAnno) for verlapping peaks detection
library(ChIPpeakAnno)
ol_XRCC4_target_genes <-findOverlapsOfPeaks(XRCC4_REP1_genes_target.gr,
                                            XRCC4_REP2_genes_target.gr,
                                            XRCC4_REP3_genes_target.gr,
                                            maxgap=1000)
ol_TOP1_target_genes <-findOverlapsOfPeaks(TOP1cc_REP1_genes_target.gr,
                               TOP1cc_REP2_genes_target.gr,
                               maxgap=1000) 

# Venn Diagram representation and statistical analysis of intersection significance
result_XRCC4_target_genes <-makeVennDiagram(ol_XRCC4_target_genes, totalTest=20203 ,
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
result_XRCC4_target_genes$p.value


# Venn Diagram representation and statistical analysis of intersection significance
result_TOP1 <-makeVennDiagram(ol_TOP1_target_genes, totalTest=20203, #20203 protein-coding genes GRCh38
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




########################################################################
# Graphical representation of annotated peaks related to gene features #
########################################################################

##library(ChIPseeker)
## This library provides a set of functions for graphical analysis of peaks

# XRCC4 ChIP-seq : peaks representation related to genes features

plotAnnoPie(XRCC4_CPT.peakAnno)
plotAnnoBar(XRCC4_CPT.peakAnno)
plotDistToTSS(XRCC4_CPT.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")
upsetplot(XRCC4_CPT.peakAnno)

plotPeakProf2(c(XRCC4_CPT_REP1.peaks,XRCC4_CPT_REP2.peaks,XRCC4_CPT_REP3.peaks), 
              upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)

# TOP1-seq: peaks representation related to genes features

plotAnnoPie(TOP1cc.peakAnno)
plotAnnoBar(TOP1cc.peakAnno)
plotDistToTSS(TOP1cc.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")
upsetplot(TOP1cc.peakAnno)

plotPeakProf2(peak = c(TOP1cc_REP1.peaks,TOP1cc_REP2.peaks), upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)

# TOP1vsXRCC4: graphical comparison

# Create a list of annotated peaks
peakAnnoList <- list(XRCC4 = XRCC4_CPT.peakAnno, TOP1 = TOP1cc.peakAnno)
# Peak feature distribution plotAnnoBar function
plotAnnoBar(peakAnnoList)

# Create a list of peaks bed files
tss_list <- list(XRCC4 = c(XRCC4_CPT_REP1.peaks,XRCC4_CPT_REP2.peaks,XRCC4_CPT_REP3.peaks), 
                 TOP1 = c(TOP1cc_REP1.peaks,TOP1cc_REP2.peaks))
# Peak distribution relative to TSS
plotPeakProf2(tss_list,
              upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)

########################################################
#   GEO enrichment analysis: GEO, KEGG, Reactome, DOSE #
########################################################

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

## Save the target genes selection of consensus peaks
write(x = target.genes.XRCC4_CPT,file = "XRCC4_CPT_target_genes.txt")
write(x = target.genes.TOP1cc,file = "TOP1cc_target_genes.txt")

# GO: gene ontology enrichment analysis (clusterProfiler)

XRCC4_CPT.enrich.go <- enrichGO(gene = target.genes_XRCC4_CPT,
                                OrgDb         = "org.Hs.eg.db",
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff = 0.05,
                                minGSSize = 10,
                                maxGSSize = 500,
                                readable = FALSE)

TOP1cc.enrich.go <- enrichGO(gene = target.genes_TOP1cc,
                             OrgDb         = "org.Hs.eg.db",
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05,
                             pvalueCutoff  = 0.05,
                             minGSSize = 10,
                             maxGSSize = 500,
                             readable = FALSE)

# Graphical representation
dotplot(XRCC4_CPT.enrich.go,showCategory = 20)
dotplot(TOP1cc.enrich.go,showCategory = 9, 
        title = "GO Enrichment Analysis") + theme(axis.text.x=element_text(size=10),
                                                  axis.text.y=element_text(size=9)) 

# KEGG Pathway enrichment analysis: clusterProfiler 
XRCC4_CPT.enrich.kegg <- enrichKEGG(gene  = target.genes_XRCC4_CPT,
                                    organism = "human",
                                    keyType = "kegg",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff = 0.05,
                                    minGSSize = 10,
                                    maxGSSize = 500)

TOP1cc.enrich.kegg <- enrichKEGG(gene  = target.genes_TOP1cc,
                                 organism = "human",
                                 keyType = "kegg",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff = 0.05,
                                 minGSSize = 10,
                                 maxGSSize = 500)

# Graphical representation
dotplot(XRCC4_CPT.enrich.kegg, showCategory = 20)
dotplot(TOP1cc.enrich.kegg, showCategory = 20)

# Reactome Pathway Enrichment Analysis
XRCC4_CPT.enrich.reactome <- enrichPathway(gene= target.genes_XRCC4_CPT,
                                          organism = "human",
                                          pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH",
                                          qvalueCutoff = 0.2,
                                          minGSSize = 10,
                                          maxGSSize = 500,
                                          readable = FALSE)

TOP1cc.enrich.reactome <- enrichPathway(gene= target.genes_TOP1cc,
                                       organism = "human",
                                       pvalueCutoff = 0.05,
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 0.05,
                                       minGSSize = 10, 
                                       maxGSSize = 500,
                                       readable = FALSE)
# Graphical representation
dotplot(XRCC4_CPT.enrich.reactome, showCategory = 10)
dotplot(TOP1cc.enrich.reactome, showCategory = 9,
        title ='Reactome') + theme(axis.text.x=element_text(size=10),
                                 axis.text.y=element_text(size=9))

