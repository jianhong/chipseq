#!/usr/bin/env Rscript
library(DiffBind)
library(ChIPpeakAnno)
library(rtracklayer)
library(ggplot2)
library(GenomicFeatures)
library(optparse)

option_list <- list(make_option(c("-d", "--design"), type="character", default=NULL, help="filename of design table", metavar="path"),
                    make_option(c("-p", "--peaks"), type="character", default=NULL, help="peak files", metavar="string"),
                    make_option(c("-g", "--gtf"), type="character", default=NULL, help="filename of gtf file", metavar="path"),
                    make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$design)){
  print_help(opt_parser)
  stop("Please provide design table file.", call.=FALSE)
}
if (is.null(opt$peaks)){
  print_help(opt_parser)
  stop("Please provide bam and peak file name.", call.=FALSE)
}
if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide gtf file.", call.=FALSE)
}
out <- "sample.csv"
## create a csv file with SampleID, Condition, Replicate, bamReads Peaks Peakcaller PeakFormat, ScoreCol, Factor, Tissue
sampleDesign <- read.csv(opt$design)
opt$peaks <- unlist(strsplit(opt$peaks, "___"))
bamReads <- opt$peaks[grepl("bam$", opt$peaks)]
names(bamReads) <- sub(".mLb.clN.*bam", "", bamReads)
Peaks <- opt$peaks[grepl("Peak$", opt$peaks)]
SampleID <- sub("_peaks.*?Peak", "", Peaks)
if(!all(names(bamReads) %in% SampleID)){
  stop("some names of bamReads is not in SampleID. names(bamReads)=", 
       paste(names(bamReads), collapse = "; "),
       "SampleID=", paste(SampleID, collapse=";"))
}
bamReads <- bamReads[SampleID]
names(Peaks) <- SampleID
rownames(sampleDesign) <- paste(sampleDesign$group, sampleDesign$replicate, sep="_R")
Condition <- sampleDesign[SampleID, "group"]
Replicate <- sampleDesign[SampleID, "replicate"]
Factor <- sampleDesign[SampleID, "antibody"]
Peakcaller <- "macs2"
PeakFormat <- sub("^.*?_peaks.(.*)$", "\\1", Peaks)

samples <- data.frame(SampleID=SampleID,
                      Condition=Condition,
                      Replicate=Replicate,
                      Factor=Factor,
                      bamReads=bamReads,
                      Peaks=Peaks,
                      Peakcaller=Peakcaller,
                      PeakFormat=PeakFormat,
                      ScoreCol=5)
pf <- "DiffBind"
dir.create(pf)
l <- lapply(Peaks, readLines)
keep <- lengths(l)>0
samples <- samples[keep, , drop=FALSE]
if(nrow(samples)>3){
  write.csv(samples, file.path(pf, "sample.csv"))
  
  chip <- dba(sampleSheet = file.path(pf, "sample.csv"))
  pdf(file.path(pf, "DiffBind.sample.correlation.pdf"), width = 9, height = 9)
  plot(chip)
  dev.off()
  pdf(file.path(pf, "DiffBind.PCA.plot.pdf"))
  dba.plotPCA(chip, DBA_CONDITION, label=DBA_ID)
  dev.off()
  png(file.path(pf, "DiffBind.sample.correlation.png"))
  plot(chip)
  dev.off()
  png(file.path(pf, "DiffBind.PCA.plot.png"))
  dba.plotPCA(chip, DBA_CONDITION, label=DBA_ID)
  dev.off()
  
  chip <- dba.count(chip, bLog=TRUE)
  saveRDS(chip, file.path(pf, "chip.rds"))
  chip.bk <- chip
  
  txdb <- makeTxDbFromGFF(opt$gtf)
  gtf <- import(opt$gtf)
  id2symbol <- function(gtf){
    if(is.null(gtf$gene_name)) return(NULL)
    x <- data.frame(id=gtf$gene_id, symbol=gtf$gene_name)
    x <- unique(x)
    x <- x[!duplicated(x$id), ]
    x <- x[!is.na(x$id), , drop=FALSE]
    if(nrow(x)==0) return(NULL)
    y <- x$symbol
    names(y) <- x$id
    y
  }
  id2symbol <- id2symbol(gtf)
  anno <- toGRanges(txdb)
  resList <- list()
  if(length(unique(Condition))>=2){
    contrasts <- combn(unique(Condition), m = 2, simplify = FALSE)
    names(contrasts) <- sapply(contrasts, function(.ele) paste(.ele[2], .ele[1], sep="-"))
    for(i in seq_along(contrasts)){
      gp1 <- contrasts[[i]][1]
      gp2 <- contrasts[[i]][2]
      if(sum(chip.bk$samples$Condition %in% c(gp1, gp2))<3){
        next
      }
      chip <- chip.bk
      chip <- dba.contrast(chip,
                           group1 = dba.mask(chip, DBA_CONDITION, gp1),
                           group2 = dba.mask(chip, DBA_CONDITION, gp2),
                           categories = DBA_CONDITION)
      chip <- dba.analyze(chip, bFullLibrarySize=FALSE)
      chip.DB <- dba.report(chip, th=1)
      
      # Annotation
      chip.anno <- annotatePeakInBatch(chip.DB, AnnotationData = anno,
                                       output = "nearestLocation",
                                       PeakLocForDistance = "middle",
                                       FeatureLocForDistance = "TSS",
                                       ignore.strand = TRUE)
      if(length(id2symbol)>0) chip.anno$symbol[!is.na(chip.anno$feature)] <- id2symbol[chip.anno$feature[!is.na(chip.anno$feature)]]
      resList[[names(contrasts)[i]]] <- chip.anno[chip.anno$FDR<0.05]
      chip.m <- as.data.frame(unname(chip.anno),
                              stringsAsFactor=FALSE)
      write.csv(chip.m, file.path(pf, paste0("DiffBind.res.", names(contrasts)[i], ".all.csv")))
      chip.m <- chip.m[chip.m$FDR<0.05, ]
      write.csv(chip.m, file.path(pf, paste0("DiffBind.res.", names(contrasts)[i], ".FDR.05.xls")))
      chip.anno.b <- chip.anno[chip.anno$FDR<0.05]
      mcols(chip.anno.b) <- DataFrame(score=chip.anno.b$Fold)
      export(chip.anno.b, file.path(pf, paste0("DiffBind.chip.res.", names(contrasts)[i], ".FDR.05.bedGraph")))
      # plots
      pdf(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".MA.plot.pdf")))
      dba.plotMA(chip)
      dev.off()
      
      pdf(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".Volcano.plot.pdf")))
      dba.plotVolcano(chip, bUsePval = TRUE)
      dev.off()
      
      png(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".MA.plot.png")))
      dba.plotMA(chip)
      dev.off()
      
      png(file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".Volcano.plot.png")))
      dba.plotVolcano(chip, bUsePval = TRUE)
      dev.off()
      
      # export counts table
      counts <- dba.peakset(chip, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)
      write.csv(counts, file.path(pf, paste0("DiffBind.", names(contrasts)[i], ".counts.csv")))
    }
  }
  resList <- if(length(resList)>1) GRangesList(resList) else if(length(resList)>0) resList[[1]]
  
  if(packageVersion("ChIPpeakAnno")>="3.23.12"){
    if(length(resList)>0){
      out <- genomicElementDistribution(resList, 
                                        TxDb = txdb,
                                        promoterRegion=c(upstream=2000, downstream=500),
                                        geneDownstream=c(upstream=0, downstream=2000),
                                        promoterLevel=list(
                                          # from 5' -> 3', fixed precedence 3' -> 5'
                                          breaks = c(-2000, -1000, -500, 0, 500),
                                          labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                     "upstream <500b", "TSS - 500b"),
                                          colors = c("#FFE5CC", "#FFCA99", 
                                                     "#FFAD65", "#FF8E32")),
                                        plot = FALSE)
      
      ggsave(file.path(pf, "genomicElementDistribuitonOfDiffBind.pdf"), plot=out$plot, width=9, height=9)
      ggsave(file.path(pf, "genomicElementDistribuitonOfDiffBind.png"), plot=out$plot)
      out <- metagenePlot(resList, txdb)
      ggsave(file.path(pf, "metagenePlotToTSSofDiffBind.pdf"), plot=out, width=9, height=9)
      ggsave(file.path(pf, "metagenePlotToTSSofDiffBind.png"), plot=out)
    }
    if(length(samples$Peaks)>0){
      peaks <- mapply(samples$Peaks, samples$PeakFormat, 
                      FUN=function(.ele, .format) toGRanges(.ele, format=.format), 
                      SIMPLIFY = FALSE)
      names(peaks) <- samples$SampleID
      peaks <- GRangesList(peaks)
      out <- genomicElementDistribution(peaks, 
                                        TxDb = txdb,
                                        promoterRegion=c(upstream=2000, downstream=500),
                                        geneDownstream=c(upstream=0, downstream=2000),
                                        promoterLevel=list(
                                          # from 5' -> 3', fixed precedence 3' -> 5'
                                          breaks = c(-2000, -1000, -500, 0, 500),
                                          labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                     "upstream <500b", "TSS - 500b"),
                                          colors = c("#FFE5CC", "#FFCA99", 
                                                     "#FFAD65", "#FF8E32")),
                                        plot = FALSE)
      
      ggsave(file.path(pf, "genomicElementDistribuitonOfEachPeakList.pdf"), plot=out$plot, width=9, height=9)
      ggsave(file.path(pf, "genomicElementDistribuitonOfEachPeakList.png"), plot=out$plot)
      
      out <- metagenePlot(peaks, txdb)
      ggsave(file.path(pf, "metagenePlotToTSSOfEachPeakList.pdf"), plot=out, width=9, height=9)
      ggsave(file.path(pf, "metagenePlotToTSSOfEachPeakList.png"), plot=out)
      
      if(length(peaks)<=5){
        ol <- findOverlapsOfPeaks(peaks)
        png(file.path(pf, "vennDiagram.all.png"))
        makeVennDiagram(ol, connectedPeaks="keepAll")
        dev.off()
      }else{
        peaks.split <- split(peaks, samples$Condition)
        null <- mapply(peaks.split, names(peaks.split), FUN=function(.peaks, .name){
          if(length(.peaks)<=5){
            ol <- findOverlapsOfPeaks(.peaks)
            png(file.path(pf, paste0("vennDiagram.", .name, ".png")))
            makeVennDiagram(ol, connectedPeaks="keepAll")
            dev.off()
          }
        })
      }
    }
  }
}



