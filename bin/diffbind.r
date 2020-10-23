#!/usr/bin/env Rscript
library(DiffBind)
library(ChIPpeakAnno)

option_list <- list(make_option(c("-d", "--design"), type="character", default=NULL, help="filename of design table", metavar="path"),
                    make_option(c("-b", "--bams"), type="character", default=NULL, help="filename after sample name in featurecount file header e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'", metavar="string"),
                    make_option(c("-p", "--peaks"), type="character", default=NULL, help="peak files", metavar="string"),
                    make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$design)){
  print_help(opt_parser)
  stop("Please provide design table file.", call.=FALSE)
}
if (is.null(opt$bams)||is.null(opt$peaks)){
  print_help(opt_parser)
  stop("Please provide bam and peak file name.", call.=FALSE)
}

out <- "sample.csv"
## create a csv file with SampleID, Condition, Replicate, bamReads Peaks Peakcaller PeakFormat, ScoreCol, Factor, Tissue
sampleDesign <- read.csv(opt$design)
bamReads <- unlist(strsplit(opt$bams, "___"))
names(bamReads) <- sub(".mLb.clN.*.bam", "", bamReads)
Peaks <- unlist(strsplit(opt$peaks, "___"))
SampleID <- sub("_peaks.*?Peak", "", Peaks)
stopifnot(identical(names(bamReads), SampleID))
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
                      Peakcaller=Peakcaller,
                      PeakFormat=PeakFormat)
write.csv(samples, "sample.csv")
