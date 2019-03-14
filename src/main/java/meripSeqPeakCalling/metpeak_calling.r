# MetPeak½øĞĞpeakcalling
Args <- commandArgs()
library(MeTPeak)

gtf <- Args[6]
ipFileDir <- Args[7]
inputFileDir <- Args[8]
outputDir <- Args[9]
experimenName <- Args[10]

setwd(outputDir)

ipFiles <- as.vector(list.files(ipFileDir, pattern="*.bam$", full.names=T))
inputFiles <- as.vector(list.files(inputFileDir, pattern="*.bam$", full.names=T))

metpeak(GENE_ANNO_GTF=gtf,IP_BAM=ipFiles,INPUT_BAM=inputFiles, EXPERIMENT_NAME="example")