## =========================================
## High Complexity Barcode Library Analysis
## 6.5.2017 CH, Yale U
## 10.3.2018 MM, Yale U
## =========================================
rm(list = ls()); gc()

library(Biobase)
library(Biostrings)
library(ShortRead)
library(ggplot2)

# !!! PLEASE SET WORKING DIRECTORY BEFORE RUNING
setwd("")
source("scripts/HCBarcodePipeline.R")

## ............................................
## ... 1. Read in Cellecta barcode sequences

ct.bcodes18.rev <-
  read.delim("data/Cellecta-CellTracker-13Kx13K-Lentiviral-Barcode-Library-Sequences-Rev.txt",
    stringsAsFactors = FALSE,
    header = FALSE
  )[, 1]
ct.bcodes18.rev <- DNAStringSet(ct.bcodes18.rev)

# Preprocess dictionary
bc18dict <- PDict(ct.bcodes18.rev, max.mismatch = 1)

## ............................................
## ... 2. Read in FASTQ sequence

# list all files in data folder
fls <- dir("data/", "*fastq.gz", full = TRUE)

# analyze samples one by one
for (a in 1:length(fls)) {
  name <- strsplit(fls[a],split="/")[[1]]
  name <- strsplit(name[length(name)],split=".fastq")[[1]][1]
  cat("Sample:", name, "\n")
  
  cat("Reading FASTQ files...\n")
  data <- readFastq(fls[a])
  
  # pre-processing
  res <- runBCPipeline(data)
  
  # save results
  saveRDS(res, file = paste0("res/", name, "_barcodes.rds"))
  
  # plot distribution of barcode abundance per sample
  pdf(file=paste0("res/bc_distr_",name,".pdf"),width=7.5,height=4)
  par(mar = c(4,4,2,0))
  plot(res, main = paste0(name," - ",res$n["nas"], " reads, ",
                          res$n["nubc"]," barcodes"), log="x")
  dev.off()
  
  # plot distribution of barcode abundance in percents per sample
  pdf(file=paste0("res/bc_distr_",name,"_perc.pdf"),width=7.5,height=4)
  par(mar = c(4,4,2,0))
  plot(res, main = paste0(name," - ",res$n["nas"], " reads, ",
                          res$n["nubc"]," barcodes"), scale = "percent", log="x")
  dev.off()
}
