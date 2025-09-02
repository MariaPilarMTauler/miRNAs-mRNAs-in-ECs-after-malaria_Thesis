#! /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources Rscript

## functions for rpkm and tpm
## from https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## from https://www.biostars.org/p/171766/
rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e9
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## ---- Load libraries
library(dplyr)
library(tidyr)

## ---- Read table "counts.txt" from the current working directory
input_file <- "counts.txt"
ftr.cnt    <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

## ---- Extract the file name (without extension) for output naming
tag <- tools::file_path_sans_ext(input_file)

## ---- Calculate RPKM
ftr.rpkm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(rpkm = rpkm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, rpkm)

write.table(
  ftr.rpkm,
  file      = paste0(tag, "_rpkm.txt"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

## ---- Calculate TPM
ftr.tpm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(tpm = tpm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, tpm)

write.table(
  ftr.tpm,
  file      = paste0(tag, "_tpm.txt"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
