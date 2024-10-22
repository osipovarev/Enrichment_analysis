#!/usr/bin/env Rscript

## Load packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(org.Hs.eg.db))

## Initiate argument parser
parser = ArgumentParser()
parser$add_argument("-gos", "--gos", type="character", help="list of GO terms to list genes")
args = parser$parse_args()

all_go <- unique(keys(org.Hs.eg.db, keytype ="GOALL"))
gos <- strsplit(args$gos, ",")[[1]]


get_genes_by_go_term <- function(gos) {
  go_genes <- select(org.Hs.eg.db, keytype="GOALL", keys=gos, columns= c("SYMBOL"))
  return(unique(go_genes$SYMBOL))
  }

for (gene in get_genes_by_go_term(gos)){
  cat(gene, "\n")
}