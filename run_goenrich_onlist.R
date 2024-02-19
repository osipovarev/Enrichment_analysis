#!/usr/bin/env Rscript

## Enrichment analysis ---------------------------------------------------------

## Load argument parser package
suppressPackageStartupMessages(library(argparse))

## Load clusterProfiler and related packages
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ggplot2))


## Initiate argument parser
parser = ArgumentParser()
parser$add_argument("-w", "--workdir", type="character", help="abs path to work (current) dir")
# parser$add_argument("-p", "--pcutoff", type="numeric", default=0.1, help="p-value cutoff on the input values; default=0.1")
parser$add_argument("-g", "--genetab", type="character", help="table (.tsv) from Degenotate MK test")
parser$add_argument("-d", "--disttop", type="character", help="list of genes with top values")
parser$add_argument("-o", "--outdir", type="character", help="output dir name (adj p-value cutoff <.01) to put files .enrichGO.tsv")
args = parser$parse_args()

## Parse arguments from command line
curr_dir = args$workdir
setwd(curr_dir)
file_name = args$genetab
enrich_outdir = args$outdir
disttop = args$disttop
# pval_thresh = args$pcutoff

file_df <- read.csv(file_name, header=TRUE, sep='\t')

## Perform gene set enrichment analysis (clusterProfiler). DoS > 0
genes = unlist(read.csv(disttop, header=FALSE, sep='\t'))

enrich_res <- enrichGO(
  gene          = genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  pAdjustMethod = "BH",
  minGSSize     = 40,
  maxGSSize     = 400,
  universe = na.omit(file_df[(file_df$mk.raw.p.value != 1), ])$gene
  # universe = file_df$gene
)
  
  # enrich_result = enrich_res@result
  # enrich_result = gofilter(enrich_res, level=4)
enrich_filt = simplify(enrich_res, cutoff=0.5, by="p.adjust", select_fun=min)
enrich_result = enrich_filt@result

enrich_out_file = paste0(enrich_outdir, 'pos.top5.dist.enrichGO.tsv')
write.table(enrich_result[enrich_result$pvalue < 0.01, ], row.names=F, quote=F, file=enrich_out_file, sep="\t")
