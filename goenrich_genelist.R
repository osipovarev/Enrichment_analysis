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
parser$add_argument("-g", "--genes", type="character", help="gene list (newline-separated)")
parser$add_argument("-u", "--universe", type="character", help="gene list (newline-separated) \
                      to be used as universe in the enrichment analysis")
parser$add_argument("-o", "--outfile", type="character", help="output file name (adj p-value cutoff <.01)")
args = parser$parse_args()

## Parse arguments from command line
curr_dir = args$workdir
setwd(curr_dir)
file_name = args$genes
enrich_outfile = args$outfile
universe_file = args$universe


genes <- read.csv(file_name, header=FALSE)[[1]]

## Perform gene set enrichment analysis (clusterProfiler)
if (is.null(universe_file)) {
  print('background gene set is not provided => using all')
  
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
    )

} else {
  print('background gene set is provided')
  universe <- read.csv(universe_file, header=FALSE)[[1]]
  
  enrich_res <- enrichGO(
    gene          = genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    pAdjustMethod = "BH",
    # minGSSize     = 40,
    # maxGSSize     = 400,
    # universe = na.omit(file_df[(file_df$mk.raw.p.value != 1), ])$gene
    universe = universe
    )
}



# enrich_filt = gofilter(enrich_res, level=5)
enrich_filt = enrich_res
enrich_filt = simplify(enrich_res, cutoff=0.5, by="p.adjust", select_fun=min)
enrich_result = enrich_filt@result

# enrich_out_file = paste0(enrich_outdir, dos, '.enrichGO.tsv')
write.table(enrich_result[enrich_result$pvalue < 0.05, ], row.names=F, quote=F, file=enrich_outfile, sep="\t")
# write.table(enrich_result, row.names=F, quote=F, file=enrich_outfile, sep="\t")
