#!/usr/bin/env Rscript

### Simplify results of enrichment analysis ----------------------------------------------

## Load argument parser package
suppressPackageStartupMessages(library(argparse))

## Load clusterProfiler and related packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GOSemSim))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(simplifyEnrichment))

parser = ArgumentParser()
parser$add_argument("-w", "--workdir", type="character", help="abs path to work (current) dir")
parser$add_argument("-p", "--pcutoff", type="numeric", default=1, help="p-value cutoff on the input values; default=1")
parser$add_argument("-e", "--enrichtab", type="character", help="table (.tsv) from enrichGO in clusterProfiler")
parser$add_argument("-o", "--outfile", type="character", help=" to put output file with representative GOs")
args = parser$parse_args()

## Parse arguments from command line
curr_dir = args$workdir
setwd(curr_dir)

file_name = args$enrichtab
out_file = args$outfile
pval_thresh = args$pcutoff

# Object with ancestor GO terms
GO_anc = as.list(GOBPANCESTOR)

# Object to match GO ID with GO definitions
GO_def = suppressMessages(AnnotationDbi::select(GO.db,
                                                keys=(keys(GO.db)),
                                                columns=c("GOID","TERM","DEFINITION"),
                                                keytype="GOID"))

# Object with semantic similarity and information content
GO_semsim = godata("org.Hs.eg.db", ont="BP")

df <- read.csv(file_name, header=TRUE, sep='\t')
df_sign <- df[df$pvalue < pval_thresh, ]

### Subir's approach
mat_all <- GO_similarity(df_sign$ID, ont = "BP", db = 'org.Hs.eg.db')
go_df <- simplifyGO(mat_all, plot = FALSE)

### Retrieve representative GO term for each cluster
repr_GO = 
  
  # Retrieve ancestral GO terms for each GO term in go_id
  GO_anc[df_sign$ID] %>% stack() %>% `colnames<-`(c("ancestral", "go_id")) %>% 
  
  # Add clusters from simplifyGO and count the number of times each ancestral GO appears in ancestral terms for each GO cluster
  left_join(dplyr::select(go_df, c("cluster", "id")), by = setNames("id", "go_id")) %>% 
  group_by(cluster, ancestral) %>% add_count() %>% ungroup %>% 
  
  # Add Information Content for each ancestral GO and calculate "importance" (most informative somewhat common ancestral GO term)
  mutate(IC = GO_semsim@IC[ancestral]) %>% 
  mutate(importance = (n*IC^2)) %>% 
  group_by(cluster) %>% 
  dplyr::slice(which.max(importance)) %>% 
  
  # Clean output
  left_join(dplyr::select(GO_def, c("GOID","TERM")), by = setNames("GOID", "ancestral")) %>% 
  dplyr::select(c("ancestral","cluster","importance","TERM"))


# Function to find matching elements
find_ids <- function(n, df) {
  col_values <- df$id[df$id %in% df$id[df$cluster == n]]
  # return(paste(elements, collapse = ","))
  return(col_values)
}

find_pvalues <- function(n, df) {
  col_values <- df$pvalue[df$ID %in% unlist(repr_GO$gos[repr_GO$cluster == n])]
  # return(paste(elements, collapse = ","))
  return(col_values)
}

# Add new columns with all p-values and all GOs 
repr_GO$gos <- sapply(repr_GO$cluster, find_ids, go_df)
repr_GO$pvalues <- sapply(repr_GO$cluster, find_pvalues, df_sign)
repr_GO$minpval <- sapply(repr_GO$pvalues, min)

# Concatenate each list of p-values and GOs into a comma-separated string
repr_GO$pvalues <- sapply(repr_GO$pvalues, function(x) format(x, digits=4))
repr_GO$pvalues <- sapply(repr_GO$pvalues, function(x) paste(x, collapse = ","))
repr_GO$gos <- sapply(repr_GO$gos, function(x) paste(x, collapse = ","))

## Write results to a file
write.table(repr_GO, out_file, sep='\t', row.names=FALSE)


