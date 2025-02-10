# Enrichment analysis 
# single-cell analysis package
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)
theme_set(theme_cowplot())
set.seed(12345)

# EnrichR
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')
# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)
# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

# Visualize enrichments #
# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2021", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)

# Marker gene overlap analysis
# It is necessary to determine the cell type
# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$spot_category
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1.5
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# Save the results in a CSV file
write.csv(markers, "markers.csv", row.names = FALSE)
write.csv(overlap_df, "overlap_df.csv", row.names = FALSE)
