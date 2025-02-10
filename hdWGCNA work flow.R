gc()

suppressMessages(library(knitr))
suppressMessages(library(rmdformats))
suppressMessages(library(stringr))
suppressMessages(library(DT))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(shiny))
suppressMessages(library(ggplot2))
suppressMessages(library(extrafont))
suppressMessages(library(cowplot))
suppressMessages(library(magick))
suppressMessages(library(scales))
suppressMessages(library(viridis))
suppressMessages(library(tidyverse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(gplots))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(msigdbr))
suppressMessages(library(enrichplot))
suppressMessages(library(forcats))
library(ggnewscale)
library(patchwork)
library(WGCNA)
library(flashClust)
library(DESeq2)
library(RColorBrewer)
library(hdWGCNA)
gc()

custom_colors <- list()
colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)
colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)
colors_custom = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
                  '#49beaa', '#611c35', '#2708a0')
custom_colors$discrete <- unique(c(colors_dutch, colors_spanish, colors_custom))


# Read the Visium 10x h5 file, including spatial data
load(".../D1_Visium_Seurat.rda")
load(".../D2_Visium_Seurat.rda")

# ----merge seurat data ----- #
data_D1$region <- 'D1'
data_D2$region <- 'D2'
# merge into one seurat object
seurat_obj <- merge(data_D1, data_D2)
seurat_obj$region <- factor(as.character(seurat_obj$region), levels=c('D1', 'D2')) # seurat_obj = mergeしたvisium data

#perform Seurat/hdWGCNA analysis
seurat_obj <- seurat_obj %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures=4000) %>%
  ScaleData() 
#data_D1 = SCTransform(data_D1, assay = "Spatial") #generally unrecommended to use SCtransform data in hdWGCNA analysis
seurat_obj = RunPCA(seurat_obj, verbose = F)

# Louvain clustering and umap
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50)
seurat_obj <- FindClusters(seurat_obj,verbose = TRUE, resolution = 0.8)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:50)

# show the UMAP
p1 <- DimPlot(seurat_obj, label=TRUE, reduction = "umap", group.by = "seurat_clusters") + NoLegend()
p1
p2 <- SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3)
p2

##perform hdWGCNA analysis
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "vis"
)

seurat_obj <- MetaspotsByGroups(
  seurat_obj,
  group.by = c("orig.ident", "seurat_clusters"),
  ident.group = "seurat_clusters",
  assay = 'Spatial',
  slot = 'counts',
  min_spots=50
)
seurat_obj  <- NormalizeMetacells(seurat_obj)
m_obj <- GetMetacellObject(seurat_obj)

##perform hdWGCNA aganist metacells
seurat_obj <- SetDatExpr(
  seurat_obj,
  group.by= "seurat_clusters",
  group_name = levels(seurat_obj@active.ident),
  assay="Spatial"
)
######

# test different soft power thresholds

seurat_obj <- TestSoftPowers(seurat_obj)
plot_list <- PlotSoftPowers(seurat_obj)

wrap_plots(plot_list, ncol=2)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name='test',
  overwrite_tom=TRUE,
  soft_power = 6 #change power to 6 that seems to be enough to acheve scale-free nature
)

# plot the dendrogram
PlotDendrogram(seurat_obj, main='Spatial hdWGCNA dendrogram')

seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "SM"
)

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
head(modules[,1:3])

# get module eigengenes and gene-module assignment tables
MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add the MEs to the seurat metadata so we can plot it with Seurat functions
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'spot_category', dot.min=0.1)
# Scaling data with a low number of groups may produce misleading results

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  xlab('') + ylab('')

p

#We can visualize the MEs directly on the spatial coordinates using SpatialFeaturePlot.
p <- SpatialFeaturePlot(
  seurat_obj,
  features = mods,
  alpha = c(0.1, 1),
  ncol = 4
)

png("MEs_featureplot.png", height=4, width=5, units='in', res=300, bg="white")
p