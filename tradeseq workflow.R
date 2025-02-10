### The tradeSeq workflow ###
## https://statomics.github.io/tradeSeq/articles/tradeSeq.html ##
gc()
setwd("/Users/niitsutakayuki/Desktop")

# set up #
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(princurve)
install.packages("BiocManager")
BiocManager::install("slingshot", update = TRUE, ask = FALSE)
# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))

# data loading, clustering #
ILDs <- readRDS("/GSE135893F.rds") # Seurat data filtered for the fibroblast cell type
DimPlot(ILDs, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.5, label.size = 5, label=FALSE, raster=FALSE)

## data extract from seurat object ##
# UMAP, clustering data
rd <- ILDs@reductions$umap@cell.embeddings # The results of UMAP
cl <- ILDs@active.ident # The results of clustering (ref. seurat_cluster)

# sds (crv) data 
# quiescentmarker =  INMT, NPNT, in this public scRNA-seq data, expressed in Myofibroblasts
sds <- slingshot(rd, clusterLabels = cl, start.clus = '4', end.clus = '3')
lin <- getLineages(rd, clusterLabels = cl, start.clus = '4', end.clus ='3')
crv <- getCurves(lin)

countMatrix <- as(ILDs@assays$SCT@counts, "dgCMatrix")
counts <- as.matrix(countMatrix)

## Fit negative binomial model ##
# refer to the fitGAM vignette  https://statomics.github.io/tradeSeq/articles/fitGAM.html
# Additional to add additional covariates to the model, speed up computation or allow for custom normalization, amongst others.
# evaluateK function 
# This takes a little time to run　
# knot from 3 to 20
options(future.globals.maxSize = 96 * 1024^3)  # 96GB
# the output from evaluateK, we refer users to the fitGAM vignette
# More knots will allow more flexibility, but also increase the risk of overfitting
set.seed(5)
# Parallel computing
BPPARAM <- BiocParallel::bpparam()
print(BiocParallel::bpparam())
BPPARAM <- BiocParallel::MulticoreParam(workers = 12)  
print(BPPARAM)
icMat <- evaluateK(counts = counts, sds = sds, k = 16:20, 
                   nGenes = 200, verbose = TRUE, BPPARAM = BPPARAM)
print(icMat)


set.seed(7)
pseudotime <- slingPseudotime(sds, na = FALSE)
cellWeights <- slingCurveWeights(sds)
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 20, verbose = FALSE, parallel=TRUE, BPPARAM = BPPARAM)

# check the convergence of each gene
table(rowData(sce)$tradeSeq$converged)

## Within-lineage comparisons ##
# Association of gene expression with pseudotime
assoRes <- associationTest(sce)
head(assoRes)
# Discovering progenitor marker genes
startRes <- startVsEndTest(sce, lineages=TRUE)
# testing against fold change threshold of 2
start2 <- startVsEndTest(sce, l2fc = log2(2))

# visualize the estimated smoothers for the most significant gene
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[1]] 
# Alternatively, color the cells in UMAP space with that gene’s expression
plotSmoothers(sce, counts, gene = sigGeneStart)
plotGeneCount(sds, counts, gene = sigGeneStart)
## Comparing specific pseudotime values within a lineage
customRes <- startVsEndTest(sce, pseudotimeValues = c(0, 200))
customRes
customRes2 <- startVsEndTest(sce)
customRes2

## Between-lineage comparisons ##
# ①Discovering differentiated cell type markers
endRes <- diffEndTest(sce)
# the most significant gene using the plotSmoothers function
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[7]] #Ref. endRes
plotSmoothers(sce, counts, sigGene)
# UMAP space with that gene’s expression
plotGeneCount(sds, counts, gene = sigGene)

## ②Discovering genes with different expression patterns ##
# this number can be changed using the nPoints argument
patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])
plotSmoothers(sce, counts, gene = rownames(patternRes)[oPat][7])
plotGeneCount(sds, counts, gene = rownames(patternRes)[oPat][4])
# when ordered by pseudotime (default: 100), genes exhibiting the most statistically significant differences in expression between lineages are identified.

library(ggplot2)
patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes <- patternRes[, c("Gene", "pattern")]
endRes$Gene <- rownames(endRes)
endRes$end <- endRes$waldStat
endRes <- endRes[, c("Gene", "end")]
compare <- merge(patternRes, endRes, by = "Gene", all = FALSE)
compare$transientScore <- 
  rank(-compare$end, ties.method = "min")^2 + rank(compare$pattern, ties.method = "random")^2

# plot
ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = "patternTest Wald Statistic (log scale)",
       y = "diffEndTest Wald Statistic (log scale)") +
  scale_color_continuous(low = "yellow", high = "red") +
  theme_classic()
# visualize the expression in UMAP space of the top gene
topTransient <- compare[which.max(compare$transientScore), "Gene"]
plotSmoothers(sce, counts, gene = topTransient)
plotGeneCount(sds, counts, gene = topTransient)
head(
  compare[order(compare$transientScore, decreasing = TRUE), "Gene"],
  n = 10
)

### Early drivers of differentiation ###
# To find a list of genes that are differentially expressed between lineages at a particular region
# A method for assessing differences in average gene expression across specific regions, 
# particularly around points where two or more lineages diverge.
plotGeneCount(curve = sds, counts = counts,
              clusters = apply(slingClusterLabels(sds), 1, which.max),
              models = sce)

earlyDERes <- earlyDETest(sce, knots = c(1, 7))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])
plotGeneCount(sds, counts, gene = rownames(earlyDERes)[oEarly][1]) 


## Differential expression in large datasets ##
# testing against fold change threshold of 2
start2 <- startVsEndTest(sce, l2fc = log2(2))
# testing against fold change threshold of 1.5
pat2 <- patternTest(sce, l2fc = log2(1.5))

## Clustering of genes according to their expression pattern ##
# Extracting fitted values to use with any clustering method #
# Predict the tradeSeq model's fitted values for a specific gene, both at the single-cell level and along a pseudotime grid
yhat <- predictCells(models = sce, gene = "SPARC")
ysmooth <- predictSmooth(models = sce, gene = "SPARC", nPoints = 40)
write.csv(yhat, file="yhat_SPARC.csv")
write.csv(ysmooth, file="ysmooth_SPARC.csv")

## Clustering using RSEC, clusterExperiment ##
library(clusterExperiment)

# Clusters gene expression patterns from a sce using 20 discrete points. It selects the first 100 genes from the dataset
nPointsClus <- 20
clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus,
                                     genes = rownames(counts)[1:100])

# Save
write.csv(assoRes, file = "assoRes.csv", row.names = TRUE)
write.csv(startRes, file = "startRes.csv", row.names = TRUE)
write.csv(customRes, file = "customRes.csv", row.names = TRUE)
write.csv(customRes2, file = "customRes2.csv", row.names = TRUE)
write.csv(endRes, file = "endRes.csv", row.names = TRUE)
write.csv(patternRes, file = "patternRes.csv", row.names = TRUE)
write.csv(compare, file = "compare.csv", row.names = TRUE)
write.csv(earlyDERes, file = "earlyDERes.csv", row.names = TRUE)
write.csv(start2, file = "start2.csv", row.names = TRUE)
write.csv(pat2, file = "pat2.csv", row.names = TRUE)
write.csv(yhat, file = "yhat.csv", row.names = TRUE)
write.csv(ysmooth, file = "ysmooth.csv", row.names = TRUE)
