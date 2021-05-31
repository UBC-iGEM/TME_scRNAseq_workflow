library(TMExplorer)
library(scater)
library(scran)
library(uwot)
library(plotly)

geo_accession_id = "GSE84465"

# load data
sce <- queryTME(geo_accession = geo_accession_id)[[1]]

# QUALITY CONTROL
sce <- addPerCellQC(sce)

libsize_drop <- isOutlier(sce$sum, nmads=3, type="lower", log=TRUE)
gene_drop <- isOutlier(sce$detected, nmads=3, type="lower", log=TRUE)

sce <- sce[, !(libsize_drop | gene_drop)]

# NORMALIZATION
set.seed(1111)

min_size <- min(150, floor(dim(sce)[2] * 0.3))
max_win <- min(101, min_size + 1)

clusters <- quickCluster(sce, assay.type="counts", method = "igraph", use.ranks = FALSE)
sce <- computeSumFactors(sce, assay.type="counts", sizes = seq(21, max_win, 5), min.mean = 0.1, clusters = clusters)

# Ensure that size factors are non-zero and non-negative before normalizing 
if(summary(sizeFactors(sce))[['Min.']] <= 0){
  stop('Negative size factors, cannot progess')
}

sce <- logNormCounts(sce)

# DIM REDUCTION 
set.seed(1010)
sce <- runPCA(sce, exprs_values = "logcounts")
sce <- runTSNE(sce, dimred = "PCA", exprs_values = "logcounts", ncomponents = 3, perplexity = 20)
sce <- runUMAP(sce, dimred = "PCA", exprs_values = "logcounts", ncomponents = 3, min_dist = 0.5, n_neighbors = 15, metric = "euclidean")


# CLUSTERING
g <- buildSNNGraph(sce, use.dimred = "PCA", k = 15, type = "number") # shared neighbor weighting
cluster <- igraph::cluster_louvain(g)$membership # extract clustering info from the graph that we just made 
sce$cluster <- factor(cluster) # add that info to the sce 

# save 
path_to_dir <- paste("/Users/amunzur/Documents/GitHub/iGEM_2021_project/data/clustered", geo_accession_id, sep = "/")
dir.create(file.path(path_to_dir))

saveRDS(sce, paste(path_to_dir, "sce.rds", sep = "/"))

















