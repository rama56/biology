library(dplyr)
library(Seurat)
library(patchwork)


# 1. SET UP SEURAT OBJECT. 
s1 <- Sys.time()
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../../../Data/filtered_gene_bc_matrices/hg19/")
e1 <- Sys.time()
print(paste("Reading data = ", e1 - s1))

s2 <- Sys.time()
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
e2 <- Sys.time()
print(paste("Create Seurat object = ", e2 - s2))

pbmc

# 2. STANDARD PRE-PROCESSING FLOW.
s3 <- Sys.time()
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
e3 <- Sys.time()
print(paste("Percentage Feature Set = ", e3 - s3))

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

s4 <- Sys.time()
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
e4 <- Sys.time()
print(paste("Subset = ", e4 - s4))

# 3. NORMALIZE DATA.
s5 <- Sys.time()
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
e5 <- Sys.time()
print(paste("Normalize data (log normalize) = ", e5 - s5))

s6 <- Sys.time()
pbmc <- NormalizeData(pbmc)
e6 <- Sys.time()
print(paste("Normalize data (default params) = ", e6 - s6))

s7 <- Sys.time()
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
e7 <- Sys.time()
print(paste("Find variable features = ", e7 - s7))

# 4. IDENTIFICATION OF HIGHLY VARIABLE FEATURES.
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 5. SCALE DATA.
all.genes <- rownames(pbmc)
s8 <- Sys.time()
pbmc <- ScaleData(pbmc, features = all.genes)
e8 <- Sys.time()
print(paste("Scale data (all genes) = ", e8 - s8))

s9 <- Sys.time()
pbmc <- ScaleData(pbmc)
e9 <- Sys.time()
print(paste("Scale data (param less) = ", e9 - s9))

s10 <- Sys.time()
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
e10 <- Sys.time()
print(paste("Scale data (vars to regress = percent.mt) = ", e10 - s10))

# 5. PERFORM LINEAR DIMENSIONAL REDUCTION.
s11 <- Sys.time()
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
e11 <- Sys.time()
print(paste("Run PCA = ", e10 - s10))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# 6. DETERMINE DIMENSIONALITY OF DATASET
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
s12 <- Sys.time()
pbmc <- JackStraw(pbmc, num.replicate = 100)
e12 <- Sys.time()
print(paste("Jackk Straw = ", e12 - s12))

s13 <- Sys.time()
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
e13 <- Sys.time()
print(paste("Score Jack Straw = ", e13 - s13))

JackStrawPlot(pbmc, dims = 1:15)

# 7. CLUSTER THE CELLS.
s14 <- Sys.time()
pbmc <- FindNeighbors(pbmc, dims = 1:10)
e14 <- Sys.time()
print(paste("Find Neighbours (dims = 1:10) ", e14 - s14))

s15 <- Sys.time()
pbmc <- FindClusters(pbmc, resolution = 0.5)
e15 <- Sys.time()
print(paste("Find clusteres = ", e15 - s15))

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# 8. NON-LINEAR DIMENSIONAL REDUCTION.
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
s16 <- Sys.time()
pbmc <- RunUMAP(pbmc, dims = 1:10)
e16 <- Sys.time()
print(paste("Run UMAP = ", e16 - s16))
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

# SAVE FILE.
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

# FIND DIFFERENTIALLY EXPRESSED FEATURES.

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
s17 <- Sys.time()
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
e17 <- Sys.time()
print(paste("Find all Markers = ", e17 - s17))

pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

s18 <- Sys.time()
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
e18 <- Sys.time()
print(paste("Find Markers = ", e18 - s18))

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# ASSIGNING CELL TYPE IDENTITY TO CLUSTERS.
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "../output/pbmc3k_final.rds")

