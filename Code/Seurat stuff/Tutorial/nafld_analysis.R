library(dplyr)
library(Seurat)
library(patchwork)
#library(future)
#plan()

run = "Run"
project = "nafld"
read_data_dir = paste("../../../Data/nafld_mice/", sep = "")
output_file_1 = paste("outputs/", project, "/", run, ".rds", sep = "")
output_file_2 = paste("outputs/", project, "/", run, "_before_clustering.rds", sep = "")

# 1. SET UP SEURAT OBJECT. 
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = read_data_dir)
pbmc.data[c("mt-Nd1", "mt-Atp8", "MS4A1"), 1:30]
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = project, min.cells = 3, min.features = 200)

# 2. STANDARD PRE-PROCESSING FLOW.
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
head(pbmc@meta.data, 5)
nafld <- subset(pbmc, subset = orig.ident == "nafld")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 3. NORMALIZE DATA.
nafld <- NormalizeData(nafld, normalization.method = "LogNormalize", scale.factor = 10000)
nafld <- NormalizeData(nafld)

# 4. IDENTIFICATION OF HIGHLY VARIABLE FEATURES.
nafld <- FindVariableFeatures(nafld, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(nafld), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nafld)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 5. SCALE DATA.
all.genes <- rownames(nafld)
nafld <- ScaleData(nafld, features = all.genes)

nafld <- ScaleData(nafld)

nafld <- ScaleData(nafld, vars.to.regress = "percent.mt")

# 5. PERFORM LINEAR DIMENSIONAL REDUCTION.
nafld <- RunPCA(nafld, features = VariableFeatures(object = nafld))

# 6. DETERMINE DIMENSIONALITY OF DATASET
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
nafld <- JackStraw(nafld, num.replicate = 100)
nafld <- ScoreJackStraw(nafld, dims = 1:20)
 
ElbowPlot(object = nafld,ndims = 30)

saveRDS(nafld, file = output_file_2)

# 7. CLUSTER THE CELLS.
nafld <- FindNeighbors(nafld, dims = 1:3)
nafld <- FindClusters(nafld, resolution = 0.5)

# 8. NON-LINEAR DIMENSIONAL REDUCTION.
nafld <- RunUMAP(nafld, dims = 1:18)
nafld <- RunTSNE(nafld, dims = 1:18)

# 7. CLUSTER THE CELLS.
nafld <- FindNeighbors(nafld, dims = 1:18)
nafld <- FindClusters(nafld, resolution = 0.5)

DimPlot(object = nafld, reduction = "umap")
DimPlot(object = nafld, reduction = "tsne")
FeaturePlot(nafld,features = c("Tmsb4x","Malat1"), slot="data", reduction = "tsne")
FeaturePlot(nafld,features = c("Gm26917","Neat1", "Glul", "Cd52", "S100a9", "S100a6"), slot="data", reduction = "tsne")

# SAVE FILE.
saveRDS(nafld, file = output_file_1)

# 9. FIND DIFFERENTIALLY EXPRESSED FEATURES.

# find all markers of cluster 2
cluster2.markers <- FindMarkers(nafld, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(nafld, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
nafld.markers <- FindAllMarkers(nafld, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

nafld.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(nafld, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# SAVE FILE.

saveRDS(nafld, file = output_file_2)

