library(dplyr)
library(Seurat)
library(patchwork)
library(future)
plan()

run = "_1"
type = "sample"
# read_data_dir = "../../../Data/nafld_mice/"
read_data_dir = "../../../Data/filtered_gene_bc_matrices/hg19/"
output_file_1 = paste("outputs/", type, run, ".rds", sep = "")
output_file_2 = paste("outputs/", type, run, "_final.rds", sep = "")
times = list()

times_vector = c()
labels_vector = c()


# 1. SET UP SEURAT OBJECT. 
# Load the PBMC dataset
s1 <- Sys.time()
pbmc.data <- Read10X(data.dir = read_data_dir)
e1 <- Sys.time()
t1 = list("Reading data", difftime(e1, s1, units="secs")[[1]])
print(t1)
times = c(times, t1)

labels_vector = c(labels_vector, t1[[1]])
times_vector = c(times_vector, t1[[2]])

# Initialize the Seurat object with the raw (non-normalized data).
s2 <- Sys.time()
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
e2 <- Sys.time()
t2 = list("Create Seurat object", difftime(e2, s2, units="secs")[[1]])
print(t2)
times = c(times, t2)

labels_vector = c(labels_vector, t2[[1]])
times_vector = c(times_vector, t2[[2]])

barplot(times_vector, main="Car Distribution", names.arg = c("a", "bc"),
        xlab="Number of Gears")

# 2. STANDARD PRE-PROCESSING FLOW.
s3 <- Sys.time()
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
e3 <- Sys.time()
t3=list("Percentage Feature Set", difftime(e3, s3, units="secs")[[1]])
print(t3)
times = c(times, t3)

s4 <- Sys.time()
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
e4 <- Sys.time()
t4=list("Subset", difftime(e4, s4, units="secs")[[1]])
print(t4)
times = c(times, t4)

# 3. NORMALIZE DATA.
s5 <- Sys.time()
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
e5 <- Sys.time()
t5=list("Normalize data (log normalize)", difftime(e5, s5, units="secs")[[1]])
print(t5)
times = c(times, t5)

s6 <- Sys.time()
pbmc <- NormalizeData(pbmc)
e6 <- Sys.time()
t6=list("Normalize data (default params)", difftime(e6, s6, units="secs")[[1]])
print(t6)
times = c(times, t6)

s7 <- Sys.time()
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
e7 <- Sys.time()
t7=list("Find variable features", difftime(e7, s7, units="secs")[[1]])
print(t7)
times = c(times, t7)

# 4. IDENTIFICATION OF HIGHLY VARIABLE FEATURES.

# 5. SCALE DATA.
all.genes <- rownames(pbmc)
s8 <- Sys.time()
pbmc <- ScaleData(pbmc, features = all.genes)
e8 <- Sys.time()
t8=list("Scale data (all genes)", difftime(e8, s8, units="secs")[[1]])
print(t8)
times = c(times, t8)

s9 <- Sys.time()
pbmc <- ScaleData(pbmc)
e9 <- Sys.time()
t9=list("Scale data (param less)", difftime(e9, s9, units="secs")[[1]])
print(t9)
times = c(times, t9)

s10 <- Sys.time()
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
e10 <- Sys.time()
t10=list("Scale data (vars to regress = percent.mt)", difftime(e10, s10, units="secs")[[1]])
print(t10)
times = c(times, t10)

# 5. PERFORM LINEAR DIMENSIONAL REDUCTION.
s11 <- Sys.time()
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
e11 <- Sys.time()
t11=list("Run PCA", difftime(e11, s11, units="secs")[[1]])
print(t11)
times = c(times, t11)

# 6. DETERMINE DIMENSIONALITY OF DATASET
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
s12 <- Sys.time()
pbmc <- JackStraw(pbmc, num.replicate = 100)
e12 <- Sys.time()
t12=list("Jack Straw", difftime(e12, s12, units="secs")[[1]])
print(t12)
times = c(times, t12)

s13 <- Sys.time()
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
e13 <- Sys.time()
t13=list("Score Jack Straw", difftime(e13, s13, units="secs")[[1]])
print(t13)
times = c(times, t13)

# 7. CLUSTER THE CELLS.
s14 <- Sys.time()
pbmc <- FindNeighbors(pbmc, dims = 1:10)
e14 <- Sys.time()
t14=list("Find Neighbours (dims = 1:10)", difftime(e14, s14, units="secs")[[1]])
print(t14)
times = c(times, t14)

s15 <- Sys.time()
pbmc <- FindClusters(pbmc, resolution = 0.5)
e15 <- Sys.time()
t15=list("Find Clusters", difftime(e15, s15, units="secs")[[1]])
print(t15)
times = c(times, t15)

# 8. NON-LINEAR DIMENSIONAL REDUCTION.
s16 <- Sys.time()
pbmc <- RunUMAP(pbmc, dims = 1:10)
e16 <- Sys.time()
t16=list("Run UMAP", difftime(e16, s16, units="secs")[[1]])
print(t16)
times = c(times, t16)

# SAVE FILE.
saveRDS(pbmc, file = output_file_1)

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
t17=list("Find all Markers", difftime(e17, s17, units="secs")[[1]])
print(t17)
times = c(times, t17)

pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

s18 <- Sys.time()
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
e18 <- Sys.time()
t18=list("Find Markers", difftime(e18, s18, units="secs")[[1]])
print(t18)
times = c(times, t18)

print(times)

saveRDS(pbmc, file = output_file_2)

