

get_seurat_zero_impute_data <- function() {
  InstallData("pbmc3k")
  data("pbmc3k")
  # Initial processing
  return(pbmc3k)
}

get_husch_blood_data <- function() {
  project <- "HUSCH_blood"
  read_data_dir <- paste("../../Data/", project, "/", sep = "")
  data_file <- paste(read_data_dir,"HU_0042_Blood_10x_gene_count.h5", sep = "")
  pbmc.data <- Read10X_h5(filename = data_file)
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = project, min.cells = 3, min.features = 200)
  return(pbmc)
}

run_initial_processes <- function(seurat_object) {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # seurat_object <- SCTransform(seurat_object) 

  # seurat_object <- RunPCA(seurat_object)
  # seurat_object <- RunUMAP(seurat_object, dims = 1:30)
  return(seurat_object)
}

recreate_seurat_obj_from_matrix <- function(matrix) {
  project = "dummy"
  seurat_object <- CreateSeuratObject(counts = matrix, project = project)
  return(seurat_object)
}

