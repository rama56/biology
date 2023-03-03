library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(dplyr)
source("data_getter.R")
source("matrix_operations.R")
source("plot_helper.R")

run = "Trial"
project = "HUSCH_blood"
output_file_1 = paste("outputs/", project, "/", run, ".rds", sep = "")

# 1. SET UP SEURAT OBJECT. 
#seurat_object = get_husch_blood_data()

data_getters = list(get_husch_blood_data, get_seurat_zero_impute_data)
for (data_getter in data_getters){
  # 1. SET UP SEURAT OBJECT.
  seurat_object= data_getter()
  # print(data_getter)
  seurat_object@project.name

  # 2. STANDARD PRE-PROCESSING FLOW.
  processed_object = run_initial_processes(seurat_object)
  
  # 3. ARTIFICIALLY DROP VALUES
  original_matrix = processed_object@assays$RNA@counts
  
  # p_list = list(0.01, 0.02, 0.05,  0.1, 0.2)
  p_list = list(0.00)
  for (p in p_list) {
    cat("Dropout probability = ", p)
    artificial_matrix = drop_values_probabilistic(original_matrix, p)
    artificial_object = recreate_seurat_obj_from_matrix(artificial_matrix)
    
    # 3. RUN DATA IMPUTATION.
    imputed_object = RunALRA(artificial_object)
    imputed_matrix = imputed_object@assays$alra@data
    
    # 4. EVALUATE PERFORMANCE.
    distance_1 = find_distance(imputed_matrix, original_matrix)
    distance_2 = find_distance(imputed_matrix, artificial_matrix)
    distance_3 = find_distance(original_matrix, artificial_matrix)
    cat("Original, Imputed = ", distance_1)
    cat("Artificial, Imputed = ", distance_2)
    cat("Original, Artificial = ", distance_3)
  }
}

# SAVE FILE.
#saveRDS(nafld, file = output_file_1)
