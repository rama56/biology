
library(dplyr)
library(Seurat)
library(patchwork)
library(future)
plan()

run = "_1"
type = "sample"
# read_data_dir = "../../../Data/nafld_mice/"
read_data_dir = "../../../Data/filtered_gene_bc_matrices/hg19/"

input_file_2 = paste("outputs/", type, run, "_final.rds", sep = "")

pbmc <- readRDS(input_file_2)

timing.comparisons <- data.frame(fxn = character(), time = numeric(), strategy = character())
plan("sequential")
start <- Sys.time()
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
end <- Sys.time()
timing.comparisons <- rbind(timing.comparisons, data.frame(fxn = "ScaleData", time = as.numeric(end -
                                                                                                  start, units = "secs"), strategy = "sequential"))

start <- Sys.time()
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)
end <- Sys.time()
timing.comparisons <- rbind(timing.comparisons, data.frame(fxn = "FindMarkers", time = as.numeric(end -
                                                                                                    start, units = "secs"), strategy = "sequential"))

plan("multiprocess", workers = 4)
start <- Sys.time()
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
end <- Sys.time()
timing.comparisons <- rbind(timing.comparisons, data.frame(fxn = "ScaleData", time = as.numeric(end -
                                                                                                  start, units = "secs"), strategy = "multiprocess"))

start <- Sys.time()
markers <- FindMarkers(pbmc, ident.1 = "NK", verbose = FALSE)
end <- Sys.time()
timing.comparisons <- rbind(timing.comparisons, data.frame(fxn = "FindMarkers", time = as.numeric(end -start, units = "secs"), strategy = "multiprocess"))
                                                                                                    
library(ggplot2)
library(cowplot)
ggplot(timing.comparisons, aes(fxn, time)) + geom_bar(aes(fill = strategy), stat = "identity", position = "dodge") +
ylab("Time(s)") + xlab("Function") + theme_cowplot()                                                                                              
