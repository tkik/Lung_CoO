library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(harmony)
library(cerebroApp)
library(biomaRt)
library(future)
library(qs)
###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
PATH ="V:/"} else {
PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
PATH_Y="N:/"} else {
PATH_Y="/C010-datasets/"}
PATH_ICGC <- "/icgc/"
plan("multiprocess", workers = 2)
options(future.globals.maxSize= 5000*1024^2)
DATA = paste0(PATH, "Reka/33_CoO_lung/scRNASeq/data")
DATA_datasets = paste0(PATH_Y, "External/Sotillo_mouse_lung_cell/scRNASeq_CoO_lung/scrnaseq_analysis/")
RESULT = "./output/"
CODE = "./code/"
CALLS = paste0(PATH_Y, "External/Sotillo_mouse_lung_cell/scRNASeq_CoO_lung/cellranger_output_30_07/")
CALLS2 = paste0(PATH_Y, "External/Sotillo_mouse_lung_cell/scRNASeq_CoO_lung/cellranger_output_22_11/")
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
humanx <- unique(genesV2[, 2])
# Print the first 6 genes found to the screen
return(humanx)
}

source(file.path(CODE,"/hi_res_clusters_functions.R"), echo=TRUE)

seu.dat <- readRDS(file=file.path(DATA, "harmony_integrated_slingshot_3rd_no_ciliated.RDS"))
seu.dat2 <- FindClusters_highres(seu.dat, resolution = 2)

seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==13] <- "Pre-tumor stage"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==12] <- "AT1/AT2-like tumor"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==11] <- "Club cells"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==10] <- "H2-K1 high club-like progenitors"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==9] <- "Proliferating cells"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==8] <- "Tumor"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==6] <- "AT2 cells"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==4] <- "Club cells"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==3] <- "Club cells"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==2] <- "Club cells"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==1] <- "Tumor"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==0] <- "Club cells"

seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==7] <- "Club-like progenitors"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==7 & seu.dat2$seurat_clusters_high_res %in% c("15", "3", "13")] <- "AT2-like progenitors"
seu.dat@meta.data$cell_type[seu.dat$seurat_clusters==5] <- "Activated cells"

DimPlot(seu.dat, group.by = "cell_type")

qsave(seu.dat, file=file.path(DATA, "harmony_integrated_slingshot_3rd_no_ciliated.qs"), shuffle_control = 0)
