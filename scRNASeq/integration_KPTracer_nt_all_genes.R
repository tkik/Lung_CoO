#Convert("data/marjanovics/adata_processed.combined.h5ad", dest = "h5seurat", overwrite = TRUE, assay = "RNA")
#KPTumors <- LoadH5Seurat("data/marjanovics/adata_processed.combined.h5seurat", assays = "RNA")


#saveRDS(KPTumors, file="P:/3_Lung_CoO/scRNASeq/data/marjanovics/KPTumors.RDS")


###integration###
library(scRNAseq)
library(zellkonverter)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggalluvial)
#library(harmony)
library(cerebroApp)
library(biomaRt)
library(future)
library(cerebroApp)
library(qs)
library(readxl)
library(ggpubr)
library(stringr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(SingleR)
library(harmony)
library(ggsci)
library(scExtras)
library(SingleCellExperiment)
library(ggimage)
library(magick)
library(ggpubr)
library(slingshot)
library(ggrastr)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(qs)
library(corrplot)
library(ggcorrplot)
library(tidyr)
library(scales)
library(xlsx)
library(stringr)
library(RColorBrewer)
# Reading in the individual datasets



mypal <- c("2"="#FFFF12", "14"="#C4E538", "1"= "#09EBE3", "0"= "#25CED1", "15"= "#FF8A5B", "16"="#EA526F", "17"="#904C77", "18"="#ECB8A5",
           "19"="#294D4A", "20"="#776472",
           "3"="#FDA7DF", "4"= "#ED4C67", "5" = "#F79F1F", "7"= "#4599ED",
           "6" = "#1833A1", "12" = "#C94FF5",
           "10"  = "#B53471", "9"="#EE5A24", "8"= '#009432',"13"= "#01040A", "11"= "#8367E0")
mypal <- mypal[order(as.numeric(names(mypal)))]

mypal2 <- mypal
names(mypal2) <- NULL


sce <- readH5AD("c:/Temp/marjanovics/adata_processed.nt.h5ad")
sce_combined <- readH5AD("c:/Temp/marjanovics/adata_processed.combined.h5ad")


#sce <- sce[rownames(sce_combined),]
intersct <- intersect(colnames(sce), colnames(sce_combined))
sce_combined <- sce_combined[,intersct]
sce <- sce[,intersct]

counts(sce) <- assay(sce, "X")

logcounts(sce) <- assay(sce, "counts")

seu.dat <- as.Seurat(sce)
# change the counts data with the normalized one
#seu.list <- SplitObject(seu.dat, split.by = "Batch_Library")
#seu.list <- lapply(X = seu.list, FUN = function(x) {
#   x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 4000)
#})
seu.dat <- NormalizeData(seu.dat, normalization.method = "LogNormalize", scale.factor = 10000)
seu.dat@assays$originalexp@data <- assay(sce, "X")
seu.dat <- FindVariableFeatures(seu.dat, selection.method = "vst", nfeatures = 4000)


original_data <- qread("c:/Temp/marjanovics/harmony_integrated_slingshot_3rd_no_ciliated.qs")
original_data$seurat_clusters_original <- original_data$seurat_clusters
lung.list <- list(original_data, seu.dat)



features <- SelectIntegrationFeatures(object.list = lung.list, nfeatures = 4000)
lung.list <- lapply(X = lung.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
lung.anchors <- FindIntegrationAnchors(object.list = lung.list, anchor.features = features, reduction = "rpca")
lung.integrated <- IntegrateData(anchorset = lung.anchors, dims = 1:30)


# lung.anchors <- FindIntegrationAnchors(object.list = lung.list, dims = 1:30, scale = TRUE, normalization.method = "LogNormalize", anchor.features=3000)
# lung.integrated <- IntegrateData(anchorset = lung.anchors, dims = 1:30)
DefaultAssay(lung.integrated) <- "integrated"
lung.integrated <- ScaleData(lung.integrated, verbose = FALSE)
lung.integrated <- RunPCA(lung.integrated, npcs = 30, verbose = TRUE)
#lung.integrated <- RunPCA(lung.integrated, npcs = 30, verbose = FALSE, assay = "integrated", features = lung.integrated@assays$integrated@var.features)
lung.integrated <- lung.integrated %>%
    RunUMAP(reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    # FindClusters(resolution = 0.5) %>%
    FindClusters(resolution = 0.5) %>%
    identity()

#saveRDS(lung.integrated, file="data/lung.integrated_rpca_KPTracer.RDS")
qsave(lung.integrated, file="c:/Temp/marjanovics/lung.integrated_rpca_KPTracer_nt_orig.RDS")

lung.integrated <- qread(file="c:/Temp/marjanovics/lung.integrated_rpca_KPTracer_nt_orig.RDS")


lung.integrated@meta.data <- lung.integrated@meta.data %>% tidyr::unite("integrtated_name", c("Cluster.Name","cell_type"), na.rm = TRUE, remove = FALSE)
library(pals)
mypal2 <- c(mypal2, polychrome(10))
mypal3 <-   polychrome(n = 30)
names(mypal2) <- NULL


DimPlot(lung.integrated, reduction = "umap", label = TRUE, pt.size = 0.1, cols=alpha(mypal2, 0.7), label.size = 4,  group.by="integrtated_name") +
    theme(axis.title = element_text(size=8))+ guides(color = guide_legend(override.aes = list(size=3, alpha=1)))


#lung.integrated@meta.data <- lung.integrated@meta.data[,c(1:13, 102, 127)]

lung.integrated$seurat_clusters_original <- recode(as.character(lung.integrated$seurat_clusters_original),
                                                   "2" = "1", "0" = "2", "7"="6", "6"="7", "11"="8", "10"="9", "9"="10", "13"="11", "8"="12", "12"="13", "1"="14")


lung.integrated$seurat_clusters_original_old <- lung.integrated$seurat_clusters_original

lung.integrated$seurat_clusters_original[is.na(lung.integrated$ALK_norm)] <- lung.integrated$Cluster.Name[is.na(lung.integrated$ALK_norm)]

lung.integrated$seurat_clusters_original <- factor(lung.integrated$seurat_clusters_original, levels=c(1:14, unique(lung.integrated$Cluster.Name)))

lung.integrated@meta.data$orig.ident <- recode(lung.integrated@meta.data$orig.ident, "control01_cas9" = "Cas9", "control01_tam" = "TAM",
                                               "eav1_2wks_1" = "2 wk-1", "eav1_2wks_2" = "2 wk-2",
                                               "eav1_4wks_1" = "4 wk-1", "eav1_4wks_2" = "4 wk-GFP",
                                               "eav1_endpoint" = "Tumour")
lung.integrated@meta.data$orig.ident <- factor(lung.integrated@meta.data$orig.ident, levels=rev(levels(as.factor(lung.integrated@meta.data$orig.ident))))


DimPlot(lung.integrated, reduction = "umap", label = TRUE, pt.size = 0.6, cols=alpha(mypal2, 0.2), label.size = 4,  group.by="seurat_clusters_original") +
    theme(axis.title = element_text(size=8))+ guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+ggtitle("Seurat clusters")
#ggsave(filename = "output/integrated_dimplot_clusters_KPTracer_nt.pdf", width=12, height=6)

DimPlot(lung.integrated, reduction = "umap", label = TRUE, pt.size = 0.6, cols=alpha(mypal2, 0.2), label.size = 4,  group.by="seurat_clusters_original_old") +
    theme(axis.title = element_text(size=8))+ guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+ggtitle("Seurat clusters")
#ggsave(filename = "output/integrated_dimplot_clusters_KPTracer_nt_colored_our.pdf", width=12, height=6)

DimPlot(lung.integrated, reduction = "umap", label = TRUE, pt.size = 0.6, cols=alpha(mypal2, 0.2), label.size = 4,  group.by="Cluster.Name") +
    theme(axis.title = element_text(size=8))+ guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+ggtitle("Seurat clusters")
#ggsave(filename = "output/integrated_dimplot_clusters_KPTracer_nt_colored_new.pdf", width=12, height=6)

DimPlot(lung.integrated, reduction = "umap", label = TRUE, pt.size = 0.6, cols=alpha(mypal2, 0.2), label.size = 4,  group.by="seurat_clusters_original") +
    theme(axis.title = element_text(size=8))+ guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+ggtitle("Seurat clusters")
#ggsave(filename = "output/integrated_dimplot_clusters_KPTracer_nt.pdf", width=12, height=6)


DimPlot(lung.integrated, reduction = "umap", label = TRUE, pt.size = 0.7, cols=alpha(mypal2, 0.2), label.size = 4,  group.by="integrtated_name") +
    theme(axis.title = element_text(size=8))+ guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+ggtitle("Original sample")+
    scale_color_manual(values = alpha(mypal2, 0.2))
#ggsave(filename = "output/integrated_dimplot_KPTracer_nt.pdf", width=12, height=6)


plot_data <- lung.integrated@meta.data
plot_data2 <- plot_data %>%
    dplyr::count(integrtated_name, seurat_clusters) %>%
    group_by(integrtated_name) %>%
    mutate(sums_i=sum(n)) %>%
    group_by(seurat_clusters)%>%
    mutate(sums=sum(n), norm.n=(1000/sums_i)*n) %>%
    group_by(integrtated_name) %>%
    mutate(norm.sums=sum(norm.n)) %>%
    mutate(norm.perc=(norm.n/norm.sums)*100)
p5 <- ggplot(plot_data2, aes(x=integrtated_name, y=norm.perc, fill=seurat_clusters)) +
    geom_bar(stat="identity", position = position_fill(), colour="black") + ylab("Percentage")+
    scale_fill_manual(values = mypal2, breaks = )+theme_bw()+ xlab("Timpeoints")+
    guides(fill=guide_legend(title="Clusters"))+
    theme(legend.position = "top", panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1),
          legend.title = element_blank(), legend.key.size = unit(2, "mm"))+guides(fill=guide_legend(nrow =2))
p5

#ggsave(p5,filename = "output/integrated_barplot_KPTracer_nt.pdf", width=12, height=6)


plot_data2 <- plot_data %>%
    dplyr::count(integrtated_name, seurat_clusters) %>%
    group_by(seurat_clusters) %>%
    mutate(sums_i=sum(n)) %>%
    group_by(integrtated_name)%>%
    mutate(sums=sum(n), norm.n=(1000/sums_i)*n) %>%
    group_by(seurat_clusters) %>%
    mutate(norm.sums=sum(norm.n)) %>%
    mutate(norm.perc=(norm.n/norm.sums)*100)
p5 <- ggplot(plot_data2, aes(x=seurat_clusters, y=norm.perc, fill=integrtated_name)) +
    geom_bar(stat="identity", position = position_fill(), colour="black") + ylab("Percentage")+
    scale_fill_manual(values = mypal2, breaks = )+theme_bw()+ xlab("Clusters")+
    guides(fill=guide_legend(title="Original clusters"))+
    theme(legend.position = "top", panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1),
          legend.title = element_blank(), legend.key.size = unit(2, "mm"))+guides(fill=guide_legend(nrow =2))
p5
#ggsave(p5,filename = "output/integrated_barplot_KPTracer_nt_by_cluster_names.pdf", width=12, height=6)

plot_data2 <- plot_data %>%
    dplyr::count(seurat_clusters_original, seurat_clusters) %>%
    group_by(seurat_clusters) %>%
    mutate(sums_i=sum(n)) %>%
    group_by(seurat_clusters_original)%>%
    mutate(sums=sum(n), norm.n=(1000/sums_i)*n) %>%
    group_by(seurat_clusters) %>%
    mutate(norm.sums=sum(norm.n)) %>%
    mutate(norm.perc=(norm.n/norm.sums)*100)
p5 <- ggplot(plot_data2, aes(x=seurat_clusters, y=norm.perc, fill=seurat_clusters_original)) +
    geom_bar(stat="identity", position = position_fill(), colour="black") + ylab("Percentage")+
    scale_fill_manual(values = mypal2, breaks = )+theme_bw()+ xlab("Clusters")+
    guides(fill=guide_legend(title="Original clusters"))+
    theme(legend.position = "top", panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1),
          legend.title = element_blank(), legend.key.size = unit(2, "mm"))+guides(fill=guide_legend(nrow =2))
p5
#ggsave(p5,filename = "output/integrated_barplot_KPTracer_nt_by_cluster_numbers.pdf", width=12, height=6)


plot_data <- lung.integrated@meta.data
plot_data2 <- plot_data %>%
    dplyr::count(seurat_clusters_original, seurat_clusters) %>%
    group_by(seurat_clusters_original) %>%
    mutate(sums_i=sum(n)) %>%
    group_by(seurat_clusters)%>%
    mutate(sums=sum(n), norm.n=(1000/sums_i)*n) %>%
    group_by(seurat_clusters_original) %>%
    mutate(norm.sums=sum(norm.n)) %>%
    mutate(norm.perc=(norm.n/norm.sums)*100)
p5 <- ggplot(plot_data2, aes(x=seurat_clusters_original, y=norm.perc, fill=seurat_clusters)) +
    geom_bar(stat="identity", position = position_fill(), colour="black") + ylab("Percentage")+
    scale_fill_manual(values = mypal2, breaks = )+theme_bw()+ xlab("Timpeoints")+
    guides(fill=guide_legend(title="Clusters"))+
    theme(legend.position = "top", panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1),
          legend.title = element_blank(), legend.key.size = unit(2, "mm"))+guides(fill=guide_legend(nrow =2))
p5
#ggsave(p5,filename = "output/integrated_barplot_KPTracer_nt_2.pdf", width=12, height=6)
#
#markers <- FindAllMarkers(seu.dat)

plot_data2 <- plot_data %>%
    dplyr::count(seurat_clusters_original, seurat_clusters) %>%
    group_by(seurat_clusters_original) %>%
    ggplot(aes(y = n, axis1 = seurat_clusters_original, axis2 = seurat_clusters)) +
    geom_alluvium(aes(fill = seurat_clusters_original), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Original clusters", "New clusters"), expand = c(.05, .05)) +
    scale_fill_manual(values = mypal2) +
    ggtitle("Integration")


p3 <- plot_data %>%
    dplyr::count(seurat_clusters_original, seurat_clusters) %>%
    group_by(seurat_clusters_original) %>%
    ggplot(aes(y = n, axis1 = seurat_clusters_original, axis2 = seurat_clusters)) +
    geom_alluvium(aes(fill = seurat_clusters_original), width = 1/12) +
    geom_stratum(width = 1/12, fill = "black", color = "grey") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Original clusters", "New clusters"), expand = c(.05, .05)) +
    scale_fill_manual(values = mypal2) +
    ggtitle("Integration")
#ggsave(p3,filename = "output/alluvial_plot_KPTracer_nt.pdf", width=12, height=18)


Krt8_cells<- as.data.frame(read_excel(file.path("P:/3_Lung_CoO/scRNASeq/data/", "41467_2020_17358_MOESM6_ESM.xlsx"),
                                      sheet = "cell_types_2", col_types = c("skip",
                                                                            "text", "numeric", "numeric", "numeric",
                                                                            "numeric", "text")))


marker_lists <- lapply(unique(Krt8_cells$cluster),function(x) Krt8_cells[Krt8_cells$cluster==x,"gene"][1:30])
names(marker_lists) <- gsub("(^[a-z]{1})", "\\U\\1", unique(Krt8_cells$cluster), perl = TRUE)
names(marker_lists) <- gsub(" cells", "", names(marker_lists))
marker_lists <- marker_lists[c("AT2", "Club", "AT1","Basal","Goblet", "Krt8 ADI", "Proliferation")]

marker_lists[["Club \n progenitor"]] <- c("Cd74", "Cd14", "H2-K1", "AW112010")

marker_lists[["Mixed AT1/AT2"]] <- c("Clic3", "Ager", "Col4a3", "Dpysl2", "Mmp11", "Akap5", "Rtkn2", "Hopx", "Radil", "Sec14l3")
#marker_lists[["Highly mixed"]] <- c("Epcam", "Cldn7", "Cldn4", "Cars", "Hmces", "Dusp10", "Slc4a11Lif", "Tigit", "Lgals3")
#marker_lists[["EMT program"]] <- c("Lamb1", "Fbln2", "Fxyd5", "Tnnt2", "Inhba", "Dbn1", "Cpe", "Axl", "Igfbp4","Vim")
marker_lists[["Activated Club"]] <- c("Cst3", "H2-Ab1", "H2-Eb1", "Psap", "Ctss", "H2-Aa", "Laptm5", "Actb", "Alox5ap",
                                      "Srgn", "Ccl17", "Lsp1", "Scgb3a1", "S100a4")
marker_lists[["Tumour"]] <-  c("Meg3", "Areg", "Pf4", "Kng2", "Id2", "Serpine1", "Rian", "Dmkn", "Tspan8","Napsa",
                               "Cald1","Ly6c1","Tmem213","Basp1","Thbs1","Igfbp7","Kcnc3","Ly6c2","Arg2","Malt1","Stbd1")
marker_lists[["DATP"]] <- c("Mif","Pold4","Mboat1","AW112010","Hif1a","Pdk4","Cxx1b","Lrrc16","Cdkn2a",
                            "Cdkn1a","Mdm2","Ccnd1","Gdf15","Trp53","Bax","Ifngr1","Ly6a","Ifr7","Cxcl16","Timp1")

marker_lists[["PATS"]] <- c("S100a6","Sfn","Tmsb10","Cldn4","Clu","AW112010","Krt19","Krt18","Anxa1","Krt8",
                            "Tpm2","Krt7","Serpinb9","Anxa2","Ifitm3","Epcam","Lgals3","Tnip3","Sox4","Anxa5",
                            "Cyr61","Cavin3","Tuba1a","S100a11","Serpinb1a","Crip1","Lgals1","Ccl20","Actn1","Mfge8","Ubd",
                            "Tubb5","Ywhah","S100a10","Hspb1","2200002D01Rik","Nfkbia","Cyba","F3","Cd24a","Cd81","Lurap1l",
                            "Myl12a","S100a14","Tpm1","Igfbp7","Marcksl1","Cfl1","Myl12b","Hsp90ab1")
#marker_lists[["Biosynthetic mixed identity"]] <- NULL
marker_lists[["Hepatic-gastric inflammation"]] <- c("Ly6c1Hp","Hp",  "Ttc36", "Ly6c2","Ly6k","Rps4l","Etfbkmt","Atp6v0a4", "Kng2","Orm1" )
marker_lists[["Highly mixed"]] <- c("Epcam","Cldn7","Cldn4","Cars","Hmces","Dusp10","Slc4a11", "Lif", "Tigit", "Lgals3")
marker_lists[["Early tumor markers"]] <- c("Lyz2", "Lamp3", "Sftpa1", "Sftpc", "Cd74", "Cxcl15", "B2m")
marker_lists[["Fate Cluster 1 markers"]] <- c("Gkn2", "Gdf15", "Yy1", "Flrt3", "Vim", "Fabp5", "Marcks", "Acta2", "Emp2")
marker_lists[["Fate Cluster 2 markers"]] <- c("Scgb1a1", "Pax9", "Prom1", "Sox2", "Itgb4", "Cd24a", "Gsto1", "Palc8", "Cldn4", "Tff1")
names(marker_lists)[names(marker_lists)=="Club \n progenitor"] <- "Club progenitor"
#seu.dat@meta.data <- seu.dat@meta.data[,1:89]

module_names <- c("Interferon", "Club", "Proliferating", "Stem-like tumour", "Immune activation/inflammation",
                  "Activated Club", "AT1-like tumour", "Normal AT2", "Regeneration-like tumour")
module_names_lev <- c("Club", "Activated Club", "Normal AT2", "Proliferating", "Interferon", "Immune activation/inflammation",
                      "Regeneration-like tumour", "AT1-like tumour", "Stem-like tumour")

k=9
modules <- read.delim( paste0("data/cNMF.gene_spectra_score.k_", k, ".dt_0_01.txt"), row.names = 1)
modules <- t(modules)
modules <- modules[-1,]
colnames(modules) <- module_names
modules <- modules[,module_names_lev]

re_ordered <- data.frame(matrix(NA, nrow=k*10, ncol=k))
rnames <- list()
for (i in 1:k){
    from <- (i*10-(k))
    to <- i*10
    re_ordered[from:to, ] <- modules[order(modules[,i], decreasing = T)[1:10],]
    rnames[[i]] <-  rownames(modules[order(modules[,i], decreasing = T)[1:10],])
    rnames[[i]] <- gsub(".", "-", rnames[[i]], fixed=T)
}

rnames_l <- rnames
rnames <- unlist(rnames)
rnames[duplicated(rnames)]<- paste0(rnames[duplicated(rnames)], "_1")
rnames[duplicated(rnames)]<- paste0(rnames[duplicated(rnames)], "_2")
rownames(re_ordered) <- rnames
colnames(re_ordered) <- module_names_lev


names(rnames_l) <-  module_names_lev

lung.integrated <- AddModuleScore(
    object = lung.integrated,
    features = rnames_l,
    ctrl = 50,
    name = names(rnames_l)
)


colnames(lung.integrated@meta.data)[(ncol(lung.integrated@meta.data)-length(rnames_l)+1):ncol(lung.integrated@meta.data)] <-
    paste0(str_sub(colnames(lung.integrated@meta.data)[(ncol(lung.integrated@meta.data)-length(rnames_l)+1):ncol(lung.integrated@meta.data)], 1, -2), "_module")

# for (module in colnames(lung.integrated@meta.data)[grep("_module", colnames(lung.integrated@meta.data))]){
#   mod_name <- module %>%
#     gsub("_modules","", x = .) %>%
#     gsub("(.*)([[:digit:]]{1}$)","\\1", x = .) %>%
#     gsub("(.*)(1$)","\\1", x = .) %>%
#     gsub(".", " ", x = ., fixed=T)
#
# }
lung.integrated <- AddModuleScore(
    object = lung.integrated,
    features = marker_lists,
    ctrl = 60,
    name = names(marker_lists)
)
a <- 1:length(colnames(lung.integrated@meta.data)[(ncol(lung.integrated@meta.data)-length(marker_lists)+1):ncol(lung.integrated@meta.data)])
b <- rep("", length(colnames(lung.integrated@meta.data)[(ncol(lung.integrated@meta.data)-length(marker_lists)+1):ncol(lung.integrated@meta.data)]))

colnames(lung.integrated@meta.data)[(ncol(lung.integrated@meta.data)-length(marker_lists)+1):ncol(lung.integrated@meta.data)] <-
    paste0(mapply(gsub, a, b, colnames(lung.integrated@meta.data)[(ncol(lung.integrated@meta.data)-length(marker_lists)+1):ncol(lung.integrated@meta.data)], USE.NAMES = FALSE), "_signatures")
#colnames(lung.integrated@meta.data)[which(colnames(lung.integrated@meta.data)=="seurat_clusters")] <- "seurat_clusters_0.6"
#lung.integrated <- lung.integrated %>%
#RunUMAP(reduction = "pca", dims = 1:30, reduction.name="umap_0.8") %>%
#  FindNeighbors(reduction = "pca", dims = 1:30) %>%
# FindClusters(resolution = 0.5) %>%
#  FindClusters(resolution = 0.9) %>%
#  identity()



#marker_lists[["stress_related"]] <- c("Egr1", "Fos","Zfp36","Fosb","Nr4a1","Nfkbia","Atf3", "Dusp1", "Ppp1r15a")
#p6 <- DotPlot(lung.integrated, features = marker_lists[1:5], cols =c("white", "#24325FFF"), dot.scale=4, ) + RotatedAxis()+
#  ylab("Seurat clusters") +
#  theme(legend.title = element_text(size=11), legend.text = element_text(size=10), axis.title.x =  element_blank(),
#        axis.text=element_text(size=8), axis.title.y =  element_text(size=8))+
#  g
source("P:/3_Lung_CoO/scRNASeq/code/getMarkerGenes.R")
spectral <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", alpha("#5e4fa2", 0.6))
plot_data <- lung.integrated@meta.data
plot_data <- cbind(plot_data, lung.integrated@reductions$umap@cell.embeddings)
plot_data <- plot_data %>%
    filter(!is.na(seurat_clusters_original))
order_factors <- function(fact_var, colname="Club2"){
    fact_var <- factor(plot_data[,fact_var], levels=plot_data %>%
                           group_by(!!sym(fact_var)) %>%
                           summarize(mean_club=mean(!!sym(colname))) %>%
                           arrange(mean_club) %>%
                           pull(!!sym(fact_var)))
}
selected_cols <- colnames(plot_data)[grep("module|signatures", colnames(plot_data))]
plots <- list()
for (cols in selected_cols){
    plot_data$seurat_clusters_original <- order_factors(fact_var = "seurat_clusters_original",  colname=cols)
    plots[[cols]] <- ggplot(plot_data)+geom_boxplot(aes_string(x="seurat_clusters_original", y=cols), fill=3)+coord_flip()+theme_bw()
}

plots <- list()
for (cols in selected_cols){
    plot_data$original_clusters <- order_factors(fact_var = "original_clusters",  colname=cols)
    plots[[cols]] <- ggplot(plot_data)+geom_boxplot(aes_string(x="original_clusters", y=cols), fill=3)+coord_flip()+theme_bw()+
        theme(text = element_text(size = 16))
}
ggarrange(plotlist = plots, ncol=2, nrow = 2)
#ggexport(plotlist = plots, ncol=2, nrow = 2, filename = "output/rebuttal/boxplots_modules.pdf")
#ggsave(plots[[7]], filename="output/rebuttal/boxplots_regeneration.pdf")
#ggsave(plots[[8]], filename="output/rebuttal/boxplots_AT1_AT2.pdf")
#ggsave(plots[[9]], filename="output/rebuttal/boxplots_stem.pdf")
plots <- list()
for (cols in selected_cols){
    plot_data$seurat_clusters_original <- order_factors(fact_var = "seurat_clusters_original",  colname=cols)
    plots[[cols]] <- ggplot(plot_data)+geom_boxplot(aes_string(x="seurat_clusters_original", y=cols), fill=3)+coord_flip()+theme_bw()
}
ggarrange(plotlist = plots, ncol=2, nrow = 2)
#ggexport(plotlist = plots, ncol=2, nrow = 2, filename = "output/boxplots_modules_KPTracer_numbers.pdf", height = 8, width = 12)
plots <- list()
for (cols in selected_cols){
    plot_data$integrtated_name <- order_factors(fact_var = "integrtated_name",  colname=cols)
    plots[[cols]] <- ggplot(plot_data)+geom_boxplot(aes_string(x="integrtated_name", y=cols), fill=3)+coord_flip()+theme_bw()
}
ggarrange(plotlist = plots, ncol=2, nrow = 2)
#ggexport(plotlist = plots, ncol=2, nrow = 2, filename = "output/boxplots_modules_KPTracer.pdf", height = 8, width = 12)
#ggsave(plots[[7]], filename="output/boxplots_regeneration_KPTracer.pdf")
#ggsave(plots[[8]], filename="output/boxplots_AT1_AT2_KPTracer.pdf")
unique(plot_data$seurat_clusters_original)
unique(plot_data$Cluster.Name)
plot_data_t <- plot_data %>%
    summarise_each(funs(t.test(.[seurat_clusters_original == "12"], .[seurat_clusters_original == "Gastric-like"])$p.value),
                   vars = Club_module:Fate.Cluster.2.markers_signatures)
plot_data_t <- list()
plot_data_t[["14 vs Late gastric, Stem.like.tumour_module"]] <- plot_data %>%
    summarise_each(funs(t.test(.[seurat_clusters_original == "14"], .[seurat_clusters_original == "Late gastric"])$p.value),
                   vars = Stem.like.tumour_module)
plot_data_t[["14 vs Late gastric, Stem.like.tumour_module"]] <- plot_data %>%
    summarise_each(funs(t.test(.[seurat_clusters_original == "14"], .[seurat_clusters_original == "Late Gastric"])$p.value),
                   vars = Stem.like.tumour_module)
plot_data_t <- list()
plot_data_t[["12 vs Gastric-like, Regeneration.like.tumour_module"]] <- plot_data %>%
    summarise_each(funs(t.test(.[seurat_clusters_original == "12"], .[seurat_clusters_original == "Gastric-like"])$p.value),
                   vars = Regeneration.like.tumour_module)
plot_data_t[["12 vs Late gastric, Regeneration.like.tumour_module"]] <- plot_data %>%
    summarise_each(funs(t.test(.[seurat_clusters_original == "12"], .[seurat_clusters_original == "Late Gastric"])$p.value),
                   vars = Regeneration.like.tumour_module)
plot_data_t[["13 vs AT1-like, AT1.like.tumour_module"]] <- plot_data %>%
    summarise_each(funs(t.test(.[seurat_clusters_original == "13"], .[seurat_clusters_original == "AT1-like"])$p.value),
                   vars = AT1.like.tumour_module)
plot_data_t[["14 vs Early gastric, Stem.like.tumour_module"]] <- plot_data %>%
    summarise_each(funs(t.test(.[seurat_clusters_original == "14"], .[seurat_clusters_original == "Early gastric"])$p.value),
                   vars = Stem.like.tumour_module)


for (module in colnames(lung.integrated@meta.data)[grep("_module", colnames(lung.integrated@meta.data))]){
    mod_name <- module %>%
        gsub("_module","", x = .) %>%
        gsub("(.*)([[:digit:]]{1}$)","\\1", x = .) %>%
        gsub("(.*)(1$)","\\1", x = .) %>%
        gsub(".", " ", x = ., fixed=T)
    print(FeaturePlot(lung.integrated, features = module,  order=T, raster = T)+
              theme(plot.title = element_text( face = "italic"), plot.background = element_rect(colour = NA))+
              scale_color_gradientn(colours = c(alpha("lightgrey", 0.2), rev(spectral))))
}

lung.integrated$group <- ifelse(lung.integrated$Cluster.Name %in% c( "Gastric−like",  "Early gastric", "AT2-like", "Late Gastric"), "Yang", ifelse(lung.integrated$seurat_clusters_original=="12", "own", NA))
m4_9 <- FindMarkers(lung.integrated, ident.1 = "own", ident.2 = "Yang", group.by = "group")
genes <- rownames(m4_9)
f_plotlist <- list()
for (i in genes[1:20]){
    f_plotlist[[i]] <-  FeaturePlot(lung.integrated, features = i,  order=T, raster = T)+
        theme(plot.title = element_text( face = "italic"), plot.background = element_rect(colour = NA))+
        scale_color_gradientn(colours = c(alpha("lightgrey", 0.2), rev(spectral)))
}
ggarrange(plotlist = rev(f_plotlist), ncol=2, nrow = 2)

lung.integrated$group <- ifelse(lung.integrated$Cluster.Name %in% c( "Gastric−like",  "Early gastric", "AT2-like", "Late Gastric"), "Yang", ifelse(lung.integrated$seurat_clusters_original=="14", "own", NA))

m4_9 <- FindMarkers(lung.integrated, ident.1 = "own", ident.2 = "Yang", group.by = "group")
genes <- rownames(m4_9)
f_plotlist <- list()
for (i in genes[1:20]){
    f_plotlist[[i]] <-  FeaturePlot(lung.integrated, features = i,  order=T, raster = T)+
        theme(plot.title = element_text( face = "italic"), plot.background = element_rect(colour = NA))+
        scale_color_gradientn(colours = c(alpha("lightgrey", 0.2), rev(spectral)))
}
ggarrange(plotlist = rev(f_plotlist), ncol=2, nrow = 2)



