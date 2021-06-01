
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
library(RColorBrewer)

####colors####

sftpc_col <- "#306F1D"
scgb1a1_col <- "#CFE46D"
hopx_col <- "#F090D4"

normal_col <- "#BFDEF5"
tumor_col <- "#8891E6"

tdTom_col = "#e84118"
GFP_col = "#4cd137"



TAM_col <- "#BFDEF5"
Cas9_col <- "#87C8C2"
wk2_col <- "#E6E460"
wk4_col <- "#EDA541"
ep_col <- "#8891E6"


###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="V:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="N:/"} else {
    PATH_Y="/C010-datasets/"}
PATH_ICGC <- "/icgc/"


DATA = paste0(PATH, "Reka/33_CoO_lung/scRNASeq/data")
RESULT = "./output/"
CODE = paste0(DATA, "/../code/")
CALLS = paste0(PATH_Y, "External/Sotillo_mouse_lung_cell/scRNASeq_CoO_lung/cellranger_output/")

source(paste0(CODE, "rasterize_Seurat.R"))
design <- image_read("C:/Users/tothr/Dropbox/Manuscript/Manuscript/Biorender Figures/Figure 5A.png")
design <- ggplot() +
  background_image(design)

seu.dat <- qread(file=file.path(DATA, "harmony_integrated_slingshot_3rd_no_ciliated.qs"))

######supplementary files, marker lists 

#seu.dat@misc$marker_genes$by_cluster$cluster <- recode(as.character(seu.dat@misc$marker_genes$by_cluster$cluster), 
#                                                       "2" = "1", "0" = "2", "7"="6", "6"="7",
#                                                       "11"="8", "10"="9", "9"="10", "13"="11", "8"="12", "12"="13",
#                                                       "1"="14")

#seu.dat@misc$marker_genes$by_cluster$cluster <- factor(seu.dat@misc$marker_genes$by_cluster$cluster, levels=1:14)
#marker_list <- seu.dat@misc$marker_genes$by_cluster %>%
#group_split(cluster)

#for (i in 1:length(marker_list)){
#  write.xlsx(as.data.frame(marker_list[[i]]), row.names = F, col.names=T, 
#             file="Supplementary_tables/marker_genes_clusters.xlsx", sheetName=as.character(marker_list[[i]]$cluster[1]), append=T)
#}



#rm(marker_list)

#seu.dat@misc$marker_genes$by_sample$sample <- recode(as.character(seu.dat@misc$marker_genes$by_sample$sample), 
#                                                       "control01_cas9" = "Cas9", "control01_tam" = "TAM", 
#                                                       "eav1_2wks_1" = "2 wk-1", "eav1_2wks_2" = "2 wk-2",
#                                                       "eav1_4wks_1" = "4 wk-1", "eav1_4wks_2" = "4 wk-GFP",
#                                                       "eav1_endpoint" = "Tumour")

#seu.dat@misc$marker_genes$by_sample$sample <- factor(seu.dat@misc$marker_genes$by_sample$sample, 
#                                                     levels=c("Cas9", "TAM", "2 wk-1", "2 wk-2", "4 wk-1", "4 wk-GFP", "Tumour"))
#marker_list <- seu.dat@misc$marker_genes$by_sample %>%
 # group_split(sample, )

#for (i in 1:length(marker_list)){
#  write.xlsx(as.data.frame(marker_list[[i]]), row.names = F, col.names=T, 
#             file="Supplementary_tables/marker_genes_samples.xlsx", sheetName=as.character(marker_list[[i]]$sample[1]), append=T)
#}

#rm(marker_list)

#seu.dat <- readRDS(file=file.path(DATA, "harmony_integrated_slingshot_3rd_no_ciliated.RDS"))

seu.dat@meta.data$orig.ident <- recode(seu.dat@meta.data$orig.ident, "control01_cas9" = "Cas9", "control01_tam" = "TAM", 
                                       "eav1_2wks_1" = "2 wk-1", "eav1_2wks_2" = "2 wk-2",
                                       "eav1_4wks_1" = "4 wk-1", "eav1_4wks_2" = "4 wk-GFP",
                                       "eav1_endpoint" = "Tumour")


seu.dat@meta.data$orig.ident  <- factor(seu.dat@meta.data$orig.ident, levels= c("TAM", "Cas9",
                                                                                "2 wk-1", "2 wk-2",
                                                                                "4 wk-1", "4 wk-GFP",  "Tumour"))
#mypal <- c(pal_rickandmorty()(12), pal_futurama()(2))
# "2" = "1", "0" = "2", "7"="6", "6"="7", "11"="8", "10"="9",
#"9"="10", "13"="11", "8"="12", "12"="13", "1"="14")

mypal <- c("2"="#FFFF12", "14"="#C4E538", "1"= "#09EBE3",
          "3"="#FDA7DF", "4"= "#ED4C67", "5" = "#F79F1F", "7"= "#4599ED", 
          "6" = "#1833A1", "12" = "#C94FF5", 
          "10"  = "#B53471", "9"="#EE5A24", "8"= '#009432',"13"= "#01040A", "11"= "#8367E0")
mypal <- mypal[order(as.numeric(names(mypal)))]  
  
mypal_clusters <- c(TAM_col, Cas9_col, alpha(wk2_col, alpha=0.7), wk2_col, alpha(wk4_col, alpha=0.7),wk4_col,  ep_col)

#####design plots ####

seu.dat@meta.data$TAM_control <- ifelse(seu.dat@meta.data$orig.ident=="TAM", 10, 0)
seu.dat@meta.data$Cas9_control <- ifelse(seu.dat@meta.data$orig.ident=="Cas9", 10, 0)
seu.dat@meta.data$two_weeks_timepoint_rep1 <- ifelse(seu.dat@meta.data$orig.ident=="2 wk-1", 10, 0)
seu.dat@meta.data$two_weeks_timepoint_rep2 <- ifelse(seu.dat@meta.data$orig.ident=="2 wk-2", 10, 0)
seu.dat@meta.data$four_weeks_timepoint_rep1 <- ifelse(seu.dat@meta.data$orig.ident=="4 wk-1", 10, 0)
seu.dat@meta.data$four_weeks_timepoint_rep2 <- ifelse(seu.dat@meta.data$orig.ident=="4 wk-GFP", 10, 0)
seu.dat@meta.data$endpoint <- ifelse(seu.dat@meta.data$orig.ident=="Tumour", 10, 0)

tc <- FeaturePlot_rast(seu.dat, features = c("TAM_control", "Cas9_control", 
                                             "two_weeks_timepoint_rep1", "two_weeks_timepoint_rep2",
                                             "four_weeks_timepoint_rep1", "four_weeks_timepoint_rep2",
                                             "endpoint"), ncol = 4, 
            pt.size = 0.05, col=c(alpha("grey", 0.2), alpha(mypal_clusters[1], 0.8)), order=T, rasterize=T, xlab = "TAM", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
  theme(axis.title.y = element_blank(), axis.text = element_text(size=8))

tc_big <- FeaturePlot_rast(seu.dat, features = c("TAM_control"), 
                           pt.size =  0.5, col=c(alpha("lightgrey", 0.6), alpha(mypal_clusters[1], 1)), order=T, rasterize=T, xlab = "TAM", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
  theme(axis.title.y = element_blank(), axis.text = element_text(size=8), plot.background=element_blank())
ggsave(plot=tc_big, filename="scRNASeq/UMAP_TAM.pdf", width = 5.5, height = 4)

c9 <- FeaturePlot_rast(seu.dat, features = c("TAM_control", "Cas9_control", 
                                             "two_weeks_timepoint_rep1", "two_weeks_timepoint_rep2",
                                             "four_weeks_timepoint_rep1", "four_weeks_timepoint_rep2",
                                             "endpoint"), ncol = 4,  
                  pt.size = 0.05, col=c(alpha("grey", 0.2), alpha(mypal_clusters[2], 0.8)), order=T, rasterize=T, xlab="Cas9", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
  theme(axis.title.y = element_blank(), axis.text = element_text(size=8))

ggsave(plot=FeaturePlot_rast(seu.dat, features = c("Cas9_control"), 
                             pt.size =  0.5, col=c(alpha("lightgrey", 0.6), alpha(mypal_clusters[2], 1)), order=T, rasterize=T, xlab = "Cas9", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
         theme(axis.title.y = element_blank(), axis.text = element_text(size=8), plot.background=element_blank()), filename="scRNASeq/UMAP_Cas9.pdf", width = 5.5, height = 4)


tw_1 <- FeaturePlot_rast(seu.dat, features = c("TAM_control", "Cas9_control", 
                                               "two_weeks_timepoint_rep1", "two_weeks_timepoint_rep2",
                                               "four_weeks_timepoint_rep1", "four_weeks_timepoint_rep2",
                                               "endpoint"), ncol = 4,  
                  pt.size = 0.05, col=c(alpha("grey", 0.2), alpha(mypal_clusters[4], 0.8)), order=T, rasterize=T, xlab="2 weeks", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
  theme(axis.title.y = element_blank(), axis.text = element_text(size=8))

ggsave(plot=FeaturePlot_rast(seu.dat, features = c("two_weeks_timepoint_rep1"), 
                             pt.size =  0.5, col=c(alpha("lightgrey", 0.6), alpha(mypal_clusters[4], 1)), order=T, rasterize=T, xlab = "2 weeks", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
         theme(axis.title.y = element_blank(), axis.text = element_text(size=8), plot.background=element_blank()), filename="scRNASeq/UMAP_2wk1.pdf", width = 5.5, height = 4)

ggsave(plot=FeaturePlot_rast(seu.dat, features = c("two_weeks_timepoint_rep2"), 
                             pt.size =  0.5, col=c(alpha("lightgrey", 0.6), alpha(mypal_clusters[4], 1)), order=T, rasterize=T, xlab = "2 weeks", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
         theme(axis.title.y = element_blank(), axis.text = element_text(size=8), plot.background=element_blank()), filename="scRNASeq/UMAP_2wk2.pdf", width = 5.5, height = 4)



fw_1 <- FeaturePlot_rast(seu.dat, features = c("TAM_control", "Cas9_control", 
                                               "two_weeks_timepoint_rep1", "two_weeks_timepoint_rep2",
                                               "four_weeks_timepoint_rep1", "four_weeks_timepoint_rep2",
                                               "endpoint"), ncol = 4, 
                  pt.size = 0.05, col=c(alpha("grey", 0.2), alpha(mypal_clusters[6], 0.8)), order=T, rasterize=T, xlab="4 weeks", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
  theme(axis.title.y = element_blank(), axis.text = element_text(size=8))

ggsave(plot=FeaturePlot_rast(seu.dat, features = c("four_weeks_timepoint_rep1"), 
                             pt.size =  0.5, col=c(alpha("lightgrey", 0.6), alpha(mypal_clusters[6], 1)), order=T, rasterize=T, xlab = "4 weeks", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
         theme(axis.title.y = element_blank(), axis.text = element_text(size=8), plot.background=element_blank()), filename="scRNASeq/UMAP_4wk1.pdf", width = 5.5, height = 4)

ggsave(plot=FeaturePlot_rast(seu.dat, features = c("four_weeks_timepoint_rep2"), 
                             pt.size =  0.5, col=c(alpha("lightgrey", 0.6), alpha(mypal_clusters[6], 1)), order=T, rasterize=T, xlab = "4 weeks", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
         theme(axis.title.y = element_blank(), axis.text = element_text(size=8), plot.background=element_blank()), filename="scRNASeq/UMAP_4wk2.pdf", width = 5.5, height = 4)


ep <- FeaturePlot_rast(seu.dat, features = c("TAM_control", "Cas9_control", 
                                             "two_weeks_timepoint_rep1", "two_weeks_timepoint_rep2",
                                             "four_weeks_timepoint_rep1", "four_weeks_timepoint_rep2",
                                             "endpoint"), ncol = 4, 
                  pt.size = 0.05, col=c(alpha("grey", 0.2), alpha(mypal_clusters[7], 0.8)), order=T, rasterize=T, xlab="Tumour", xlab.size=10) & NoLegend() + theme( plot.title = element_blank()) + 
  theme(axis.title.y  = element_blank(), axis.text.x  = element_text(size=6))


ggsave(plot=FeaturePlot_rast(seu.dat, features = c("endpoint"), 
                             pt.size = 0.5, col=c(alpha("lightgrey", 0.6), alpha(mypal_clusters[7], 1)), order=T, rasterize=T, xlab = "Tumour", xlab.size=10) & NoLegend() + theme( plot.title = element_blank() ) +
         theme(axis.title.y = element_blank(), axis.text = element_text(size=8), plot.background=element_blank()), filename="scRNASeq/UMAP_tumour.pdf", width = 5.5, height = 4)



                                                              
p_1 <- wrap_plots(list(tc[[1]],c9[[2]],tw_1[[3]],fw_1[[5]],ep[[7]], 
                       ggplot()+ theme_void(), ggplot()+ theme_void(), tw_1[[4]], fw_1[[6]]), ncol=5)

#p1 <- DimPlot(seu.dat, reduction = "umap", label = F, pt.size = 0.1, group.by = "orig.ident", cols = alpha(mypal[c(11,4,5,6,8)], 0.4))


colnames(seu.dat@meta.data)[colnames(seu.dat@meta.data)=="log_ALK_norm"] <- "ALK\n translocation"

######marker genes ####
p2 <- FeaturePlot_rast(seu.dat, features = c("Sftpc", "Scgb1a1", "Scgb3a2", "Hopx", "Foxj1", "ALK\n translocation"), ncol = 3, 
                  pt.size = 0.1, col=c(alpha("lightgrey", 0.2), "#24325FFF"), rasterize = T)& NoAxes() & NoLegend()+
  theme(plot.title = element_text( face = "italic"))
#
#p2 <- CombinePlots(plots = p2)
#p2#p2_1 <- FeaturePlot(seu.dat, group_by("log_alk"), col=c("grey", "#24325FFF"))& NoAxes() & NoLegend()


######clusters#####


seu.dat$seurat_clusters <- recode(as.character(seu.dat$seurat_clusters), "2" = "1", "0" = "2", "7"="6", "6"="7", "11"="8", "10"="9", "9"="10", "13"="11", "8"="12", "12"="13", "1"="14")

seu.dat$seurat_clusters <- factor(seu.dat$seurat_clusters, levels=c(1:14))

seu.dat <- SetIdent(seu.dat, value="seurat_clusters")
                                  
p3 <- DimPlot_rast(seu.dat, reduction = "umap", label = TRUE, pt.size = 0.1, cols=alpha(mypal, 0.4), label.size = 4, rasterize=TRUE, group.by="seurat_clusters") +
  theme(axis.title = element_text(size=8))+ guides(color = guide_legend(override.aes = list(size=3, alpha=1)))

scheme <- DimPlot_rast(seu.dat, reduction = "umap", label = TRUE, pt.size = 0.4, cols=rep("grey", 14), label.size = 6, rasterize=TRUE, group.by="seurat_clusters") +
     theme(axis.title = element_text(size=8))+ theme(legend.position = "none")

ggsave(plot=scheme, filename="scRNASeq/UMAP_scheme.pdf", width = 4, height = 3)
## cell types

#p4 <- DimPlot_rast(seu.dat, reduction = "umap", label = F, group.by = "pruned.labels", pt.size = 0.1, cols=alpha(mypal, 0.4), rasterize=TRUE) + 
#  theme(legend.position="right", legend.text = element_text(size=8), legend.key.size = unit(3, "mm"), axis.title = element_text(size=8))
#+
#  guides(col=guide_legend(ncol =4, override.aes = list(size=2)))

#####cluster composition####
plot_data <- seu.dat@meta.data[,c("seurat_clusters", "orig.ident")]
plot_data2 <- plot_data %>% 
  dplyr::count(orig.ident, seurat_clusters) %>%
group_by(orig.ident) %>% 
  mutate(sums_i=sum(n)) %>% 
  group_by(seurat_clusters)%>%
mutate(sums=sum(n), norm.n=(1000/sums_i)*n) %>% 
  mutate(norm.sums=sum(norm.n)) %>% 
  mutate(norm.perc=(norm.n/norm.sums)*100)

plot_data3 <- plot_data %>% 
  dplyr::count(orig.ident, seurat_clusters) %>%
  group_by(orig.ident) %>% 
  mutate(sums_i=sum(n)) %>% 
  group_by(seurat_clusters)%>%
  mutate(sums=sum(n), norm.n=(1000/sums_i)*n) %>% 
  group_by(orig.ident) %>%
  mutate(norm.sums=sum(norm.n)) %>% 
  mutate(norm.perc=(norm.n/norm.sums)*100)


p5 <- ggplot(plot_data2, aes(x=orig.ident, y=norm.perc, fill=seurat_clusters)) +
geom_bar(stat="identity", position = position_fill(), colour="black") + ylab("Percentage")+
scale_fill_manual(values = mypal, breaks = )+theme_bw()+ xlab("Timpeoints")+
guides(fill=guide_legend(title="Clusters"))+
theme(legend.position = "top", panel.grid = element_blank(), legend.title = element_blank(), legend.key.size = unit(2, "mm"))+guides(fill=guide_legend(nrow =2))

p5_2 <- ggplot(plot_data3, aes(x=seurat_clusters, y=norm.perc, fill=orig.ident)) +
  geom_bar(stat="identity", position = position_fill(), colour="black") + ylab("Percentage")+
  scale_fill_manual(values = mypal_clusters, breaks = )+theme_bw()+ xlab("Clusters")+
  guides(fill=guide_legend(title="Timepoints"))+
  theme(legend.position = "top", panel.grid = element_blank(), legend.title = element_blank(), 
        legend.key.size = unit(2, "mm"))+guides(fill=guide_legend(nrow =2), legend.text = element_text(size=10))



###########Modules and scores#############

Krt8_cells<- as.data.frame(read_excel(file.path(DATA, "Krt8_cell_populations_2020_17358_MOESM6_ESM.xlsx"), 
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


#marker_lists[["stress_related"]] <- c("Egr1", "Fos","Zfp36","Fosb","Nr4a1","Nfkbia","Atf3", "Dusp1", "Ppp1r15a")
p6 <- DotPlot(seu.dat, features = marker_lists[1:5], cols =c("white", "#24325FFF"), dot.scale=4) + RotatedAxis()+
  ylab("Seurat clusters") + 
      theme(legend.title = element_text(size=11), legend.text = element_text(size=10), axis.title.x =  element_blank(), 
            axis.text=element_text(size=8), axis.title.y =  element_text(size=8))+ 
  guides(size = guide_legend(title = "Percent \nExpressed"), color = guide_colorbar(title = "Average \nExpression"))

p6_2 <- DotPlot(seu.dat, features = marker_lists[6:length(marker_lists)], cols =c("white", "#24325FFF"), dot.scale=4) + RotatedAxis()+
  ylab("Seurat clusters") + 
  theme(legend.title = element_text(size=11), legend.text = element_text(size=10), axis.title.x =  element_blank(), 
        axis.text=element_text(size=8), axis.title.y =  element_text(size=8))+ 
  guides(size = guide_legend(title = "Percent \nExpressed"), color = guide_colorbar(title = "Average \nExpression"))


names(marker_lists)[names(marker_lists)=="Club \n progenitor"] <- "Club progenitor"
#seu.dat@meta.data <- seu.dat@meta.data[,1:89]
seu.dat <- AddModuleScore(
  object = seu.dat,
  features = marker_lists,
  ctrl = 100,
  name = names(marker_lists)
)

colnames(seu.dat@meta.data)[(ncol(seu.dat@meta.data)-length(marker_lists)+1):ncol(seu.dat@meta.data)] <- 
  paste0(colnames(seu.dat@meta.data)[(ncol(seu.dat@meta.data)-length(marker_lists)+1):ncol(seu.dat@meta.data)], "_signatures")

####module scores, combined####

mean_sig_clusters <- data.frame(module=NA, cluster=NA, mean_z=NA) # define empty obj
for (modules in colnames(seu.dat@meta.data)[grep("signatures", colnames(seu.dat@meta.data))]){
for (seurat_clusters in  levels(seu.dat)) {
  mean_sig_clusters <- rbind(mean_sig_clusters, data.frame(module=modules,cluster=seurat_clusters, mean_z=mean(seu.dat@meta.data[which(Cells(seu.dat) %in%  WhichCells(seu.dat,ident=seurat_clusters)), modules])))
}
}
mean_sig_clusters <- mean_sig_clusters[!is.na(mean_sig_clusters$module),]
mean_sig_clusters$mean_z <- as.numeric(mean_sig_clusters$mean_z)
mean_sig_clusters$cluster <- factor(mean_sig_clusters$cluster, levels=rev(1:14)) 
#ggplot(mean_sig_clusters, aes(module, cluster))+geom_tile(mapping = aes(fill=mean_z))

mean_sig_clusters$module <- mean_sig_clusters$module %>%
  gsub("_signatures","", x = .) %>%
  gsub("(.*)([[:digit:]]{1}$)","\\1", x = .) %>%
  gsub("(.*)(1$)","\\1", x = .) %>%
  gsub(".", " ", x = ., fixed=T) %>%
  gsub("AT$", "AT1", x=.)

mean_sig_clusters$module <-   factor(mean_sig_clusters$module, levels=unique(mean_sig_clusters$module))

p7 <- ggplot(mean_sig_clusters, aes(module, cluster))+geom_tile(mapping = aes(fill=mean_z))+
  scale_fill_gradient2(low="darkblue", high="darkred", mid="white", midpoint=0, name="Module score")+
  ylab("Clusters")+
  xlab("Signatures")+ theme_bw()+
  #guides(fill=guide_legend(title="Module score"))+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

p7_split_1 <-  ggplot(mean_sig_clusters[mean_sig_clusters$module %in% c("AT2", "Club", "AT1","Basal","Goblet"),], aes(module, cluster))+geom_tile(mapping = aes(fill=mean_z))+
  scale_fill_gradient2(low="darkblue", high="darkred", mid="white", midpoint=0, name="Module score")+
  ylab("Clusters")+
  xlab("Signatures")+ theme_bw()+
  #guides(fill=guide_legend(title="Module score"))+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
p7_split_2 <-  ggplot(mean_sig_clusters[!(mean_sig_clusters$module %in% c("AT2", "Club", "AT1","Basal","Goblet", "Tumour")),], aes(module, cluster))+geom_tile(mapping = aes(fill=mean_z))+
  scale_fill_gradient2(low="darkblue", high="darkred", mid="white", midpoint=0, name="Module score")+
  ylab("Clusters")+
  xlab("Signatures")+ theme_bw()+
  #guides(fill=guide_legend(title="Module score"))+
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave(p7_split_1, filename = paste0("scRNASeq/Figure_signatures_1.pdf"), height = 5, width=5)
ggsave(p7_split_2, filename = paste0("scRNASeq/Figure_signatures_2.pdf"), height = 5, width=5)

###Supplemetary files###

#names(marker_lists) <- gsub("/", "_", names(marker_lists), fixed = T)

#for (i in 1:length(marker_lists)){
#write.xlsx(marker_lists[[i]], row.names = F, col.names=F, 
#          file="Supplementary_tables/marker_genes_signatures.xlsx", sheetName=names(marker_lists)[i], append=T)
#}

#####cell types####

#seu.dat$cell_types <- recode(seu.dat@meta.data$seurat_clusters, "0" = "Club cells", "1"="Tumor", "2"= "Club cells", "3" = "Club cells", "4"="Club cells", 
#                             "5" =  "Primed progenitors", "6"  = "AT2 cells", "7"= "Club/AT2 progenitor", "8"= "High plasticity tumor", "9"= "Proliferating",
#                             "10" = "H2-K1 high Club-like cells", "11" = "Stress induced Club cells", "12" ="AT1/AT2 bipotent progenitors", "13" = "Pre-cancer")

name_vec <- rep(NA, length(unique(seu.dat$cell_type)))
names(name_vec) <- levels(factor(seu.dat$cell_type))
name_vec["Tumor"] <- mypal[1]
name_vec["Proliferating cells"] <- mypal[10]
name_vec["Tumor"] <- mypal[2]
name_vec["Club cells"] <- mypal[3]
name_vec["AT2 cells"] <- mypal[7]
name_vec["Activated cells"] <- mypal[6]
name_vec["H2-K1 high club-like progenitors"] <- mypal[11]
name_vec["Pre-tumor stage"] <- mypal[14]
name_vec["AT1/AT2-like tumor"] <- mypal[13]
name_vec["AT2-like progenitors"] <- mypal[8]
name_vec["Club cells"] <- mypal[5]
name_vec["Club-like progenitors"] <- mypal[3]


p8 <- DimPlot_rast(seu.dat, reduction = "umap",label = F, pt.size = 0.1, cols=alpha(name_vec, 0.4), label.size = 4, rasterize=TRUE, group.by="cell_type") +
  theme(axis.title = element_text(size=8), legend.key.size = unit(2, "mm"))+ guides(color = guide_legend(override.aes = list(size=3, alpha=1)))


#####celltype composition####

plot_data <- seu.dat@meta.data[,c("cell_type", "orig.ident")]

plot_data <- plot_data %>% 
  dplyr::count(cell_type, orig.ident) %>%
  group_by(orig.ident) %>% 
  mutate(sums_i=sum(n)) %>% 
  #group_by(cell_types)%>%
  mutate(sums=sum(n), norm.n=(1000/sums_i)*n) %>% 
  mutate(norm.sums=sum(norm.n)) %>% 
  mutate(norm.perc=(norm.n/norm.sums)*100)

p9 <- ggplot(plot_data, aes(x=orig.ident, y=norm.perc, fill=cell_type)) +
  geom_bar(stat="identity", position = position_fill(), colour="black") + ylab("Percentage")+
  scale_fill_manual(values =name_vec, breaks = )+theme_bw()+ xlab("Timpeoints")+
  guides(fill=guide_legend(title="Clusters"))+ 
  theme(legend.position = "top", panel.grid = element_blank(), legend.title = element_blank(),legend.text = element_text(size=10),
        legend.key.size = unit(2, "mm"), axis.text.x = element_text(angle=45,vjust=1, hjust=1))+guides(fill=guide_legend(nrow =4))

#####velocity and PAGA#####

velocity <- image_read("V:/Reka/33_CoO_lung/Figures/dyn_velocity.png")
velocity <- ggplot() +
  background_image(velocity) + theme_classic()

PAGA_velo <- image_read("V:/Reka/33_CoO_lung/Figures/scvelo_PAGA_velo.png")
PAGA_velo <- ggplot() +
  background_image(PAGA_velo)


##final plotting####
#plots to show 
#design
#p_1 -design plots
#p2 - marker genes
#p3 - clusters
#p5_2 - cluster composition
#p7 modules
#p8 cell types
#p9 cell type composition

# 
# layout <- c(
#   area(t = 1, b = 5, l = 1, r = 9),
#   area(t = 6, b = 7, l = 3, r = 9), 
#   area(t = 3, b = 7, l = 12, r = 16),
#   area(t = 1, b = 15, l = 16, r = 17),
#   area(t = 11, b = 15, l = 1, r = 8),
#   area(t = 12, b = 15, l = 9, r = 11), 
#   area(t = 12, b = 15, l = 12, r = 16),
#   area(t = 16, b = 20, l = 1, r = 10),
#   area(t = 17, b = 20, l = 12, r = 15), 
#   area(t = 16, b = 22, l = 16, r = 17)
#  # area(t = 22, b = 23, l = 10, r = 15)
# )
# 
# patch <- p2
# patch <- patch + plot_layout(tag_level = 'new')
# patch_1 <- design  + plot_layout(tag_level = 'new')
# #p <- (design + p_1 +patch+plot_spacer()+p3+p5+p4+p7+p6+p8+p9) + plot_layout(design = layout)
# p <- (design + p_1 +patch+plot_spacer()+p3+p5_2+p7+p8+p9+plot_spacer()) + plot_layout(design = layout)


g1 <- plot_grid(design, p_1, align = "v", ncol=1, rel_widths = c(1.3, 0.8), axis="l")
g2 <- plot_grid(g1, NULL, plot_grid(NULL, p2, NULL, rel_heights=c(0.3, 1, 0.3), ncol=1), NULL, align = "h", rel_widths = c(1,0.2, 0.7, 0.2), axis = "tb", nrow=1)

g3 <- plot_grid(p3, NULL, p5_2, NULL, plot_grid(NULL, p7, NULL, align="v", rel_heights=c(0.1, 1, 0.1), ncol=1), align = "h", rel_widths = c(1.5,0.1, 1,0.1, 0.8), axis="bt", nrow=1)
g4 <- plot_grid(p8,NULL, plot_grid(NULL, p9, NULL, align="v", rel_heights=c(0.1, 1, 0.1), ncol=1), NULL,  align="v", rel_widths = c(3,0.3, 1, 0.7), axis="t", nrow=1)
g5 <- plot_grid(velocity, NULL, PAGA_velo, align = "h", axis="tb", rel_widths  = c(1,0.1,1), nrow=1)

g6 <- plot_grid(g2, g3, g4, g5, align = "v", axis="lr", rel_heights = c(1.3, 1,1, 1), ncol=1)



ggsave(g6, filename = "scRNASeq/Figure_3_v3.tiff", width = 14, height = 18)
ggsave(g6, filename = "scRNASeq/Figure_3_v3.pdf", width = 14, height = 18)

plot_list <- list("p_1"=p_1, "p2"=p2, "p3"= p3,"p5"= p5,"p5_2"= p5_2,"p6"= p6,"p7" = p7, "p8" = p8, "p9"= p9, "PAGA_velo"= PAGA_velo,"velocity"= velocity)

for (plots in names(plot_list)){
  
  ggsave(plot_list[[plots]], filename = paste0("scRNASeq/Figure_", plots, ".tiff"))
  ggsave(plot_list[[plots]], filename = paste0("scRNASeq/Figure_", plots, ".pdf"))
}


#####additional plots####


#### dotplot tumor ####
ggsave( DotPlot(seu.dat, features = marker_lists[["Tumour"]], cols =c("white", "#24325FFF"), dot.scale=4) + RotatedAxis()+
  ylab("Seurat clusters") + 
  theme(legend.title = element_text(size=11), legend.text = element_text(size=10), axis.title.x =  element_blank(), 
        axis.text=element_text(size=8), axis.title.y =  element_text(size=8))+ 
  guides(size = guide_legend(title = "Percent \nExpressed"), color = guide_colorbar(title = "Average \nExpression")), 
  filename= "scRNASeq/dotplot_tumor.pdf", width = 5, height = 3)

####module violins###

for (names in colnames(seu.dat@meta.data)[grep("signatures", colnames(seu.dat@meta.data))]){

ggsave(Seurat::VlnPlot(seu.dat, features = names, cols=mypal, pt.size = 0), filename = paste0("scRNASeq/violin_", names, ".tiff")) 
  ggsave(Seurat::VlnPlot(seu.dat, features = names, cols=mypal, pt.size = 0), filename = paste0("scRNASeq/violin_", names, ".pdf")) 
  
  }

ggsave(Seurat::VlnPlot(seu.dat, features = "ALK_norm", cols=mypal, pt.size = 0), filename = paste0("scRNASeq/violin_", "ALK_norm", ".tiff")) 
ggsave(Seurat::VlnPlot(seu.dat, features = "ALK_norm", cols=mypal, pt.size = 0), filename = paste0("scRNASeq/violin_", "ALK_norm", ".pdf")) 


spectral <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", alpha("#5e4fa2", 0.6))

for (module in colnames(seu.dat@meta.data)[grep("signatures", colnames(seu.dat@meta.data))]){
mod_name <- module %>%
  gsub("_signatures","", x = .) %>%
  gsub("(.*)([[:digit:]]{1}$)","\\1", x = .) %>%
  gsub("(.*)(1$)","\\1", x = .) %>%
  gsub(".", " ", x = ., fixed=T) 




ggsave(FeaturePlot_rast(seu.dat, features = module, order=T, rasterize = T)+
         scale_color_gradientn(colours = c(alpha("lightgrey", 0.2), rev(spectral)))+
         theme( plot.background = element_rect(colour = NA))+
         ggtitle(mod_name), filename = paste0("scRNASeq/feature_", module, ".tiff"), width=8, height = 6)

ggsave(FeaturePlot_rast(seu.dat, features = module, order=T, rasterize = T)+
         scale_color_gradientn(colours = c(alpha("lightgrey", 0.2), rev(spectral)))+
         theme( plot.background = element_rect(colour = NA))+
         ggtitle(mod_name), filename = paste0("scRNASeq/feature_", module, ".pdf"), width=8, height = 6)

}
    
####markers#### 
genes <- c("Nkx2-1", "Lyz2", "Cd274", "Ly6i", "Ly6a", "Iigp1", "Igtp", 
           "Isg15", "Ifitm3", "Sftpc", "Scgb1a1", "Scgb3a2", "Hopx", 
           "Foxj1", "Chia1", "Ager", "AW112010", "H2-K1", "Cd74")
for (gene in genes){
  

  #ggsave(FeaturePlot_rast(seu.dat, features = gene, rasterize=T, 
  #col=c(alpha("lightgrey", 0.2),  alpha("#440154FF", 0.6), alpha("#20A387FF", .6), alpha("#FDE725FF", 0.6)), order=T)+
  #theme(plot.title = element_text( face = "italic"), plot.background = element_rect(colour = NA)), 
  #filename = paste0("scRNASeq/feature_", gene, ".tiff"), width = 8, height = 6)
  ggsave(FeaturePlot_rast(seu.dat, features = gene, rasterize=T, order=T)+
           theme(plot.title = element_text( face = "italic"), plot.background = element_rect(colour = NA))+
           scale_color_gradientn(colours = c(alpha("lightgrey", 0.2), rev(spectral))), 
         filename = paste0("scRNASeq/feature_", gene, ".pdf"), width = 8, height = 6)
  }

#####cell cycle###

ggsave(DimPlot_rast(seu.dat, group.by = "Phase", rasterize = T),
       filename = paste0("scRNASeq/cell_cycle.pdf"), width = 8, height = 6)



######modules#####



modules <- read.delim(file.path(DATA, paste0("cNMF/cNMF.spectra.k_9.dt_0_01.consensus.txt")), row.names = 1)

modules <- t(modules)
colnames(modules) <- paste0("module_", 1:9)
tam_gene_scores <- read.delim("V:/Reka/33_CoO_lung/scRNASeq/data/cNMF/tam_scores.txt", row.names=1)

combined_modules <- merge(tam_gene_scores, modules, by="row.names", all=F)
rownames(combined_modules) <- combined_modules$Row.names
combined_modules <- combined_modules[,-1]
M <- cor(combined_modules)

module_names <- c("Interferon", "Club", "Proliferating", "Stem-like tumour", "Immune activation/inflammation", 
                  "Activated Club", "AT1-like tumour", "Normal AT2", "Regeneration-like tumour")

colnames(M) <- rownames(M)<- c("AT1/AT2", "High cycling", "Hepatic-gastric/inflammation", "Biosynthetic mixed identity", "Stress response", "Highly mixed", 
                 "EMT", "GI epithelium-like", "Embryonic-liver like", "Gastric-like", "AT2-like", 
                 module_names)

M2 <- M[1:11,12:20] 
M2 <- M2[,c("Club", "Activated Club", "Normal AT2", "Proliferating", "Interferon", "Immune activation/inflammation", 
            "Regeneration-like tumour", "AT1-like tumour", "Stem-like tumour")]

module_names_lev <- c("Club", "Activated Club", "Normal AT2", "Proliferating", "Interferon", "Immune activation/inflammation", 
                  "Regeneration-like tumour", "AT1-like tumour", "Stem-like tumour")
cor_plot_data <- M2 %>% as.data.frame %>% 
  tibble::rownames_to_column() %>% 
  gather(key="Modules", value="Correlation", -rowname) %>% 
  mutate(Modules=factor(Modules, levels=module_names_lev)) %>%
  mutate(rowname=factor(rowname, levels=rev(c("AT1/AT2", "High cycling", "Hepatic-gastric/inflammation", "Biosynthetic mixed identity", "Stress response", "Highly mixed", 
                           "EMT", "GI epithelium-like", "Embryonic-liver like", "Gastric-like", "AT2-like")))) 
cor_plot <- cor_plot_data %>%  
  ggplot(aes(rowname, Modules, fill=Correlation))+geom_tile()+xlab("Modules, Marjanovic et al.")+ylab("Expression modules")+
  coord_flip()+theme(panel.background = element_rect(fill=NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_gradient2(low=muted("blue"), mid="white", high=muted("red"))

 #cor_plot <- ggcorrplot(M, method = "square")
 ggsave(cor_plot, filename = paste0("scRNASeq/corrplot_modules.tiff"))
 ggsave(cor_plot, filename = paste0("scRNASeq/corrplot_modules.pdf"))
 
 
 k=9
 modules <- read.delim(file.path(DATA, paste0("cNMF/cNMF.gene_spectra_score.k_", k, ".dt_0_01.txt")), row.names = 1)
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
 
 scaled_re <- scale(re_ordered)
 pdf("scRNASeq/modules_markers.pdf", height = 9, width = 7)
 ht <- Heatmap(scaled_re, cluster_columns = F, cluster_rows=F, col=colorRamp2(breaks = c(-4, 0, 4),
                                                                              colors=c("darkblue", "white", "darkred")),
               row_names_gp = gpar(fontsize=7, angle=45), name="z-score",
               rect_gp = gpar(col = "darkgrey", lwd=0.5),column_names_gp = gpar(fontsize=10), column_names_rot = 45)
 draw(ht, padding = unit(c(1, 2, 0.2, 0.2), "cm"))
 dev.off()

 names(rnames_l) <-  module_names_lev
 
 seu.dat <- AddModuleScore(
   object = seu.dat,
   features = rnames_l,
   ctrl = 100,
   name = names(rnames_l)
 )
 
 colnames(seu.dat@meta.data)[(ncol(seu.dat@meta.data)-length(rnames_l)+1):ncol(seu.dat@meta.data)] <- 
   paste0(colnames(seu.dat@meta.data)[(ncol(seu.dat@meta.data)-length(rnames_l)+1):ncol(seu.dat@meta.data)], "_modules")
 
 for (module in colnames(seu.dat@meta.data)[grep("_modules", colnames(seu.dat@meta.data))]){
   mod_name <- module %>%
     gsub("_modules","", x = .) %>%
     gsub("(.*)([[:digit:]]{1}$)","\\1", x = .) %>%
     gsub("(.*)(1$)","\\1", x = .) %>%
     gsub(".", " ", x = ., fixed=T) 
   

   #ggsave(
     FeaturePlot_rast(seu.dat, features = module, order=T, rasterize = T)+
            scale_color_gradientn(colours = c(alpha("lightgrey", 0.2), rev(spectral)))+
            theme( plot.background = element_rect(colour = NA))+
            ggtitle(mod_name)
  #        ,filename = paste0("scRNASeq/feature_", module, ".pdf"), width=8, height = 6)
   
 }
 
 ####module scores, combined####
 
 mean_sig_clusters <- data.frame(module=NA, cluster=NA, mean_z=NA) # define empty obj
 for (modules in colnames(seu.dat@meta.data)[grep("_modules", colnames(seu.dat@meta.data))]){
   for (seurat_clusters in  levels(seu.dat)) {
     mean_sig_clusters <- rbind(mean_sig_clusters, data.frame(module=modules,cluster=seurat_clusters, mean_z=mean(seu.dat@meta.data[which(Cells(seu.dat) %in%  WhichCells(seu.dat,ident=seurat_clusters)), modules])))
   }
 }
 mean_sig_clusters <- mean_sig_clusters[!is.na(mean_sig_clusters$module),]
 mean_sig_clusters$mean_z <- as.numeric(mean_sig_clusters$mean_z)
 mean_sig_clusters$cluster <- factor(mean_sig_clusters$cluster, levels=rev(1:14)) 
 #ggplot(mean_sig_clusters, aes(module, cluster))+geom_tile(mapping = aes(fill=mean_z))
 
 #mean_sig_clusters$module <- mean_sig_clusters$module %>%
 #  gsub("(.*)([[:digit:]]{1}$)","\\1", x = .) %>%
 #  gsub("(.*)(1$)","\\1", x = .) %>%
 #  gsub(".", " ", x = ., fixed=T) 
 mean_sig_clusters$module <- gsub("_modules", "", mean_sig_clusters$module)
 mean_sig_clusters$module <- gsub("[[:digit:]]{1}$", "", mean_sig_clusters$module)
 mean_sig_clusters$module <- gsub(".", " ", mean_sig_clusters$module, fixed=T)
 mean_sig_clusters$module <- gsub("immune", "Immune", mean_sig_clusters$module, fixed=T)
 mean_sig_clusters$module <- gsub(" like", "-like", mean_sig_clusters$module, fixed=T)
 mean_sig_clusters$module <- gsub("activation ", "activatiion/", mean_sig_clusters$module, fixed=T)
 
 
 mean_sig_clusters$module <-   factor(mean_sig_clusters$module, levels=unique(mean_sig_clusters$module))
 #mean_sig_clusters <- mean_sig_clusters[-1,]
 p7 <- ggplot(mean_sig_clusters, aes(module, cluster))+geom_tile(mapping = aes(fill=mean_z))+
   scale_fill_gradient2(low="darkblue", high="darkred", mid="white", midpoint=0, name="Module score")+
   ylab("Clusters")+
   xlab("Signatures")+ theme_minimal()+
   theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
 ggsave(p7, filename =  paste0("scRNASeq/modules.tiff"), width = 8, height = 8)
 ggsave(p7, filename =  paste0("scRNASeq/modules.pdf"), width = 8, height = 8)
 
 # correlation modules - signatures
 
 
 
# marker_lists[["Combined Krt8 ADI/DATP/PATS"]] <- c(marker_lists[["Krt8 ADI"]], marker_lists[["DATP"]], marker_lists[["PATS"]])
 
# seu.dat <- AddModuleScore(
#   object = seu.dat,
#   features = marker_lists[c("Combined Krt8 ADI/DATP/PATS", "Club")],
#   ctrl = 100,
#   name = c("Combined Krt8 ADI/DATP/PATS", "uninportant")
# )
# colnames(seu.dat@meta.data)[colnames(seu.dat@meta.data)=="Combined.Krt8.ADI.DATP.PATS1"] <- "Combined.Krt8.ADI.DATP.PATS1_signatures"
 
 M <- cor(seu.dat@meta.data[grep("modules|signatures", colnames(seu.dat@meta.data))])
 M2 <- M[grep("signatures", rownames(M)),grep("modules", colnames(M))]
 rownames(M2) <- rownames(M2) %>% 
  gsub("_signatures", "", .) %>%
  gsub("[[:digit:]]{1}$", "", .) %>%
  gsub(".", " ", ., fixed=T) %>%
  gsub("AT1 AT2", "AT1/AT2", .,  fixed=T) %>% 
   gsub("ADI DATP PATS", "ADI/DATP/PATS", .)
  
  
 rownames(M2)[!grepl("AT2|AT1", rownames(M2))] <-   gsub("[[:digit:]]{1}$", "", rownames(M2)[!grepl("AT2|AT1", rownames(M2))])
 
 colnames(M2) <- colnames(M2) %>% 
   gsub("_modules", "", .) %>%
   gsub("[[:digit:]]{1}$", "", .) %>%
   gsub(".", " ", .,  fixed=T) %>%
   gsub("AT1 AT2", "AT1/AT2", .,  fixed=T) %>%
   gsub("Immune activation inflammation", "Immune activation/inflammation", .,  fixed=T) %>%
   gsub(" like", "-like", .)
 
 M2 <- M2[!(rownames(M2) %in% c("Mixed AT1/AT2", "Tumour")),]
 
 cor_plot_data2 <- M2 %>% as.data.frame %>% 
   tibble::rownames_to_column() %>% 
   gather(key="Modules", value="Correlation", -rowname) %>% 
   mutate(Modules=factor(Modules, levels=rev(colnames(M2)))) %>%
   mutate(rowname=factor(rowname, levels=rev(rownames(M2)))) 
 cor_plot2 <- cor_plot_data2 %>%
   ggplot(aes(rowname, factor(Modules, levels = rev(levels(Modules))), fill=Correlation))+geom_tile()+xlab("Published signatures")+ylab("Expression modules")+
   coord_flip()+theme(panel.background = element_rect(fill=NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
   scale_fill_gradient2(low=muted("blue"), mid="white", high=muted("red"))
 
 
 ggsave(cor_plot2, filename =  paste0("scRNASeq/corplot_modules_signatures.tiff"), width = 8, height = 8)
 ggsave(cor_plot2, filename =  paste0("scRNASeq/corplot_modules_signatures.pdf"), width = 8, height = 8)
 
 
 cor_plot3 <- cor_plot_data2 %>%
   filter(rowname %in% c("AT2", "Club", "AT1",  "Basal", "Goblet")) %>%
   ggplot(aes(rowname, factor(Modules, levels = rev(levels(Modules))), fill=Correlation))+geom_tile()+xlab("Published signatures")+ylab("Expression modules")+
   coord_flip()+theme(panel.background = element_rect(fill=NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
   scale_fill_gradient2(low=muted("blue"), mid="white", high=muted("red"))
 ggsave(cor_plot3, filename =  paste0("scRNASeq/corplot_modules_signatures_sel.pdf"), width = 8, height = 8)
 
 p2 <- cor_plot_data2 %>%
   filter(rowname %in% c("Krt8 ADI", "DATP",  "PATS")) %>%
   ggplot(aes(rowname, factor(Modules, levels = rev(levels(Modules))), fill=Correlation))+geom_tile()+xlab("Published signatures")+ylab("Expression modules")+
   coord_flip()+theme(panel.background = element_rect(fill=NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
   scale_fill_gradient2(low=muted("blue"), mid="white", high=muted("red"), limits=c(-1, 1))#+theme(aspect.ratio=4/11)
 
 
 cor_plot <- cor_plot_data %>%  
   filter(rowname %in% c("AT1/AT2", "Highly mixed",  "Embryonic-liver like")) %>%
   ggplot(aes(rowname, Modules, fill=Correlation))+geom_tile()+xlab("Modules, Marjanovic et al.")+ylab("Expression modules")+
   coord_flip()+theme(panel.background = element_rect(fill=NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
   scale_fill_gradient2(low=muted("blue"), mid="white", high=muted("red"))
 
 p1 <- cor_plot+scale_fill_gradient2(low=muted("blue"), mid="white", high=muted("red"), limits=c(-1, 1))+
   theme(legend.position = "none", axis.title.x=element_blank(),axis.text.x=element_blank())#+theme(aspect.ratio=1)

 
 p1/p2 + plot_layout(guides = "collect")
 
 ggsave( filename =  paste0("scRNASeq/corplot_modules_signatures_combined_selected.pdf"), 
         width = 7, height = 5.5)
 
 
 
 
 
 
 ### enrichment ###
 k=9
 

 modules <- read.delim(file.path(DATA, paste0("cNMF/cNMF.gene_spectra_score.k_", k, ".dt_0_01.txt")), row.names = 1)
 modules <- t(modules)
 modules <- modules[-1,]
 colnames(modules) <- module_names
 
 library(org.Mm.eg.db)
 suppressPackageStartupMessages(library(clusterProfiler))
i <- "Interferon signaling"
 
    gene.df <- bitr(rownames(modules)[order(modules[,i], decreasing = T)][1:200], fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Mm.eg.db)

   
   enr <- enrichGO(gene          = gene.df$ENTREZID,
                             OrgDb         = org.Mm.eg.db,
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05,
                             readable      = TRUE)
   enr@result$module <- i
   p1 <- dotplot(enr, showCategory=30) + ggtitle(paste0("GO dotplot for module ", i))+
     scale_color_gradient(low  = "darkblue", high="grey")
   ggsave(p1, filename = "scRNASeq/enrichment_interferon.pdf", height = 7, width = 10)
   
   i <- "Immune activation/inflammation"
   
   gene.df <- bitr(rownames(modules)[order(modules[,i], decreasing = T)][1:200], fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Mm.eg.db)
   
   
   enr <- enrichGO(gene          = gene.df$ENTREZID,
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
   enr@result$module <- i
   #enr@result$p.adjust <- -log10(enr@result$p.adjust)
   p1 <- dotplot(enr, showCategory=30) + ggtitle(paste0("GO dotplot for module ", i))+
     scale_color_gradient(low  = "darkblue", high="grey")
   ggsave(p1, filename = "scRNASeq/enrichment_immune.pdf", height = 7, width = 10)
   
  
   i <- "Tumour"
   
   gene.df <- bitr(rownames(modules)[order(modules[,i], decreasing = T)][1:200], fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Mm.eg.db)
   
   
   enr <- enrichKEGG(gene         = gene.df$ENTREZID,
                     organism     = 'mouse',
                     pvalueCutoff = 0.05)
   
   enr@result$module <- i
   
   
   wpgmtfile <- paste0(DATA, "/wikipathways-20210310-gmt-Mus_musculus.gmt")
   wp2gene <- read.gmt(wpgmtfile)
   wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
   wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
   wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
   
   ewp <- enricher(gene.df$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
   ewp@result$Count <- as.numeric(gsub("(.)/[[:digit:]]+", "\\1", ewp@result$GeneRatio))
   ewp@result$module <- i
   enr@result <- rbind(enr@result, ewp@result)
   
   p1 <- dotplot(enr, showCategory=30) + ggtitle(paste0("KEGG and Wikipathways dotplot for module ", i))+
     scale_color_gradient(low  = "darkblue", high="grey")
   ggsave(p1, filename = "scRNASeq/enrichment_tumour.pdf", height = 7, width = 10)
   
   i <- "Regeneration-like state"
   
   gene.df <- bitr(rownames(modules)[order(modules[,i], decreasing = T)][1:200], fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Mm.eg.db)
   
   
   enr <- enrichKEGG(gene         = gene.df$ENTREZID,
                     organism     = 'mouse',
                     pvalueCutoff = 0.05)
   
   enr@result$module <- i
   
   
   wpgmtfile <- paste0(DATA, "/wikipathways-20210310-gmt-Mus_musculus.gmt")
   wp2gene <- read.gmt(wpgmtfile)
   wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
   wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
   wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
   
   ewp <- enricher(gene.df$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
   ewp@result$Count <- as.numeric(gsub("(.)/[[:digit:]]+", "\\1", ewp@result$GeneRatio))
   ewp@result$module <- i
   enr@result <- rbind(enr@result, ewp@result)
   
   p1 <- dotplot(enr, showCategory=30) + ggtitle(paste0("KEGG and Wikipathways dotplot for ", i))+
     scale_color_gradient(low  = "darkblue", high="grey")
   ggsave(p1, filename = "scRNASeq/enrichment_regeneration.pdf", height = 7, width = 10)
   