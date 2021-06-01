###Figure 2, methylation


library(Gviz)
library(pheatmap)
library(matrixStats)
library(ComplexHeatmap)
library(knitr)
library(kableExtra)
library(tidyverse)
library(circlize)
library(ggplotify)
library(patchwork)
library(ggrepel)
library(ggsci)
library(ggfortify)
library(MeDeCom)
library(cowplot)
library(org.Mm.eg.db)
library(fgsea)
library(methrix)
####colors####

sftpc_col <- "#306F1D"
scgb1a1_col <- "#CFE46D"
hopx_col <- "#F090D4"

normal_col <- "#BFDEF5"
tumor_col <- "#8891E6"

tdTom_col = "#e84118"
GFP_col = "#4cd137"


AT2_col <- "#096F51"
club_col <- "#CBE557"
AT1_col <- "#C44EF0"
ciliated_col <- "#FF9300"
  

###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="V:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="N:/"} else {
    PATH_Y="/C010-datasets/"}

DATA = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/data/norfp/")
DATA_meth = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/data/")
CODE = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/code/")
RESULT = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/output/")
CALLS = paste0(PATH_Y, "External/2018-10-Sotillo/data/methylDackel/")


########

gseaplot2 <- function(geneSetID, title = "", color = "green",
                      base_size = 11, rel_heights = c(1.5, 0.5, 1), subplots = 1:3,
                      pvalue_table = FALSE, ES_geom = "line", ranks, geneSet, pvaluetable) {
  
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  gsdata2 <- do.call(rbind, lapply(geneSetID, function(id) gsInfo2(geneSetID=id, geneset = geneSet, res = ranks)))
  p <- ggplot(gsdata2, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand = c(0, 0))
  
  es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description),
                        size = 0.8)
  
  p.res <- p + es_layer + theme(legend.position = "top", legend.box="vertical", legend.margin=margin(),
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  
  p.res <- p.res + ylab("Enrichment Score") +
    theme(axis.title.y = element_text(size=8),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), plot.margin = margin(t = 0.2,
                                                              r = 0.2, b = 0, l = 0.2, unit = "cm"))
  p.res <- p.res + geom_hline(yintercept=0, linetype = "dashed", color="grey")
  
  i <- 0
  
  for (term in unique(gsdata2$Description)) {
    idx <- which(gsdata2$ymin != 0 & gsdata2$Description ==
                   term)
    gsdata2[idx, "ymin"] <- i
    gsdata2[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  
  #browser()
  
  p2 <- ggplot(gsdata2, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin,
                                                            ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) +
    theme_classic(base_size) + theme(legend.position = "none",
                                     plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
                                     axis.ticks = element_blank(), axis.text = element_blank(),
                                     axis.line.x = element_blank()) + scale_x_continuous(expand = c(0,
                                                                                                    0)) + scale_y_continuous(expand = c(0, 0))
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x,
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Correlation Coef.") + xlab("Rank in Ordered Dataset") +
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2,
                               l = 0.2, unit = "cm"), axis.title.y = element_text(size=8))+scale_x_continuous(expand = c(0, 0))
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(labels= paste0(c(pvaluetable[geneSetID,"pathway"]$pathway), ", p value= ", round(c(pvaluetable[geneSetID , "pval"])$pval, digits = 3)),
                                        values = color, guide=guide_legend(nrow=2))
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- pvaluetable[geneSetID , c("pathway", "pval",
                                    "padj")]
    pd <- pd[order(pd[, 1], decreasing = FALSE), ]
    rownames(pd) <- pd$pathway
    pd <- pd[, -1]
    pd <- round(pd, 4)
    tp <- enrichplot:::tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp, xmin = quantile(p.res$data$x,
                                            0.5), xmax = quantile(p.res$data$x, 0.95), ymin = quantile(p.res$data$runningScore,
                                                                                                       0.75), ymax = quantile(p.res$data$runningScore,
                                                                                                                              0.9))
  }
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(),
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  plot_grid(plotlist = plotlist, ncol = 1, align = "v",
            rel_heights = rel_heights)
  
  
}



gsInfo2 <-  function (geneset=marker_lists_short, geneSetID, res=lmc1_rank)
{
  geneList <- res
  if (is.numeric(geneSetID))
    geneSetID <- names(geneset)[geneSetID]
  geneSet <- geneset[[geneSetID]]
  exponent <- 1
  df <- DOSE:::gseaScores(res, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- res
  df$Description <- geneSetID
  return(df)
}

######
#### PCA plot####

mat  <- as.data.frame(readRDS(paste0(DATA, "medecom/bivalent_enhancers.RDS")))

rowsds <- rowSds(as.matrix(mat[,-(1:3)]), na.rm=T)
mat <- mat[order(rowsds, decreasing = T)[1:100000],]
mat <- mat[complete.cases(mat[,-(1:3)]),]

#mat <- mat[,colnames(mat)!="MsHOPX_control01"]

anno <- data.frame(row.names=colnames(mat)[-(1:3)], 
                   Construct=ifelse(grepl("CCSP", colnames(mat)[-(1:3)]), "Scgb1a1", ifelse(grepl("SPC", colnames(mat)[-(1:3)]), "Sftpc", "Hopx")), 
                   "Sample type" = ifelse(grepl("control", colnames(mat)[-(1:3)]), "Normal", "Tumour"))


df <- data.frame(colors=c( scgb1a1_col, sftpc_col, hopx_col), name=c( "Scgb1a1", "Sftpc", "Hopx"))



pca_res <- prcomp(t(as.matrix(mat[,-(1:3)])))
pca_plot <- autoplot(pca_res, data = anno, colour="Construct", shape="Sample.type", size=8)+theme_bw()+
  theme(legend.text = element_text(size=10))+
  scale_color_manual(name="Line", values = c(hopx_col, scgb1a1_col, sftpc_col))+ 
  stat_ellipse(aes(fill = Sample.type), geom = "polygon", type="norm")+
  scale_fill_manual(name="N/T", values = c(alpha(normal_col, 0.3), alpha(tumor_col, 0.3)))+labs(shape="N/T")



#### B. MeDeCom plot ####

medecom.result2 <-  readRDS( file = paste0(DATA, "medecom/poised_enhancers_new_samples.RDS"))

K_sel <- 4
lambda_sel <- 0.001
proportions <- MeDeCom::getProportions(medecom.result2, K=K_sel, lambda=lambda_sel)
colnames(proportions) <- colnames(mat[,-(1:3)])
rownames(proportions) <- paste0("LMC", c(1,2,4,3)) 


LMCs <- MeDeCom::getLMCs(medecom.result2, K=K_sel, lambda=lambda_sel)
colnames(LMCs) <- paste0("LMC", 1:K_sel)
#LMCs <- cbind(mat[,1:3], LMCs)


annotation <- data.frame(staining=gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\1", colnames(mat[,-(1:3)])),
                         type= gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\2", colnames(mat[,-(1:3)])),
                         origin = ifelse(grepl("rfp", colnames(mat[,-(1:3)])), "RFP", "GFP"))
rownames(annotation) <- colnames(mat[,-(1:3)])


p <- pheatmap(proportions,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12, annotation_col = annotation)

new_groups <- list(SPC=colnames(proportions)[proportions["LMC3",]>0.1 & grepl("tumor", colnames(proportions))],
                   CCSP=colnames(proportions)[proportions["LMC1",]>0.1 & grepl("tumor", colnames(proportions))])


proportions_orig <- proportions 
proportions <- proportions[rownames(proportions)!="LMC4",]

anno <- data.frame(row.names=colnames(proportions), 
                 Construct=ifelse(grepl("CCSP", colnames(proportions)), "Scgb1a1", ifelse(grepl("SPC", colnames(proportions)), "Sftpc", "Hopx")), 
                 "Sample type" = ifelse(grepl("control", colnames(proportions)), "Normal", "Tumour"))
anno$Sample.type <- factor(anno$Sample.type)

#mycols <- colorRamp2(breaks=c(0,0.5,1), colors=c("#24325FFF", "#fdfde6", "#e84118"))
mycols <- colorRamp2(breaks=c(0,0.5,1), colors=c("#002157", "#f5f5f5",  "#BC3C29FF"))

# Normal SPC
 ha1 <- HeatmapAnnotation("Construct" = anno$Construct[which(anno$Construct=="Sftpc" & anno$Sample.type=="Normal")],
                              "Sample type" = anno$Sample.type[anno$Construct=="Sftpc" & anno$Sample.type=="Normal"],
                              col = list("Construct"= c(Sftpc = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
                                         "Sample type"=c(Normal = normal_col,  "Tumour" = tumor_col)), 
                              annotation_legend_param = list(Construct = list(nrow = 1), 
                                                             "Sample type" = list(nrow = 1)), 
                          annotation_name_gp = gpar(col="white"))
 ht1 <- Heatmap(proportions[,rownames(anno)[anno$Construct=="Sftpc" & anno$Sample.type=="Normal"]], name = "Proportions", cluster_rows = FALSE,   
               show_row_names = F, 
               top_annotation = ha1, show_column_names = F,
               col = mycols, heatmap_legend_param = list(direction = "horizontal"), column_title = "Sftpc",  column_title_gp = gpar(fontsize=10))

# Normal CCSP
 
 ha2 <- HeatmapAnnotation("Construct" = anno$Construct[which(anno$Construct=="Scgb1a1" & anno$Sample.type=="Normal")],
                          "Sample type" = anno$Sample.type[anno$Construct=="Scgb1a1" & anno$Sample.type=="Normal"],
                          col = list("Construct"= c(Sftpc = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
                                     "Sample type"=c(Normal = normal_col,  "Tumour" = tumor_col)), 
                          annotation_legend_param = list(Construct = list(nrow = 1), 
                                                         "Sample type" = list(nrow = 1)),  
                          annotation_name_gp = gpar(col="white"))
 ht2 <- Heatmap(proportions[,rownames(anno)[anno$Construct=="Scgb1a1" & anno$Sample.type=="Normal"]], name = "Proportions", cluster_rows = FALSE,   
                show_row_names = F, 
                top_annotation = ha2, show_column_names = F,
                col = mycols, heatmap_legend_param = list(direction = "horizontal"), column_title = "Scgb1a1",  column_title_gp = gpar(fontsize=10))


# Normal HOPX
 
 ha3 <- HeatmapAnnotation("Construct" = anno$Construct[which(anno$Construct=="Hopx" & anno$Sample.type=="Normal" & rownames(anno)!="MsHOPX_control01")],
                          "Sample type" = anno$Sample.type[anno$Construct=="Hopx" & anno$Sample.type=="Normal" & rownames(anno)!="MsHOPX_control01"],
                          col = list("Construct"= c(Sftpc = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
                                     "Sample type"=c(Normal = normal_col,  "Tumour" = tumor_col)),
                          annotation_legend_param = list(Construct = list(nrow = 1), 
                                                         "Sample type" = list(nrow = 1)), 
                          annotation_name_gp = gpar(col="white"))
 ht3 <- Heatmap(proportions[,rownames(anno)[anno$Construct=="Hopx" & anno$Sample.type=="Normal" & rownames(anno)!="MsHOPX_control01"]], name = "Proportions", cluster_rows = FALSE,   
                show_row_names = F, 
                top_annotation = ha3, show_column_names = F,
                col = mycols, heatmap_legend_param = list(direction = "horizontal"), column_title = "Hopx",  column_title_gp = gpar(fontsize=10))
 
 # Green tumors with SPC origin
 
 ha4 <- HeatmapAnnotation("Construct" = anno$Construct[rownames(anno) %in% new_groups$SPC],
                          "Sample type" = anno$Sample.type[rownames(anno) %in% new_groups$SPC],
                          col = list("Construct"= c(Sftpc = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
                                     "Sample type"=c(Normal = normal_col,  "Tumour" = tumor_col)), 
                          annotation_legend_param = list(Construct = list(nrow = 1), 
                                                         "Sample type" = list(nrow = 1)), 
                          annotation_name_gp = gpar(col="white"), show_legend = F)
 
 ht4 <- Heatmap(proportions[,rownames(anno)[rownames(anno) %in% new_groups$SPC]], name = "Proportions", cluster_rows = FALSE,   
                show_row_names = TRUE, 
                top_annotation = ha4, show_column_names = F,
                col = mycols, heatmap_legend_param = list(direction = "horizontal"), column_title = "AT2 origin",  
                column_title_gp = gpar(fontsize=10))


# Green tumors with SPC origin 
 
 ha5 <- HeatmapAnnotation("Construct" = anno$Construct[rownames(anno) %in% new_groups$CCSP[!grepl("rfp", new_groups$CCSP)]],
                          "Sample type" = anno$Sample.type[rownames(anno) %in% new_groups$CCSP[!grepl("rfp", new_groups$CCSP)]],
                          col = list("Construct"= c(Sftpc = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
                                     "Sample type"=c(Normal = normal_col,  "Tumour" = tumor_col)),
                          annotation_legend_param = list(Construct = list(nrow = 1), 
                                                         "Sample type" = list(nrow = 1)), 
                          annotation_name_gp = gpar(col="white"), show_legend = F)
 ht5 <- Heatmap(proportions[,rownames(anno)[rownames(anno) %in% new_groups$CCSP[!grepl("rfp", new_groups$CCSP)]]], name = "Proportions", cluster_rows = FALSE,   
                show_row_names = TRUE, 
                top_annotation = ha5, show_column_names = F,
                col = mycols, heatmap_legend_param = list(direction = "horizontal"), column_title = "Club origin",  column_title_gp = gpar(fontsize=10))
 
 ###### Red  tumors with SPC or CCSP origin #####
 # 
 # ha6 <- HeatmapAnnotation("Construct" = anno$Construct[rownames(anno) %in% c(new_groups$CCSP[grepl("rfp", new_groups$CCSP)],"MsSPC_tumor03.rfp")],
 #                          "Sample type" = anno$Sample.type[rownames(anno) %in% c(new_groups$CCSP[grepl("rfp", new_groups$CCSP)],"MsSPC_tumor03.rfp")],
 #                          "Color" = anno$Color[rownames(anno) %in% c(new_groups$CCSP[grepl("rfp", new_groups$CCSP)],"MsSPC_tumor03.rfp")],
 #                                                    col = list("Construct"= c(Sftpc = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
 #                                                                 "Sample type"=c(Normal = normal_col,  "Tumour" = tumor_col), 
 #                                                                 "Color"=c("tdTomato" = tdTom_col, GFP = GFP_col)),
 #                          annotation_legend_param = list(Construct = list(nrow = 1), 
 #                                                         "Sample type" = list(nrow = 1), 
 #                                                         "Color" = list(nrow=2)), 
 #                          annotation_name_gp = gpar(col="white"))
 # 
 # 
 # ht6 <- Heatmap(proportions[,rownames(anno) %in% c(new_groups$CCSP[grepl("rfp", new_groups$CCSP)],"MsSPC_tumor03.rfp")], name = "Proportions", cluster_rows = FALSE,   
 #                show_row_names = TRUE, 
 #                top_annotation = ha6, show_column_names = F,
 #                col = mycols, heatmap_legend_param = list(direction = "horizontal"), column_title = "tdTomato tumors,\nClub and AT2",  column_title_gp = gpar(fontsize=8))
 # 
 # Red  tumors with unknown origin
 # 
 # ha7 <- HeatmapAnnotation("Construct" = anno$Construct[rownames(anno) %in% new_groups$HOPX],
 #                          "Sample type" = anno$Sample.type[rownames(anno) %in% new_groups$HOPX],
 #                          "Color" = anno$Color[rownames(anno) %in% new_groups$HOPX],
 #                          col = list("Construct"= c(Sftpc = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
 #                                     "Sample type"=c(Normal = normal_col,  "Tumour" = tumor_col), 
 #                                     "Color"=c("tdTomato" = tdTom_col, GFP = GFP_col)),
 #                          annotation_legend_param = list(Construct = list(nrow = 1), 
 #                                                         "Sample type" = list(nrow = 1), 
 #                                                         "Color" = list(nrow=2)), 
 #                          annotation_name_gp = gpar(col="white"))
 # 
 # 
 # ht7 <- Heatmap(proportions[,rownames(anno) %in% new_groups$HOPX], name = "Proportions", cluster_rows = FALSE,   
 #                show_row_names = TRUE, 
 #                top_annotation = ha7, show_column_names = F,
 #                col = mycols, heatmap_legend_param = list(direction = "horizontal"), column_title = "tdTomato tumours, \nunknown", column_title_gp = gpar(fontsize=8))
 # 
 # 
 
 ###all
 #proportions_orig <- proportions_orig[c(1,2,4,3),]
 ha_all <- HeatmapAnnotation("Construct" = anno$Construct,
                          "Sample type" = anno$Sample.type,
                          col = list("Construct"= c(Sftpc = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
                                     "Sample type"=c(Normal = normal_col,  "Tumour" = tumor_col)), 
                          annotation_legend_param = list(Construct = list(nrow = 1), 
                                                         "Sample type" = list(nrow = 1)), 
                          annotation_name_gp = gpar(col="white"))
 ht_all <- Heatmap(proportions_orig, name = "Proportions", cluster_rows = FALSE,   
                show_row_names = T, 
                top_annotation = ha_all, show_column_names = F,
                col = mycols, heatmap_legend_param = list(direction = "horizontal"),  column_title_gp = gpar(fontsize=10))
 
 
 ht_list_all = ggplotify::as.ggplot(grid.grabExpr(draw(ht_all , column_title = "Methylation deconvolution", merge_legend = F, 
                                                   heatmap_legend_side = "bottom", annotation_legend_side = "bottom"), wrap = T))
 
 ggsave(ht_list_all, filename = "methylation/deconvolution_all.pdf", height = 5, width = 7)
ht_all_3 <- Heatmap(proportions_orig[c("LMC1", "LMC2", "LMC3"),], name = "Proportions", cluster_rows = FALSE,   
                   show_row_names = T, 
                   top_annotation = ha_all, show_column_names = F,
                   col = mycols, heatmap_legend_param = list(direction = "horizontal"),  column_title_gp = gpar(fontsize=10))
ht_list_all_3 = ggplotify::as.ggplot(grid.grabExpr(draw(ht_all_3, column_title = "Methylation deconvolution", merge_legend = F, 
                                                      heatmap_legend_side = "bottom", annotation_legend_side = "bottom"), wrap = T))
ggsave(ht_list_all_3, filename = "methylation/deconvolution_all_v2.pdf", height = 5, width = 7)

 
 ht_list <- (ht1+ht2+ht3)
 ht_list2 <- (ht4+ht5)
 
# ht_list2 <- (ht4+ht5+ht6+ht7)
 
 
 

#ht_list <- draw(ht, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
ht_list = ggplotify::as.ggplot(grid.grabExpr(draw(ht_list , column_title = "Normal GFP samples", merge_legend = F, 
                                                  heatmap_legend_side = "bottom", annotation_legend_side = "bottom"), wrap = T))


ht_list2 = ggplotify::as.ggplot(grid.grabExpr( draw(ht_list2 , column_title = "Tumour samples", merge_legend = F, 
                                                    show_heatmap_legend=F), wrap = T))


#ht_list2 <- ht_list2+theme( panel.background = element_rect(fill = "white"),
 #               plot.margin = margin(l = 1, unit = "cm"))
####C gene plots####

#Sftpc#


dmrs <- readRDS(file.path(DATA, "DMRs_new_groups.RDS"))
dmr_1 <- makeGRangesFromDataFrame(dmrs$MsCCSP_MsSPC_control_comparison, keep.extra.columns = T)

source(paste0(CODE, 'gviz_region_plot_figures.R'))
plot_data <- readRDS(paste0(DATA, "/", "Sftpc", "cand_data.RDS"))
mypal = c(scgb1a1_col, alpha(scgb1a1_col, 0.5), sftpc_col,  alpha(sftpc_col, 0.5))

sftpc <- ggplotify::as.ggplot(grid.grabExpr( gviz_region_plot(plot_data[[1]], anno = plot_data[[2]], location=plot_data[[3]], genome="mm10", plot_data[[4]], 
                 dmrs = plot_data[[5]], highlight = subsetByOverlaps(dmr_1, plot_data[[3]])+20), wrap.grobs = T))
ggsave(plot=sftpc,filename = "methylation/sftpc_gene_viz.pdf", width = 7, height = 5)

#Scgb1a1#

plot_data <- readRDS(paste0(DATA, "/", "Scgb1a1", "cand_data.RDS"))

Scgb1a1 <- ggplotify::as.ggplot(grid.grabExpr( gviz_region_plot(plot_data[[1]], anno = plot_data[[2]], location=plot_data[[3]], genome="mm10", plot_data[[4]], 
                          dmrs = plot_data[[5]], highlight = subsetByOverlaps(dmr_1, plot_data[[3]])+20), wrap.grobs = T))
ggsave(plot=Scgb1a1,filename = "methylation/Scgb1a1_gene_viz.pdf", width = 7, height = 5)

#Lyz2#

plot_data <- readRDS(paste0(DATA, "/", "Lyz2", "cand_data.RDS"))

Lyz2 <- ggplotify::as.ggplot(grid.grabExpr( gviz_region_plot(plot_data[[1]], anno = plot_data[[2]], location=plot_data[[3]], genome="mm10", plot_data[[4]], 
                                                                dmrs = plot_data[[5]], highlight = subsetByOverlaps(dmr_1, plot_data[[3]])+20), wrap.grobs = T))
ggsave(plot=Lyz2,filename = "methylation/Lyz2_gene_viz.pdf", width = 7, height = 5)

#Lyz2#

plot_data <- readRDS(paste0(DATA, "/", "Scgb3a2", "cand_data.RDS"))

Scgb3a2 <- ggplotify::as.ggplot(grid.grabExpr( gviz_region_plot(plot_data[[1]], anno = plot_data[[2]], location=plot_data[[3]], genome="mm10", plot_data[[4]], 
                                                             dmrs = plot_data[[5]], highlight = subsetByOverlaps(dmr_1, plot_data[[3]])+20), wrap.grobs = T))
ggsave(plot=Scgb3a2,filename = "methylation/Scgb3a2_gene_viz.pdf", width = 7, height = 5)


###Pathway enrichment ####

# enrich <- list.files(path= paste0(RESULT, "/enrichment_res/"), pattern= "go_enrichment.*comparison.RDS", full.names = T)
# names <- gsub("(go_enrichment_)(loss|gain)_(.*)_comparison.RDS", "\\3_\\2", 
#               list.files(path= paste0(RESULT, "/enrichment_res/"), pattern= "go_enrichment.*comparison.RDS"))
# enrich_res <- lapply(enrich, readRDS)
# names(enrich_res) <- names



# 
# ####D Motif enrichment analysis#### 
# 
# ##from homer_new_groups
#  known_motif_results_merged <- readRDS( file = file.path(DATA, "figures/homer_results_merged.RDS"))
# # 
#  known_motif_results_merged <- known_motif_results_merged[rowSums(known_motif_results_merged[,-(1:2)])>quantile(rowSums(known_motif_results_merged[,-(1:2)]))[4],]
#  known_motif_results_merged$family[known_motif_results_merged$family=="box"] <- "Homeobox"
#  known_motif_results_merged$family[known_motif_results_merged$family=="forkhead"] <- "Forkhead"
#  known_motif_results_merged$family[known_motif_results_merged$family=="TEA"] <- "TEAD"
# # 
# # 
#  levels <- unique(known_motif_results_merged$family)
#  mypal <- pal_rickandmorty()(12)
#  mypal <- c(mypal, "#4cd137")
# # 
#  annotation_colors <- mypal
#  names(annotation_colors) <- levels
# # 
# 
# # 
#  ha <- HeatmapAnnotation("TF family"=known_motif_results_merged$family, which = "row", col=list("TF family"=annotation_colors),  annotation_name_gp  = gpar(fontsize = 10),
#                          annotation_legend_param  = list("TF family" = list(direction="horizontal", nrow=2, labels_gp = gpar(fontsize = 10))))
#  h1 <- HeatmapAnnotation(Direction=ifelse(grepl("gain", colnames(known_motif_results_merged[,-(1:2)])), "Gain in second", "Loss in second"), which = "column", 
#                          col=list(Direction=c("Gain in second"= "#e84118", "Loss in second"="#4cd137")), annotation_name_gp  = gpar(fontsize = 10),
#                          annotation_legend_param  = list("Direction" = list(direction="horizontal", nrow=1, labels_gp = gpar(fontsize = 10))))
# # 
#  col_labels <- colnames(known_motif_results_merged)[-(1:2)]
# # 
#  col_labels[grepl("MsCCSP_tumor_control_comparison", col_labels)] <-
#    "Scgb1a1 normal vs. tumour"
#  col_labels[grepl("MsSPC_tumor_control_comparison", col_labels)] <-
#    "Sftpc normal vs. tumour"
#  col_labels[grepl("MsCCSP_MsSPC_tumor", col_labels)] <-
#    "Scgb1a1 tumour vs. Sftpc tumour"
#  col_labels[grepl("MsCCSP_MsSPC_control", col_labels)] <-
#    "Scgb1a1 normal vs. Sftpc normal"
# # 
#  row_labels <- known_motif_results_merged$MotifName
#  row_labels[!(row_labels %in% c("SOX2", "RUNX", "NKX2.1", "FOXA1"))] <- ""
# # 
#  ht <- ComplexHeatmap::Heatmap(known_motif_results_merged[,-(1:2)], name = "-Log10 p",
#                          col = colorRamp2(c( 0,  max(known_motif_results_merged[,-(1:2)])), c( "white", "#5e4fa2")), 
#                          left_annotation = ha, top_annotation = h1,
#                          show_row_names = T, row_labels = row_labels, row_names_gp =gpar(fontsize = 8),
#                          show_column_names = T, column_labels = col_labels, column_names_gp =gpar(fontsize = 10), 
#                          heatmap_legend_param = list(direction = "horizontal"))
# # 
#  motifs <- ggplotify::as.ggplot(grid.grabExpr(draw(ht, merge_legend = F, heatmap_legend_side = "top", annotation_legend_side = "top"), wrap = T))
#ggsave(plot = motifs, filename = 
 #        "methylation/TFBS_enrichment.pdf", width = 8,  height=7)

##### LMC correlation ####

source(file.path(PATH,  "Reka/33_CoO_lung/scRNASeq/code/get_markers.R"))
#mypal = pal_igv()(27)
df <- readRDS(file=file.path(DATA, "medecom/control_correlation_LMCs.Rmd"))

gsea_rank <-  df %>%
  arrange(desc(abs(LMC3)),desc(LMC3_p)) %>%
  mutate(symbol=replace(symbol, duplicated(symbol), NA))%>%
  dplyr::filter(!is.na(symbol)) %>%
  arrange(LMC3)


LMC3_rank <- gsea_rank$LMC3
names(LMC3_rank) <- gsea_rank$symbol


fgseaRes <- fgsea(pathways = marker_lists_short,
                  stats    = LMC3_rank,
                  minSize  = 3,
                  maxSize  = 500, nperm=1000)

LMC3_plot <- gseaplot2(geneSetID = c(4,5,7,8), color = c(AT1_col,AT2_col, ciliated_col, club_col),  base_size = 11, rel_heights = c(1.5, 0.8), subplots = c(1, 3),  pvalue_table = FALSE, ES_geom = "line",
          title = "LMC3 correlation", ranks = LMC3_rank, geneSet = marker_lists_short, pvaluetable = fgseaRes)



gsea_rank <-  df %>%
  arrange(desc(abs(LMC1)),desc(LMC1_p)) %>%
  mutate(symbol=replace(symbol, duplicated(symbol), NA))%>%
  dplyr::filter(!is.na(symbol)) %>%
  arrange(LMC1)


lmc1_rank <- gsea_rank$LMC1
names(lmc1_rank) <- gsea_rank$symbol

fgseaRes <- fgsea(pathways = marker_lists_short,
                  stats    = lmc1_rank,
                  minSize  = 3,
                  maxSize  = 500, nperm=1000)

LMC1_plot <- gseaplot2(geneSetID = c(4,5,7,8), color =  c(AT1_col,AT2_col, ciliated_col, club_col),  base_size = 11, rel_heights = c(1.5, 0.8), subplots = c(1,3),  pvalue_table = FALSE, ES_geom = "line",
          title = "LMC1 correlation", ranks = lmc1_rank, geneSet = marker_lists_short, pvaluetable = fgseaRes)


# 
# gsea_rank <-  df %>%
#   arrange(desc(abs(LMC3)),desc(LMC3_p)) %>%
#   mutate(symbol=replace(symbol, duplicated(symbol), NA))%>%
#   dplyr::filter(!is.na(symbol)) %>%
#   arrange(LMC3)
# 
# 
# lmc3_rank <- gsea_rank$LMC3
# names(lmc3_rank) <- gsea_rank$symbol
# 
# fgseaRes <- fgsea(pathways = marker_lists_short,
#                   stats    = lmc3_rank,
#                   minSize  = 3,
#                   maxSize  = 500, nperm=1000)
# 
# LMC3_plot <-  gseaplot2(geneSetID = c(4,5,7,8), color = c("green", "red", "blue", "black"),  base_size = 11, rel_heights = c(1.5, 0.5), subplots = c(1,3),  pvalue_table = FALSE, ES_geom = "line",
#           title = "LMC3 correlation", ranks = lmc3_rank, geneSet = marker_lists_short, pvaluetable = fgseaRes)

# 
# volcano_plot_correlation <- function(data, x, y, label){
#   data<-  data %>%
#     arrange(desc(abs(get(x))),desc(get(y))) %>%
#     mutate(new_labels=replace(get(label), duplicated(get(label)), NA))
#   g <- ggplot(data, aes_string(x=x, y=y, label="new_labels")) + theme_bw()+ xlab(paste0("Correlation with ", x))+ ylab("-log10(p value)")+
#     geom_point()+geom_label_repel(max.overlaps = Inf, )+ggtitle(paste0("Best correlating genes, ", x))
#   
#   g
# }
# 
# df <- readRDS(file=file.path(DATA, "medecom/control_correlation_LMCs.Rmd"))
# 
# df %>%
#   arrange(desc(abs(LMC1)),desc(LMC1_p)) %>% 
# mutate(symbol=replace(symbol, duplicated(symbol), NA)) %>%
#   
# volcano_1 <- df %>%
#   arrange(desc(abs(LMC1)),desc(LMC1_p)) %>% 
#   mutate(AT2_markers=replace(AT2_markers, duplicated(AT2_markers), NA)) %>%
#   mutate(club_markers=replace(club_markers, duplicated(club_markers), NA)) %>%
#   mutate(AT1_markers=replace(AT1_markers, duplicated(AT1_markers), NA)) %>%
#   dplyr::select(LMC1, LMC1_p, AT2_markers, club_markers, AT1_markers) %>%
#   tidyr::gather(key="markers", value="labels",   AT2_markers, club_markers, AT1_markers) %>%
#   mutate(markers = recode(markers, AT2_markers = "AT2 markers", AT1_markers = "AT1 markers", club_markers = "Club markers")) %>%
#   ggplot( aes(x=LMC1, y=LMC1_p, label=labels)) + theme_bw()+ facet_grid(cols=vars(markers))+
#   xlab(paste0("Correlation coefficient"))+ ylab("-log10(p value)")+
#   geom_point()+geom_label_repel(max.overlaps = Inf, )+ggtitle("Correlation with LMC1")
# 
# volcano_2 <- df %>%
#   arrange(desc(abs(LMC2)),desc(LMC2_p)) %>% 
#   mutate(AT2_markers=replace(AT2_markers, duplicated(AT2_markers), NA)) %>%
#   mutate(club_markers=replace(club_markers, duplicated(club_markers), NA)) %>%
#   mutate(AT1_markers=replace(AT1_markers, duplicated(AT1_markers), NA)) %>%
#   dplyr::select(LMC2, LMC2_p, AT2_markers, club_markers, AT1_markers) %>%
#   tidyr::gather(key="markers", value="labels",   AT2_markers, club_markers, AT1_markers) %>%
#   mutate(markers = recode(markers, AT2_markers = "AT2 markers", AT1_markers = "AT1 markers", club_markers = "Club markers")) %>%
#   ggplot( aes(x=LMC2, y=LMC2_p, label=labels)) + theme_bw()+ facet_grid(cols=vars(markers))+
#   xlab(paste0("Correlation coefficient"))+ ylab("-log10(p value)")+
#   geom_point()+geom_label_repel(max.overlaps = Inf, )+ggtitle("Correlation with LMC2")
# 
# volcano_3 <- df %>%
#   arrange(desc(abs(LMC3)),desc(LMC3_p)) %>% 
#   mutate(AT2_markers=replace(AT2_markers, duplicated(AT2_markers), NA)) %>%
#   mutate(club_markers=replace(club_markers, duplicated(club_markers), NA)) %>%
#   mutate(AT1_markers=replace(AT1_markers, duplicated(AT1_markers), NA)) %>%
#   dplyr::select(LMC3, LMC3_p, AT2_markers, club_markers, AT1_markers) %>%
#   tidyr::gather(key="markers", value="labels",   AT2_markers, club_markers, AT1_markers) %>%
#   mutate(markers = recode(markers, AT2_markers = "AT2 markers", AT1_markers = "AT1 markers", club_markers = "Club markers")) %>%
#   ggplot( aes(x=LMC3, y=LMC3_p, label=labels)) + theme_bw()+ facet_grid(cols=vars(markers))+
#   xlab(paste0("Correlation coefficient"))+ ylab("-log10(p value)")+
#   geom_point()+geom_label_repel(max.overlaps = Inf, )+ggtitle("Correlation with LMC3")
# 


#############
#DMR heatmaps

res <- readRDS(paste0(DATA_meth, "no_snps_methrix.RDS"))

res <- methrix::subset_methrix(res, samples = rownames(res@colData)[-which(res@colData$full_name=="MsCCSP_control01")])
res <- methrix::subset_methrix(res, samples =  rownames(res@colData)[-which(res@colData$full_name=="MsCCSP_control05")])
res <- methrix::subset_methrix(res, samples =  rownames(res@colData)[-which(res@colData$full_name=="MsHOPX_control01")])
res <- methrix::subset_methrix(res, samples =  rownames(res@colData)[-grep("rfp", res@colData$full_name)])


new_groups <-  readRDS(paste0(DATA, "new_tumor_groups_HOPX.rds"))



DMRs <- readRDS(file.path(DATA, "DMRs_new_groups.RDS"))
#DMLs <- readRDS(file.path(DATA, "DMLs_new_groups.RDS"))
comparisons_all <- names(DMRs)

regions_all <- list()
for (comp in comparisons_all) {
  
  regions <- makeGRangesFromDataFrame(DMRs[[comp]][order(abs(DMRs[[comp]]$diff.Methy), decreasing = T),], keep.extra.columns = T)
 # regions$cells <- substr(comp, 1, 3)
  regions_all[[comp]] <- regions
  
}


samples <- rownames(res@colData)
samples_final <- samples[grep("cont", samples)]
samples_final <- c(samples_final, new_groups$SPC, new_groups$CCSP)

samples_final <- as.data.frame(samples_final)
samples_final$Origin <- NA
samples_final$Origin[grep("CCSP", samples_final$samples_final)] <- "Scgb1a1"
samples_final$Origin[grep("SPC", samples_final$samples_final)] <- "Sftpc"
samples_final$Origin[grep("HOPX", samples_final$samples_final)] <- "Hopx"
#samples_final$Origin[samples_final$samples_final %in% new_groups$SPC] <- "AT2 origin"
#samples_final$Origin[samples_final$samples_final %in% new_groups$CCSP] <- "Club origin"
samples_final$Type <- NA
samples_final$Type[grep("control", samples_final$samples_final)] <- "Normal"
samples_final$Type[grep("tumor", samples_final$samples_final)] <- "Tumour"




res2 <- subset_methrix(res, samples = samples_final$samples_final)


mat1 <- as.data.frame(get_region_summary(res2, regions_all[["MsCCSP_MsSPC_control_comparison"]], overlap_type = "any"))
#mat1$cells <- regions_all$cells
#cells <- regions_all$cells[complete.cases(mat1)]
plot_mat <- mat1[complete.cases(mat1),-(1:5)]
plot_mat <- plot_mat[,samples_final$samples_final]


samples_final$Type <- factor(samples_final$Type)

samples_final$order <- paste0(samples_final$Origin, samples_final$Type)
samples_final$order[(nrow(samples_final)-4):nrow(samples_final)] <- "Club"
order <- factor(samples_final$order, 
                levels=c("HopxNormal","Scgb1a1Normal","SftpcNormal",  
                          "SftpcTumour", "Club"))
labels <- c("Hopx", "Scgb1a1", "Sftpc", "AT2", "Club")
lab_color <-c(hopx_col, scgb1a1_col, sftpc_col, sftpc_col, scgb1a1_col)

ha_all <- HeatmapAnnotation("Origin" = anno_block(gp = gpar(fill =  lab_color), labels =labels, 
                                                  labels_gp = gpar( fontsize = 10)),
                            "Lineage" = samples_final$Origin,
                            "Normal/Tumor" = samples_final$Type,
                            col = list("Lineage"= c("Sftpc" = sftpc_col, Scgb1a1 = scgb1a1_col, Hopx=hopx_col), 
                                       "Normal/Tumor"=c(Normal = normal_col,  "Tumour" = tumor_col)), 
                            annotation_legend_param = list("Lineage" = list(nrow = 1), 
                                                           "Normal/Tumor" = list(nrow = 1)), 
                            annotation_name_gp = gpar(col="black", fontsize = 10), 
                            show_legend = c(FALSE, FALSE))


meth_ht <- Heatmap(plot_mat, name = "Methylation",cluster_rows = TRUE, 
        show_row_names = FALSE, show_column_names = FALSE, cluster_columns = F,
        column_title = "Cell type specific methylation",
        top_annotation = ha_all, 
        #heatmap_legend_param = list(direction = "horizontal"), 
        col = mycols, 
        column_split = order)
draw(meth_ht, heatmap_legend_side = "right", annotation_legend_side = "bottom")

meth_heatmap <- ggplotify::as.ggplot(grid.grabExpr(draw(meth_ht, 
                                                        heatmap_legend_side = "right", 
                                                        padding = unit(c(1, 2, 2, 2), "cm")), wrap = F))
ggsave(plot = meth_heatmap, "methylation/heatmap.pdf", width = 8, height=7)


###########

g <- plot_grid(plot_grid(ht_list, plot_grid(ht_list2, NULL, nrow=2, rel_heights = c(1,0.23)), 
                         nrow=1, rel_widths = c(1,0.7)), NULL, nrow = 2, rel_heights = c(1, 0.2))
g1 <- plot_grid(pca_plot, g, rel_widths = c(1,1,0.6), nrow=1)

g2_1 <- plot_grid(sftpc, Scgb1a1,  rel_heights  = c(1,1), ncol=1)

g2_2 <- plot_grid(LMC1_plot, LMC3_plot, nrow = 2)


g_f <- plot_grid(g1, plot_grid(g2_2, NULL, g2_1, NULL, nrow=1, rel_widths = c(1,0.3,1.2, 0.3)), 
                 nrow=2, rel_heights=c(0.7, 1))


#ggsave(pca_plot, filename = "Figure_2_1.tiff")
#ggsave(ht_list, filename = "Figure_2_2.tiff")
#ggsave(ht_list2, filename = "Figure_2_3.tiff")
#ggsave(sftpc, filename = "Figure_2_4.tiff")
#ggsave(Scgb1a1, filename = "Figure_2_5.tiff")
#ggsave(g2_2, filename = "Figure_2_6.tiff")


ggsave(g_f, filename = "methylation/Figure_2_norfp.pdf", width = 17, height = 15)
ggsave(g_f, filename = "methylation/Figure_2_norfp.tiff", width = 17, height = 15)



