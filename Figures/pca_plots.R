
library(methrix)
library(dmrseq)
library(bsseq)
library(HDF5Array)
library(plotly)
library(autoplotly)
library(pheatmap)
library(ggsci)
library(annotatr)
library(rtracklayer)

###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="V:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="N:/"} else {
    PATH_Y="/C010-datasets/"}

DATA = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/data/")
CODE = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/code/")
RESULT = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/output/")
CALLS = paste0(PATH_Y, "External/2018-10-Sotillo/data/methylDackel/")
SAVE_PLOT <- paste0(PATH, "Reka/33_CoO_lung/Figures/methylation/pca/")


################################################


methrix_pca_full <- function (m, var = "top", top_var = 1000, ranges = NULL, 
                              pheno = NULL, do_plot = TRUE, n_pc = 2) 
{
  var_select <- match.arg(var, c("top", "rand"))
  if (!is(m, "methrix")) {
    stop("A valid methrix object needs to be supplied.")
  }
  if (!(is.numeric(top_var) & is.numeric(n_pc))) {
    stop("Either top_var or n_pc variables are not numeric.")
  }
  if (!is.null(ranges)) {
    message("GenomicRanges will be used for the PCA")
    meth_sub <- subset_methrix(m = m, regions = ranges)
    meth_sub <- methrix::get_matrix(m = meth_sub, type = "M", 
                                    add_loci = FALSE)
  }
  if (is.null(top_var)) {
    message("All CpGs in the dataset will be used for the PCA")
    if (is.null(ranges)) {
      meth_sub <- get_matrix(m = m, type = "M", add_loci = FALSE)
    }
  }
  else {
    if (!is.numeric(top_var)) {
      stop("top_var must be numeric.")
    }
    top_var <- as.integer(as.character(top_var))
    if (var_select == "rand") {
      if (!is.null(ranges)) {
        message("Random CpGs within provided GRanges will be used for the PCA")
        ids <- sample(x = seq_along(meth_sub), replace = FALSE, 
                      size = min(top_var, nrow(meth_sub)))
      }
      else {
        message("Random CpGs will be used for the PCA")
        ids <- sample(x = seq_along(m), replace = FALSE, 
                      size = as.integer(as.character(min(top_var, 
                                                         nrow(m)))))
      }
      meth_sub <- get_matrix(m = m[ids, ], type = "M", 
                             add_loci = FALSE)
    }
    else {
      if (!is.null(ranges)) {
        if (is_h5(m)) {
          sds <- DelayedMatrixStats::rowSds(meth_sub, 
                                            na.rm = TRUE)
        }
        else {
          sds <- matrixStats::rowSds(meth_sub, na.rm = TRUE)
        }
        meth_sub <- meth_sub[order(sds, decreasing = TRUE)[seq_len(min(top_var, 
                                                                       nrow(meth_sub)))], ]
      }
      else {
        meth_sub <- methrix::get_matrix(m = order_by_sd(m)[seq_len(min(top_var, 
                                                                       nrow(m)))], type = "M", add_loci = FALSE)
      }
    }
  }
  meth_sub <- meth_sub[complete.cases(meth_sub), , drop = FALSE]
  if (nrow(meth_sub) == 0) {
    stop("Zero loci available post NA removal :(")
  }
  meth_pca <- prcomp(x = t(meth_sub), retx = TRUE)
  n_pc <- ncol(meth_pca$x)
  pc_vars <- meth_pca$sdev^2/sum(meth_pca$sdev^2)
  names(pc_vars) <- colnames(meth_pca$x)
  pc_vars <- round(pc_vars, digits = 2)
  if (do_plot) {
    par(mar = c(3, 4, 1, 4))
    b = barplot(height = pc_vars, names.arg = NA, col = "#2c3e50", 
                ylim = c(0, 1), las = 2, axes = FALSE, ylab = "Variance Explained")
    points(x = b, y = cumsum(pc_vars), type = "l", 
           lty = 2, lwd = 1.2, xpd = TRUE, col = "#c0392b")
    points(x = b, y = cumsum(pc_vars), type = "p", 
           pch = 19, xpd = TRUE, col = "#c0392b")
    mtext(text = paste0("PC", 1:length(pc_vars)), side = 1, 
          at = b, las = 2, line = 0.5, cex = 0.8)
    axis(side = 2, at = seq(0, 1, 0.1), line = 0, las = 2, 
         cex.axis = 0.8)
    axis(side = 4, at = seq(0, 1, 0.1), line = 0, las = 2, 
         cex.axis = 0.8)
    legend(x = "topleft", legend = "Cumulative", 
           col = "#c0392b", pch = 19, lwd = 1, cex = 0.75, 
           bty = "n")
  }
  results <- list(PC_matrix = meth_pca$x, var_explained = pc_vars)
  if (do_plot) {
    plot_pca(pca_res = results, m = m, col_anno = pheno)
  }
  gc(verbose = FALSE)
  return(meth_pca)
}


###########################


res_no_snps <- readRDS(paste0(DATA, "no_snps_methrix.RDS"))


res_no_snps <- methrix::subset_methrix(res_no_snps, 
                                       samples = rownames(res_no_snps@colData)[-which(res_no_snps@colData$full_name %in% c("MsCCSP_control01", 
                                                                                                                           "MsCCSP_control05", 
                                                                                                                           "MsHOPX_control01"))])
res_no_snps <- methrix::subset_methrix(res_no_snps, 
                                       samples = rownames(res_no_snps@colData)[-grep("rfp", rownames(res_no_snps@colData))])

Cellpaper <- c("#313695","#4575b4", "#74add1", "#abd9e9","#e0f3f8","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026")

#mypal = pal_nejm()(7)
#ann_colors = list(
#  sample_name = c(control = mypal[1], tumor = mypal[2], "tumor-rfp"=mypal[3]),
#  cell_type = c(MsCCSP = mypal[4], MsSPC = mypal[5], MsHOPX=mypal[6], MsKrt5=mypal[7])
#)

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


source(file.path(CODE, 'plot_pca_manual.R'))



res_no_snps@colData$cell_type[res_no_snps@colData$cell_type=="MsCCSP"] <- "Scgb1a1"
res_no_snps@colData$cell_type[res_no_snps@colData$cell_type=="MsSPC"] <- "Sftpc"
res_no_snps@colData$cell_type[res_no_snps@colData$cell_type=="MsHOPX"] <- "Hopx"
res_no_snps@colData$sample_name[res_no_snps@colData$sample_name=="control"] <- "Normal"
res_no_snps@colData$sample_name[res_no_snps@colData$sample_name=="tumor"] <- "Tumor"



pca_res <- methrix_pca_full(res_no_snps, top_var = 30000, pheno = "cell_type", do_plot = F)
p1 <- autoplot(pca_res, data = as.data.frame(res_no_snps@colData), shape="sample_name", colour="cell_type", size=4)+
  theme_bw()+
  scale_color_manual(name="Line", values = c(hopx_col, scgb1a1_col, sftpc_col))+ 
  stat_ellipse(aes(fill = sample_name), geom = "polygon", type="norm")+
  scale_fill_manual(name="N/T", values = c(alpha(normal_col, 0.2), alpha(tumor_col, 0.2)))+
  ggtitle("All sites")+labs(shape="N/T")
ggsave(plot = p1, filename = file.path(SAVE_PLOT, "all_sites.pdf"), width = 7, height = 5)


#annots_gr_0 = import.bed(con=file.path(PATH_Y, "/External/2018-10-Sotillo/data/chromhmm/lung_0_mm10_15_posterior.bed"), extraCols = c(Type="character", First="numeric", Second="character", Position1="integer", Position2="integer", Color="character"))

annots_gr_0 = import.bed(con=file.path(PATH_Y, "/External/2018-10-Sotillo/data/chromhmm/P0_lung_15_segments.bed.gz"), extraCols = c(V4="character", Type="character"))

tssa <- subset_methrix(res_no_snps, regions = annots_gr_0[annots_gr_0$Type=="Pr-A",])
pca_res <- methrix_pca_full(tssa, top_var = 100000, pheno = "cell_type", do_plot = T)
p2 <- autoplot(pca_res, data = as.data.frame(res_no_snps@colData), shape="sample_name", colour="cell_type", size=4)+
  theme_bw()+
  scale_color_manual(name="Line", values = c(hopx_col, scgb1a1_col, sftpc_col))+ 
  stat_ellipse(aes(fill = sample_name), geom = "polygon", type="norm")+
  scale_fill_manual(name="N/T", values = c(alpha(normal_col, 0.2), alpha(tumor_col, 0.2)))+
  ggtitle("Active promoters")+labs(shape="N/T")
ggsave(plot = p2, filename = file.path(SAVE_PLOT, "active_promoters.pdf"), width = 7, height = 5)


#res_no_snps_XY_no_MsCCSP_tumor01 <- res_no_snps_XY[,res_no_snps_XY@colData$full_name!="MsCCSP_tumor01"]

#tssbiv <- subset_methrix(res_no_snps_XY_no_MsCCSP_tumor01, regions = annots_gr_0[annots_gr_0$Type=="TssBiv",])
tssbiv <- subset_methrix(res_no_snps, regions = annots_gr_0[annots_gr_0$Type=="Pr-B",])
pca_res <- methrix_pca_full(tssbiv, top_var = 100000, pheno = "cell_type", do_plot = T)
p3 <- autoplot(pca_res, data = as.data.frame(res_no_snps@colData), shape="sample_name", colour="cell_type", size=4)+
  theme_bw()+
  scale_color_manual(name="Line", values = c(hopx_col, scgb1a1_col, sftpc_col))+ 
  stat_ellipse(aes(fill = sample_name), geom = "polygon", type="norm")+
  scale_fill_manual(name="N/T", values = c(alpha(normal_col, 0.2), alpha(tumor_col, 0.2)))+
  ggtitle("Bivalent promoters")+labs(shape="N/T")
ggsave(plot = p3, filename = file.path(SAVE_PLOT, "bivalent_promoters.pdf"), width = 7, height = 5)



enhancer <- subset_methrix(res_no_snps, regions = annots_gr_0[annots_gr_0$Type %in% c("En-Sd", "En-Sp"),])
#enhancer <- subset_methrix(res_no_snps_XY, regions = annots_gr_0[annots_gr_0$Type=="Enh",])
pca_res <- methrix_pca_full(enhancer, top_var = 100000, pheno = "cell_type", do_plot = F)
p4 <- autoplot(pca_res, data = as.data.frame(res_no_snps@colData), shape="sample_name", colour="cell_type", size=4)+
  theme_bw()+
  scale_color_manual(name="Line", values = c(hopx_col, scgb1a1_col, sftpc_col))+ 
  stat_ellipse(aes(fill = sample_name), geom = "polygon", type="norm")+
  scale_fill_manual(name="N/T", values = c(alpha(normal_col, 0.2), alpha(tumor_col, 0.2)))+
  ggtitle("Strong enhancers")+labs(shape="N/T")
ggsave(plot = p4, filename = file.path(SAVE_PLOT, "strong_enhancers.pdf"), width = 7, height = 5)


#enhancer <- subset_methrix(res_no_snps_XY, regions = annots_gr_0[annots_gr_0$Type=="EnhPois1" | annots_gr_0$Type=="EnhPois2" ,])
enhancer <- subset_methrix(res_no_snps, regions = annots_gr_0[annots_gr_0$Type %in% c("En-Pd", "En-Pp") ,])
pca_res <- methrix_pca_full(enhancer, top_var = 100000, pheno = "cell_type", do_plot = F)
p5 <- autoplot(pca_res, data = as.data.frame(res_no_snps@colData), shape="sample_name", colour="cell_type", size=4)+
  theme_bw()+
  scale_color_manual(name="Line", values = c(hopx_col, scgb1a1_col, sftpc_col))+ 
  stat_ellipse(aes(fill = sample_name), geom = "polygon", type="norm")+
  scale_fill_manual(name="N/T", values = c(alpha(normal_col, 0.2), alpha(tumor_col, 0.2)))+
  ggtitle("Poised enhancers")+labs(shape="N/T")
ggsave(plot = p5, filename = file.path(SAVE_PLOT, "poised_enhancers.pdf"), width = 7, height = 5)





