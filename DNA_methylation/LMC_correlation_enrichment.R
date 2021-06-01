

library(DOSE)
library(ernrichplot)
library(fgsea)
library(ggsci)


####colors####

sftpc_col <- "#20854EFF"
scgb1a1_col <- "#E18727FF"
hopx_col <- "#6F99ADFF"

normal_col <- "#7876B1FF"
tumor_col <- "#EE4C97FF"

tdTom_col = "#e84118"
GFP_col = "#4cd137"


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

################functions#########

CapStr <- function(y) {
      paste(toupper(substring(y, 1,1)), substring(y, 2),
                         sep="")
   }


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
                        size = 1)

  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8),
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))

  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
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
  p.pos <- p.pos + ylab("Correlation Coefficient") + xlab("Rank in Ordered Dataset") +
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2,
                               l = 0.2, unit = "cm"))+scale_x_reverse(expand = c(0, 0))

  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)

  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(labels= paste0(c(pvaluetable[geneSetID,"pathway"]$pathway), ", p value= ", round(c(pvaluetable[geneSetID , "pval"])$pval, digits = 3)),
                                        values = color)
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


###################

source(file.path(PATH,  "Reka/33_CoO_lung/scRNASeq/code/get_markers.R"))

df <- readRDS(file=file.path(DATA, "medecom/control_correlation_LMCs.Rmd"))

gsea_rank <-  df %>%
  arrange(desc(abs(LMC1)),desc(LMC1_p)) %>%
  mutate(symbol=replace(symbol, duplicated(symbol), NA))%>%
  dplyr::filter(!is.na(symbol)) %>%
  arrange(LMC1)


lmc1_rank <- gsea_rank$LMC1
names(lmc1_rank) <- gsea_rank$symbol


fgseaRes <- fgsea(pathways = marker_lists_short,
                  stats    = lmc1_rank,
                  minSize  = 15,
                  maxSize  = 500, nperm=1000)

gseaplot2(geneSetID = c(4,5,7,8), color = c("green", "red", "blue", "black"),  base_size = 11, rel_heights = c(1.5, 0.5, 1), subplots = 1:3,  pvalue_table = FALSE, ES_geom = "line",
          title = "LMC1 correlation", ranks = lmc1_rank, geneSet = marker_lists_short, pvaluetable = fgseaRes)



gsea_rank <-  df %>%
  arrange(desc(abs(LMC2)),desc(LMC2_p)) %>%
  mutate(symbol=replace(symbol, duplicated(symbol), NA))%>%
  dplyr::filter(!is.na(symbol)) %>%
  arrange(LMC2)


lmc2_rank <- gsea_rank$LMC2
names(lmc2_rank) <- gsea_rank$symbol

fgseaRes <- fgsea(pathways = marker_lists_short,
                  stats    = lmc2_rank,
                  minSize  = 15,
                  maxSize  = 500, nperm=1000)

gseaplot2(geneSetID = c(4,5,7,8), color = c("green", "red", "blue", "black"),  base_size = 11, rel_heights = c(1.5, 0.5, 1), subplots = 1:3,  pvalue_table = FALSE, ES_geom = "line",
          title = "LMC2 correlation", ranks = lmc2_rank, geneSet = marker_lists_short, pvaluetable = fgseaRes)



gsea_rank <-  df %>%
  arrange(desc(abs(LMC3)),desc(LMC3_p)) %>%
  mutate(symbol=replace(symbol, duplicated(symbol), NA))%>%
  dplyr::filter(!is.na(symbol)) %>%
  arrange(LMC3)


lmc3_rank <- gsea_rank$LMC3
names(lmc3_rank) <- gsea_rank$symbol

fgseaRes <- fgsea(pathways = marker_lists_short,
                  stats    = lmc3_rank,
                  minSize  = 15,
                  maxSize  = 500, nperm=1000)

gseaplot2(geneSetID = c(4,5,7,8), color = c("green", "red", "blue", "black"),  base_size = 11, rel_heights = c(1.5, 0.2, 0.5), subplots = 1:3,  pvalue_table = FALSE, ES_geom = "line",
          title = "LMC3 correlation", ranks = lmc3_rank, geneSet = marker_lists_short, pvaluetable = fgseaRes)


#fgseaRes <- fgsea(pathways = marker_lists_short,
#stats    = lmc1_rank,
#eps      = 0.0,
#minSize  = 15,
#maxSize  = 500)
mypal = pal_igv()(26)


lung_cell_atlas_other_markers <- read.delim("V:/Reka/33_CoO_lung/CoO_Lung_Cancer/data/lung_cell_atlas_other_markers.txt")
new_list <- lapply(1:ncol(lung_cell_atlas_other_markers), function(x)  CapStr(tolower(lung_cell_atlas_other_markers[,x])))
new_list <- lapply(new_list, function(x) x[1:28])
new_list <- lapply(new_list, function(x) x[x %in% gsea_rank$symbol])
names(new_list) <- colnames(lung_cell_atlas_other_markers)


gseaplot2(geneSetID = c(1:13), color = mypal[1:13],  base_size = 11, rel_heights = c(1.5, 0.2, 0.5), subplots = 1:3,  pvalue_table = FALSE, ES_geom = "line",
title = "LMC3 correlation", ranks = lmc3_rank, geneSet = new_list)


fgseaRes <- fgsea(pathways = new_list,
stats    = lmc3_rank,
minSize  = 15,
maxSize  = 500, nperm=1000)

gseaplot2(geneSetID = c(3:5), color = mypal[1:3],  base_size = 11, rel_heights = c(1.5, 0.2, 0.5), subplots = 1:3,  pvalue_table = FALSE, ES_geom = "line",
                     title = "LMC3 correlation", ranks = lmc3_rank, geneSet = new_list, pvaluetable = fgseaRes)
