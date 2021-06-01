
library(methrix)
library(data.table)
library(HDF5Array)
library(SummarizedExperiment)
library(ggplot2)
library(rtracklayer)
library(limma)
library(knitr)
library(pheatmap)
library(MeDeCom)
library(ggsci)

###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="V:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="N:/"} else {
    PATH_Y="/C010-datasets/"}

DATA = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/data/norfp/")
DATA_meth = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/data/")
RESULT = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/output/")
CALLS = paste0(PATH_Y, "External/2018-10-Sotillo/data/methylDackel/")
CODE = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/code/")
HOMER_RESULTS <- paste0(RESULT, "homer_res/enhancers/")


#
mat  <- as.data.frame(readRDS(paste0(DATA, "medecom/bivalent_enhancers.RDS")))
rowsds <- rowSds(as.matrix(mat[,-(1:3)]), na.rm=T)
mat <- mat[order(rowsds, decreasing = T)[1:100000],]
mat <- mat[complete.cases(mat),]
#mat_control <- mat[,grep("control", colnames(mat[,-(1:3)]))]

medecom.result2 <- MeDeCom::runMeDeCom(as.matrix(mat[,-(1:3)]), 2:6, c(0, 10^(-5:-1)), NINIT = 10, NFOLDS = 10, ITERMAX = 300, NCORES = 8)
saveRDS(medecom.result2, file = paste0(DATA, "medecom/poised_enhancers_new_samples.RDS"))

medecom.result2 <-  readRDS( file = paste0(DATA, "medecom/poised_enhancers_new_samples.RDS"))
MeDeCom::plotParameters(MeDeComSet = medecom.result2)


K_sel <- 4
lambda_sel <- 0.001
proportions <- MeDeCom::getProportions(medecom.result2, K=K_sel, lambda=lambda_sel)
colnames(proportions) <- colnames(mat[,-(1:3)])

LMCs <- MeDeCom::getLMCs(medecom.result2, K=K_sel, lambda=lambda_sel)
colnames(LMCs) <- paste0("LMC", 1:K_sel)
LMCs <- cbind(mat[,1:3], LMCs)



annotation <- data.frame(staining=gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\1", colnames(mat[,-(1:3)])),
                         type= gsub("(Ms[[:alnum:]]+)_(tumor|control)[[:digit:]]+(.rfp)?", "\\2", colnames(mat[,-(1:3)])),
                         origin = ifelse(grepl("rfp", colnames(mat[,-(1:3)])), "RFP", "GFP"))
rownames(annotation) <- colnames(mat[,-(1:3)])


p <- pheatmap(proportions,  show_colnames = T, legend = F, annotation_legend = F, fontsize = 12, annotation_col = annotation)

new_groups <- list(SPC=colnames(proportions)[proportions["LMC4",]>0.1 & grepl("tumor", colnames(proportions))],
  CCSP=colnames(proportions)[proportions["LMC1",]>0.1 & grepl("tumor", colnames(proportions))])

saveRDS(proportions, file=file.path(DATA, "medecom_proportions_all_HOPX.RDS"))
saveRDS(new_groups, paste0(DATA, "new_tumor_groups_HOPX.RDS"))
