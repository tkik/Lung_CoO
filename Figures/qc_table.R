

#library(RnBeads.mm10)
library(RnBeads)
library(data.table)
library(ggplot2)
library(plotly)
library(xlsx)

###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="V:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="N:/"} else {
    PATH_Y="/C010-datasets/"}

DATA = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/data/")
RESULT = paste0(PATH, "Reka/33_CoO_lung/CoO_Lung_Cancer/output/")
CALLS = paste0(PATH_Y, "External/2018-10-Sotillo/data/methylDackel/")
SNPS = paste0(PATH_Y, "External/2018-10-Sotillo/data/snps/")



load(file=file.path(DATA, "QC_data.RData"))


final <- global_methylation_sum[global_methylation_sum$context=="CH",c("sample_name", "conversion_rate")]
final <- cbind(final, global_methylation_sum[global_methylation_sum$context=="CG",c("coverage", "ratio")])
colnames(final) <- c("Sample Name", "Conversion rate", "Global covergage", "Global methylation")

write.xlsx(as.data.frame(final),  row.names = F, col.names=T,file="Supplementary_tables/QC_data.xlsx")
