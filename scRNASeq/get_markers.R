library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(SingleR)
library(harmony)
library(readxl)
library(cerebroApp)

###########libraries and functions#############
if (grepl("Windows", Sys.getenv("OS"))){
  PATH ="V:/"} else {
    PATH ="/C010-projects/"}
if (grepl("Windows", Sys.getenv("OS"))){
  PATH_Y="N:/"} else {
    PATH_Y="/C010-datasets/"}
PATH_ICGC <- "/icgc/"


DATA_sc = paste0(PATH, "Reka/33_CoO_lung/scRNASeq/data")



Krt8_cells<- read_excel(paste0(DATA_sc, "/Krt8_cell_populations_2020_17358_MOESM6_ESM.xlsx"),
                        sheet = "cell_types_2", col_types = c("skip",
                                                              "text", "numeric", "numeric", "numeric",
                                                              "numeric", "text"))


marker_lists <-  Krt8_cells %>%
  group_split(cluster) %>%
  lapply(function(x) pull(x, gene))


names(marker_lists) <-  sort(unique(Krt8_cells$cluster))


marker_lists$activated_club_fig_5 <-  c("Cst3", "H2-Ab1", "Cd74", "H2-Eb1", "Psap", "Ctss", "Vim", "H2-Aa", "Laptm5", "Actb", "Alox5ap", "Srgn", "Ccl17", "Lsp1", "Scgb3a1", "S100a4")

marker_lists$AT2_sig_genes <- c("Sftpc", "Sftpb", "Sftpd", "Pgc", "Cldn18", "Aqp4", "Scgb3a1", "Abca3", "Gata6", "Nkx2-1", "Sfta3", "Igfbp2", "Hopx", "Napsa", "Foxa2", "Ager", "Lamp1")

marker_lists$Focal <- c("Itgb1", "Cdc42", "Col4a1", "Rock2", "Flna", "Lamc2", "Rhoa", "Col6a1", "Ctnnb1")

marker_lists$fourwk_vs_end <- c("Scgb1a1", "Scgb3a2", "Lyz2", "Cyp2f2", "B2m", "H2-K1")



p53_pathway <- c("Btg2", "Fos", "Ifi30", "Ctsd", "Ier3")


Myc_Targets_V1 <- c("Ywhae", "Npm1", "HNRNPA3", "HSP90AB1", "RPL34", "CNBP", "RPL22", "HSPE1", "YWHAQ", "PCBP1", "CANX", "HNRNPA2B1", "SERBP1", "SRSF2", "ERH", "PABPC1", "RPS2", "EIF4G2")

VEGFA_VEGFRm<- c("YWHAE", "ITGB1", "HSP90AA1", "ROCK1", "ROCK2", "HSPB1", "IQGAP1", "RHOA", "CDC42", "PPP3CA", "NCL", "TXNIP", "CTNNB1", "DNAJB9", "IGFBP7", "EZR", "PAK2", "HSPA1A")

Myc_Targets_V1 <-gsub("(^[a-z])", "\\U\\1", tolower(Myc_Targets_V1), perl = T)
VEGFA_VEGFRm <- gsub("(^[a-z])", "\\U\\1", tolower(VEGFA_VEGFRm), perl = T)

tumor_spec <- c("Meg3", "Areg", "Pf4", "Kng2", "Id2", "Serpine1", "Cdkn1a", "Rian", "Dmkn", "Tspan8",
                "Napsa",
                "Cald1",
                "Ly6c1",
                "Tmem213",
                "Basp1",
                "Thbs1",
                "Igfbp7",
                "Kcnc3",
                "Ly6c2",
                "Tnfrsf12a",
                "Arg2",
                "Malt1",
                "Stbd1")

marker_lists <- append(marker_lists, list(p53_pathway=p53_pathway, Myc_Targets_V1=Myc_Targets_V1, VEGFA_VEGFRm=VEGFA_VEGFRm, Myc_Targets_V1=Myc_Targets_V1, tumor_spec=tumor_spec))
marker_lists_short <- lapply(marker_lists, function(x) x[1:min(30, length(x))])



marker_lists_short[["Club progenitor"]] <- c("Cd74", "Cd14", "H2-K1", "AW112010")

marker_lists_short[["Mixed AT1_AT2"]] <- c("Clic3", "Ager", "Col4a3", "Dpysl2", "Mmp11", "Akap5", "Rtkn2", "Hopx", "Radil", "Sec14l3")
marker_lists_short[["Highly mixed"]] <- c("Epcam", "Cldn7", "Cldn4", "Cars", "Hmces", "Dusp10", "Slc4a11Lif", "Tigit", "Lgals3")
marker_lists_short[["EMT program"]] <- c("Lamb1", "Fbln2", "Fxyd5", "Tnnt2", "Inhba", "Dbn1", "Cpe", "Axl", "Igfbp4","Vim")
marker_lists_short[["Activated Club cells"]] <- c("Cst3", "H2-Ab1", "H2-Eb1", "Psap", "Ctss", "H2-Aa", "Laptm5", "Actb", "Alox5ap",
                                                  "Srgn", "Ccl17", "Lsp1", "Scgb3a1", "S100a4")


