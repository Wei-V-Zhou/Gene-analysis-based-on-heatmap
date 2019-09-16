#=======================================================================================#
# The usages of all packages:                                                           #
#       gdata: provide various R programming tools for data manipulation                #
#      scater: provide a class and numerous functions for the QC and data-proccess      #
#      gplots: load a heatmap packages:heatmap.2                                        #
#     ggplot2: create an elegant data visualization using the grammar of graphics       #
#    pheatmap: draw clustered heatmaps where one has better control over  parameters    #
# geneplotter: graphics related functions                                               #
# SingleCellExperiment: store single-cell sequencing data for SingleCellExperiment class#
#=======================================================================================#

##################
# Load libraries #
##################
# clear objectives and garbage collection
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)
# load packages
if(T){
  library(gdata)
  library(gplots)
  library(scater)
  library(export)
  library(ggplot2)
  library(pheatmap)
  library(geneplotter)
}

# Load variation-stablized RNA-seq data
if(!file.exists("m4T1.primary.emt.Rdata")){
  load("m4T1.D4.sceset.qc.new.Rdata")
  exprs <- log2(calculateCPM(counts(sce.all)) +1)
  exprs_norm <- t(scale(t(exprs)))
  save(ann, exprs_norm, file = "m4T1.primary.emt.Rdata")
}
load("m4T1.primary.emt.Rdata")
n <- exprs_norm
n[n > 2] = 2
n[n < -1.5] = -1.5
vitro <- n
vitro.ann <- ann

# Re-arrange sample order for presentation
# vitro.ann$fig.3.plot.ind <- c(1,2,3,7,8,9,10,11,12,19,20,21,22,23,24,13,14,15,16,17,18,4,5,6)

# EMT-related genes curated from the literature
if(F){
  # from https://github.com/Xiang-HF-Zhang/Dichotomous-of-innate-immune-landscape
  emt.e <- c("Cdh1", "Epcam", "Cldn2", "Cldn4", "Krt7", "Krt8", "Krt18", "Krt19", "Krt20", "Esrp1")
  emt.m <- c("Pdgfrb", "Zeb1","Fap", "Cdh2", "Cdh11", "Col1a1", "Fn1", "Twist1", "Snai1", "Snai2")
  emt.em <- c(emt.e, emt.m)
  vitro.emt <- vitro[which(is.element(rownames(vitro), emt.em)), ]
}
if(F){
  # from 10.1038/bjc.2014.80
  emt.e <- c("Cdh1","Krt18")
  emt.m <- c("Zeb1","Zeb2","Snai1","Snai2","Cdh2","Vim","Twist1")
  emt.em <- c(emt.e, emt.m)
  vitro.emt <- vitro[which(is.element(rownames(vitro), emt.em)), ]
}
if(F){
  # from 10.1038/s41586-019-1526-3
  emt.em1 <- c("Spp1","Nt5e","Cyr61","Pmepa1","Tnfrsf12a","Areg","Cxcl1","Sdc4","Rhob","Hmga2")
  vitro.emt <- vitro[which(is.element(rownames(vitro), emt.em1)), ]
}
if(F){
  # from 10.1172/JCI36183 (Note: Cdh1, Tjp1, Col4a1 for attenuated markers)
  emt.em2 <- c("Cdh1","Cdh11","Itgb3","Sdc1","S100a4","Acta1","Vim","Ctnnb1","Col3a1","Tjp1",
               "Fn1","Snai1","Snai2","Zeb1","Twist1","Lef1","Ets1","Col1a1","Foxc2","Col4a1")
  vitro.emt <- vitro[which(is.element(rownames(vitro), emt.em2)), ]
}
# EMT-related genes choosed from the literature and GSEA
if(T){
  emt.e <- c("Cdh1", "Epcam", "Cldn4", "Krt7", "Krt19", "Esrp1")
  emt.m <- c("Zeb1", "Lef1", "Vim", "Cfp", "Csf1r", "Dab2", "Fcgr2b","Snx20", "Spi1")
  # emt.m <- c("Zeb1", "Col1a1", "Col3a1","Itgb3", "Lef1", "Vim", "Angptl2", "Aox1")
  # emt.group <- c("Cd19", "Cd68", "Ly6g")
  emt.tg1 <- c("Htra1")
  emt.tg2 <- c("Sdc4")
  emt.tg3 <- c("Jag1")
  emt.em <- c(emt.e, emt.m)
  vitro.emt <- vitro[which(is.element(rownames(vitro), emt.em)), ]
}

# Choose the corresponding datasets to analyse
vitro.data <- rbind(vitro.emt, t(ann))
colnames(vitro.data) <- vitro.data[nrow(vitro.data), ]
{
  emt.CL <- as.matrix(vitro.data[ , colnames(vitro.data)=="CellLine"])
  emt.PT <- as.matrix(vitro.data[ , colnames(vitro.data)=="PrimaryTumor"])
  emt.D4 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D4"])
  emt.D10 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D10"])
  emt.D16 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D16"])
}
{
  colnames(emt.CL) <- emt.CL[(nrow(emt.CL)-1), ]
  emt.CL.1 <- emt.CL[c(-nrow(emt.CL),-(nrow(emt.CL)-1)), ]
  colnames(emt.PT) <- emt.PT[(nrow(emt.PT)-1), ]
  emt.PT.1 <- emt.PT[c(-nrow(emt.PT),-(nrow(emt.PT)-1)), ]
  colnames(emt.D4) <- emt.D4[(nrow(emt.D4)-1), ]
  emt.D4.1 <- emt.D4[c(-nrow(emt.D4),-(nrow(emt.D4)-1)), ]
  colnames(emt.D10) <- emt.D10[(nrow(emt.D10)-1), ]
  emt.D10.1 <- emt.D10[c(-nrow(emt.D10),-(nrow(emt.D10)-1)), ]
  colnames(emt.D16) <- emt.D16[(nrow(emt.D16)-1), ]
  emt.D16.1 <- emt.D16[c(-nrow(emt.D16),-(nrow(emt.D16)-1)), ]
}

# Compare the PrimaryTumor, D4 and D16
if(T){
  emt.PT_D4_D16 <- as.matrix(cbind(emt.PT.1, emt.D4.1, emt.D16.1))
  PT_D4_D16 <- apply(emt.PT_D4_D16, 2, as.numeric)
  row.names(PT_D4_D16) <- row.names(emt.PT_D4_D16)
  anno <- vitro.ann[colnames(PT_D4_D16), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col) <- "Species"
  
  # B cell cluster
  {
    rm(vitro.emt, vitro.data)
    emt.group <- "Cd19"
    vitro.emt <- vitro[which(is.element(rownames(vitro), emt.group)), ]
    # Choose the corresponding datasets to analyse
    vitro.data <- rbind(vitro.emt, t(ann))
    colnames(vitro.data) <- vitro.data[nrow(vitro.data), ]
    {
      emt.CL <- as.matrix(vitro.data[ , colnames(vitro.data)=="CellLine"])
      emt.PT <- as.matrix(vitro.data[ , colnames(vitro.data)=="PrimaryTumor"])
      emt.D4 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D4"])
      emt.D10 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D10"])
      emt.D16 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D16"])
    }
    {
      colnames(emt.CL) <- emt.CL[(nrow(emt.CL)-1), ]
      emt.CL.1 <- emt.CL[c(-nrow(emt.CL),-(nrow(emt.CL)-1)), ]
      colnames(emt.PT) <- emt.PT[(nrow(emt.PT)-1), ]
      emt.PT.1 <- emt.PT[c(-nrow(emt.PT),-(nrow(emt.PT)-1)), ]
      colnames(emt.D4) <- emt.D4[(nrow(emt.D4)-1), ]
      emt.D4.1 <- emt.D4[c(-nrow(emt.D4),-(nrow(emt.D4)-1)), ]
      colnames(emt.D10) <- emt.D10[(nrow(emt.D10)-1), ]
      emt.D10.1 <- emt.D10[c(-nrow(emt.D10),-(nrow(emt.D10)-1)), ]
      colnames(emt.D16) <- emt.D16[(nrow(emt.D16)-1), ]
      emt.D16.1 <- emt.D16[c(-nrow(emt.D16),-(nrow(emt.D16)-1)), ]
    }
  }
  Cd19.PT_D4_D16 <- t(c(emt.PT.1, emt.D4.1, emt.D16.1))
  row.names(Cd19.PT_D4_D16) <- emt.group
  anno <- vitro.ann[colnames(PT_D4_D16), ]
  anno$CellType <- "Cd19"
  anno_Cd19 <- colnames(Cd19.PT_D4_D16)[Cd19.PT_D4_D16 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Cd19, 2] = "Others"
  annotation_col_Cd19 = data.frame(factor(anno_df$CellType))
  row.names(annotation_col_Cd19) <- anno_df[ , 1]
  colnames(annotation_col_Cd19) <- "Cd19"
  
  # Macrophage cluster
  {
    rm(vitro.emt, vitro.data)
    emt.group <- "Cd68"
    vitro.emt <- vitro[which(is.element(rownames(vitro), emt.group)), ]
    # Choose the corresponding datasets to analyse
    vitro.data <- rbind(vitro.emt, t(ann))
    colnames(vitro.data) <- vitro.data[nrow(vitro.data), ]
    {
      emt.CL <- as.matrix(vitro.data[ , colnames(vitro.data)=="CellLine"])
      emt.PT <- as.matrix(vitro.data[ , colnames(vitro.data)=="PrimaryTumor"])
      emt.D4 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D4"])
      emt.D10 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D10"])
      emt.D16 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D16"])
    }
    {
      colnames(emt.CL) <- emt.CL[(nrow(emt.CL)-1), ]
      emt.CL.1 <- emt.CL[c(-nrow(emt.CL),-(nrow(emt.CL)-1)), ]
      colnames(emt.PT) <- emt.PT[(nrow(emt.PT)-1), ]
      emt.PT.1 <- emt.PT[c(-nrow(emt.PT),-(nrow(emt.PT)-1)), ]
      colnames(emt.D4) <- emt.D4[(nrow(emt.D4)-1), ]
      emt.D4.1 <- emt.D4[c(-nrow(emt.D4),-(nrow(emt.D4)-1)), ]
      colnames(emt.D10) <- emt.D10[(nrow(emt.D10)-1), ]
      emt.D10.1 <- emt.D10[c(-nrow(emt.D10),-(nrow(emt.D10)-1)), ]
      colnames(emt.D16) <- emt.D16[(nrow(emt.D16)-1), ]
      emt.D16.1 <- emt.D16[c(-nrow(emt.D16),-(nrow(emt.D16)-1)), ]
    }
  }
  Cd68.PT_D4_D16 <- t(c(emt.PT.1, emt.D4.1, emt.D16.1))
  row.names(Cd68.PT_D4_D16) <- emt.group
  anno <- vitro.ann[colnames(PT_D4_D16), ]
  anno$CellType <- "Cd68"
  anno_Cd68 <- colnames(Cd68.PT_D4_D16)[Cd68.PT_D4_D16 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Cd68, 2] = "Others"
  annotation_col_Cd68 = data.frame(factor(anno_df$CellType))
  row.names(annotation_col_Cd68) <- anno_df[ , 1]
  colnames(annotation_col_Cd68) <- "Cd68"
  
  # Neutrophil cluster
  {
    rm(vitro.emt, vitro.data)
    emt.group <- "Ly6g"
    vitro.emt <- vitro[which(is.element(rownames(vitro), emt.group)), ]
    # Choose the corresponding datasets to analyse
    vitro.data <- rbind(vitro.emt, t(ann))
    colnames(vitro.data) <- vitro.data[nrow(vitro.data), ]
    {
      emt.CL <- as.matrix(vitro.data[ , colnames(vitro.data)=="CellLine"])
      emt.PT <- as.matrix(vitro.data[ , colnames(vitro.data)=="PrimaryTumor"])
      emt.D4 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D4"])
      emt.D10 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D10"])
      emt.D16 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D16"])
    }
    {
      colnames(emt.CL) <- emt.CL[(nrow(emt.CL)-1), ]
      emt.CL.1 <- emt.CL[c(-nrow(emt.CL),-(nrow(emt.CL)-1)), ]
      colnames(emt.PT) <- emt.PT[(nrow(emt.PT)-1), ]
      emt.PT.1 <- emt.PT[c(-nrow(emt.PT),-(nrow(emt.PT)-1)), ]
      colnames(emt.D4) <- emt.D4[(nrow(emt.D4)-1), ]
      emt.D4.1 <- emt.D4[c(-nrow(emt.D4),-(nrow(emt.D4)-1)), ]
      colnames(emt.D10) <- emt.D10[(nrow(emt.D10)-1), ]
      emt.D10.1 <- emt.D10[c(-nrow(emt.D10),-(nrow(emt.D10)-1)), ]
      colnames(emt.D16) <- emt.D16[(nrow(emt.D16)-1), ]
      emt.D16.1 <- emt.D16[c(-nrow(emt.D16),-(nrow(emt.D16)-1)), ]
    }
  }
  Ly6g.PT_D4_D16 <- t(c(emt.PT.1, emt.D4.1, emt.D16.1))
  row.names(Ly6g.PT_D4_D16) <- emt.group
  anno <- vitro.ann[colnames(PT_D4_D16), ]
  anno$CellType <- "Ly6g"
  anno_Ly6g <- colnames(Ly6g.PT_D4_D16)[Ly6g.PT_D4_D16 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Ly6g, 2] = "Others"
  annotation_col_Ly6g = data.frame(factor(anno_df$CellType))
  row.names(annotation_col_Ly6g) <- anno_df[ , 1]
  colnames(annotation_col_Ly6g) <- "Ly6g"
  
  annotat_col = data.frame(c(annotation_col, annotation_col_Cd19,
                             annotation_col_Cd68, annotation_col_Ly6g))
  row.names(annotat_col) = anno_df[ , 1]
  annot_color = list(Species = c(PrimaryTumor = "#82B7FF", BoneMet_D4 = "#FF9289", BoneMet_D16 = "#00D65C"),
                     Cd19 = c(Cd19 = "red", Others = "white"),
                     Cd68 = c(Cd68 = "red", Others = "white"),
                     Ly6g = c(Ly6g = "red", Others = "white"))
  dev.off()
  dev.new()
  pheatmap(PT_D4_D16, treeheight_col = 0, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2",
           color = colorRampPalette(c("green", "black", "red"))(50), annotation_col = annotat_col, annotation_colors = annot_color)
  # img=pheatmap(PT_D4_D16, treeheight_col = 30, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2",
  #            color = colorRampPalette(c("green", "black", "red"))(50), annotation_col = annotation_col)
  # graph2ppt(x = img, file = "EMT heatmap", width = 830, height = 445)
}

# Choose the Cd19 groups for annotation colorbar
{
  rm(vitro.emt, vitro.data)
  emt.group <- "Cd19"
  vitro.emt <- vitro[which(is.element(rownames(vitro), emt.group)), ]
  # Choose the corresponding datasets to analyse
  vitro.data <- rbind(vitro.emt, t(ann))
  colnames(vitro.data) <- vitro.data[nrow(vitro.data), ]
  {
    emt.CL <- as.matrix(vitro.data[ , colnames(vitro.data)=="CellLine"])
    emt.PT <- as.matrix(vitro.data[ , colnames(vitro.data)=="PrimaryTumor"])
    emt.D4 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D4"])
    emt.D10 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D10"])
    emt.D16 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D16"])
  }
  {
    colnames(emt.CL) <- emt.CL[(nrow(emt.CL)-1), ]
    emt.CL.1 <- emt.CL[c(-nrow(emt.CL),-(nrow(emt.CL)-1)), ]
    colnames(emt.PT) <- emt.PT[(nrow(emt.PT)-1), ]
    emt.PT.1 <- emt.PT[c(-nrow(emt.PT),-(nrow(emt.PT)-1)), ]
    colnames(emt.D4) <- emt.D4[(nrow(emt.D4)-1), ]
    emt.D4.1 <- emt.D4[c(-nrow(emt.D4),-(nrow(emt.D4)-1)), ]
    colnames(emt.D10) <- emt.D10[(nrow(emt.D10)-1), ]
    emt.D10.1 <- emt.D10[c(-nrow(emt.D10),-(nrow(emt.D10)-1)), ]
    colnames(emt.D16) <- emt.D16[(nrow(emt.D16)-1), ]
    emt.D16.1 <- emt.D16[c(-nrow(emt.D16),-(nrow(emt.D16)-1)), ]
  }
  Cd19.PT_D4_D16 <- t(c(emt.PT.1, emt.D4.1, emt.D16.1))
  row.names(Cd19.PT_D4_D16) <- emt.group
  anno <- vitro.ann[colnames(PT_D4_D16), ]
  anno$CellType <- "Cd19"
  anno_Cd19 <- colnames(Cd19.PT_D4_D16)[Cd19.PT_D4_D16 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Cd19, 2] = NA
  annotation_col = data.frame(anno_df$CellType)
  # ann_colors = list(CellType = "#7570B3")
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(PT_D4_D16, treeheight_col = 30, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2",
           color = colorRampPalette(c("green", "black", "red"))(50), annotation_col = annotation_col)
  
}
# Choose the Cd68 groups for annotation colorbar
{
  rm(vitro.emt, vitro.data)
  emt.group <- "Cd68"
  vitro.emt <- vitro[which(is.element(rownames(vitro), emt.group)), ]
  # Choose the corresponding datasets to analyse
  vitro.data <- rbind(vitro.emt, t(ann))
  colnames(vitro.data) <- vitro.data[nrow(vitro.data), ]
  {
    emt.CL <- as.matrix(vitro.data[ , colnames(vitro.data)=="CellLine"])
    emt.PT <- as.matrix(vitro.data[ , colnames(vitro.data)=="PrimaryTumor"])
    emt.D4 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D4"])
    emt.D10 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D10"])
    emt.D16 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D16"])
  }
  {
    colnames(emt.CL) <- emt.CL[(nrow(emt.CL)-1), ]
    emt.CL.1 <- emt.CL[c(-nrow(emt.CL),-(nrow(emt.CL)-1)), ]
    colnames(emt.PT) <- emt.PT[(nrow(emt.PT)-1), ]
    emt.PT.1 <- emt.PT[c(-nrow(emt.PT),-(nrow(emt.PT)-1)), ]
    colnames(emt.D4) <- emt.D4[(nrow(emt.D4)-1), ]
    emt.D4.1 <- emt.D4[c(-nrow(emt.D4),-(nrow(emt.D4)-1)), ]
    colnames(emt.D10) <- emt.D10[(nrow(emt.D10)-1), ]
    emt.D10.1 <- emt.D10[c(-nrow(emt.D10),-(nrow(emt.D10)-1)), ]
    colnames(emt.D16) <- emt.D16[(nrow(emt.D16)-1), ]
    emt.D16.1 <- emt.D16[c(-nrow(emt.D16),-(nrow(emt.D16)-1)), ]
  }
  Cd68.PT_D4_D16 <- t(c(emt.PT.1, emt.D4.1, emt.D16.1))
  row.names(Cd68.PT_D4_D16) <- emt.group
  anno <- vitro.ann[colnames(PT_D4_D16), ]
  anno$CellType <- "Cd68"
  anno_Cd68 <- colnames(Cd68.PT_D4_D16)[Cd68.PT_D4_D16 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Cd68, 2] = NA
  annotation_col = data.frame(anno_df$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(PT_D4_D16, treeheight_col = 30, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2",
           color = colorRampPalette(c("green", "black", "red"))(50), annotation_col = annotation_col)
  
}
# Choose the Ly6g groups for annotation colorbar
{
  rm(vitro.emt, vitro.data)
  emt.group <- "Ly6g"
  vitro.emt <- vitro[which(is.element(rownames(vitro), emt.group)), ]
  # Choose the corresponding datasets to analyse
  vitro.data <- rbind(vitro.emt, t(ann))
  colnames(vitro.data) <- vitro.data[nrow(vitro.data), ]
  {
    emt.CL <- as.matrix(vitro.data[ , colnames(vitro.data)=="CellLine"])
    emt.PT <- as.matrix(vitro.data[ , colnames(vitro.data)=="PrimaryTumor"])
    emt.D4 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D4"])
    emt.D10 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D10"])
    emt.D16 <- as.matrix(vitro.data[ , colnames(vitro.data)=="BoneMet_D16"])
  }
  {
    colnames(emt.CL) <- emt.CL[(nrow(emt.CL)-1), ]
    emt.CL.1 <- emt.CL[c(-nrow(emt.CL),-(nrow(emt.CL)-1)), ]
    colnames(emt.PT) <- emt.PT[(nrow(emt.PT)-1), ]
    emt.PT.1 <- emt.PT[c(-nrow(emt.PT),-(nrow(emt.PT)-1)), ]
    colnames(emt.D4) <- emt.D4[(nrow(emt.D4)-1), ]
    emt.D4.1 <- emt.D4[c(-nrow(emt.D4),-(nrow(emt.D4)-1)), ]
    colnames(emt.D10) <- emt.D10[(nrow(emt.D10)-1), ]
    emt.D10.1 <- emt.D10[c(-nrow(emt.D10),-(nrow(emt.D10)-1)), ]
    colnames(emt.D16) <- emt.D16[(nrow(emt.D16)-1), ]
    emt.D16.1 <- emt.D16[c(-nrow(emt.D16),-(nrow(emt.D16)-1)), ]
  }
  Ly6g.PT_D4_D16 <- t(c(emt.PT.1, emt.D4.1, emt.D16.1))
  row.names(Ly6g.PT_D4_D16) <- emt.group
  anno <- vitro.ann[colnames(PT_D4_D16), ]
  anno$CellType <- "Ly6g"
  anno_Ly6g <- colnames(Ly6g.PT_D4_D16)[Ly6g.PT_D4_D16 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Ly6g, 2] = NA
  annotation_col = data.frame(anno_df$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(PT_D4_D16, treeheight_col = 30, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2",
           color = colorRampPalette(c("green", "black", "red"))(50), annotation_col = annotation_col)
  
}





# Compare the CellLine and Primary Tumor
if(T){
  emt.CL_PT <- as.matrix(cbind(emt.CL.1, emt.PT.1))
  CL_PT <- apply(emt.CL_PT, 2, as.numeric)
  row.names(CL_PT) <- row.names(emt.CL_PT)
  anno <- vitro.ann[colnames(CL_PT), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(CL_PT,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}
# Compare the D4, D10 and D16
if(T){
  emt.D10_D16 <- as.matrix(cbind(emt.D10.1, emt.D16.1))
  D10_D16 <- apply(emt.D10_D16, 2, as.numeric)
  row.names(D10_D16) <- row.names(emt.D10_D16)
  anno <- vitro.ann[colnames(D10_D16), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(D10_D16,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}
# Compare the D4, D10 and D16
if(T){
  emt.D4_D10_D16 <- as.matrix(cbind(emt.D4.1, emt.D10.1, emt.D16.1))
  D4_D10_D16 <- apply(emt.D4_D10_D16, 2, as.numeric)
  row.names(D4_D10_D16) <- row.names(emt.D4_D10_D16)
  anno <- vitro.ann[colnames(D4_D10_D16), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(D4_D10_D16,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}
# Compare the PT, D4, D10 and D16
if(T){
  emt.PT_D4_D10_D16 <- as.matrix(cbind(emt.PT.1, emt.D4.1, emt.D10.1, emt.D16.1))
  PT_D4_D10_D16 <- apply(emt.PT_D4_D10_D16, 2, as.numeric)
  row.names(PT_D4_D10_D16) <- row.names(emt.PT_D4_D10_D16)
  anno <- vitro.ann[colnames(PT_D4_D10_D16), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(PT_D4_D10_D16,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}
# Compare the CellLine, D4, D10 and D16
if(T){
  emt.CL_D4_D10_D16 <- as.matrix(cbind(emt.CL.1, emt.D4.1, emt.D10.1, emt.D16.1))
  CL_D4_D10_D16 <- apply(emt.CL_D4_D10_D16, 2, as.numeric)
  row.names(CL_D4_D10_D16) <- row.names(emt.CL_D4_D10_D16)
  anno <- vitro.ann[colnames(CL_D4_D10_D16), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(CL_D4_D10_D16,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}


# Compare the CellLine and D4
if(T){
  emt.CL_D4 <- as.matrix(cbind(emt.CL.1, emt.D4.1))
  CL_D4 <- apply(emt.CL_D4, 2, as.numeric)
  row.names(CL_D4) <- row.names(emt.CL_D4)
  anno <- vitro.ann[colnames(CL_D4), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(CL_D4,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}
# Compare the CellLine and D10
if(T){
  emt.CL_D10 <- as.matrix(cbind(emt.CL.1, emt.D10.1))
  CL_D10 <- apply(emt.CL_D10, 2, as.numeric)
  row.names(CL_D10) <- row.names(emt.CL_D10)
  anno <- vitro.ann[colnames(CL_D10), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(CL_D10,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}
# Compare the CellLine and D16
if(T){
  emt.CL_D16 <- as.matrix(cbind(emt.CL.1, emt.D16.1))
  CL_D16 <- apply(emt.CL_D16, 2, as.numeric)
  row.names(CL_D16) <- row.names(emt.CL_D16)
  anno <- vitro.ann[colnames(CL_D16), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(CL_D16,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}

# Compare the PrimaryTumor and D4
if(T){
  emt.PT_D4 <- as.matrix(cbind(emt.PT.1, emt.D4.1))
  PT_D4 <- apply(emt.PT_D4, 2, as.numeric)
  row.names(PT_D4) <- row.names(emt.PT_D4)
  anno <- vitro.ann[colnames(PT_D4), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(PT_D4, treeheight_col = 30, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2",
           color = colorRampPalette(c("green", "black", "red"))(50), annotation_col = annotation_col)
}
# Compare the PrimaryTumor and D10
if(T){
  emt.PT_D10 <- as.matrix(cbind(emt.PT.1, emt.D10.1))
  PT_D10 <- apply(emt.PT_D10, 2, as.numeric)
  row.names(PT_D10) <- row.names(emt.PT_D10)
  anno <- vitro.ann[colnames(PT_D10), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(PT_D10,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}
# Compare the PrimaryTumor and D16
if(T){
  emt.PT_D16 <- as.matrix(cbind(emt.PT.1, emt.D16.1))
  PT_D16 <- apply(emt.PT_D16, 2, as.numeric)
  row.names(PT_D16) <- row.names(emt.PT_D16)
  anno <- vitro.ann[colnames(PT_D16), ]
  annotation_col = data.frame(anno$CellType)
  row.names(annotation_col) <- anno[ , 1]
  colnames(annotation_col)<-"Species"
  pheatmap(PT_D16,treeheight_col = 30,treeheight_row = 0,show_colnames = F,clustering_method="ward.D",
           color = colorRampPalette(c("green", "black", "red"))(50),annotation_col = annotation_col)
}


#============================#
#       Musician: Resonance  #
#           Date: 2019/09/06 #
# Revised author: Resonance  #
#           Time: 2019/09/15 #
#============================#