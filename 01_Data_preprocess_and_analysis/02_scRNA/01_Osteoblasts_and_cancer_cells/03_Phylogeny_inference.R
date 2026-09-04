library(here) # project-root-relative paths; run scripts from repository root
##################### Phylogeny inference (related to Figure 2D,E, Figures S2E,S3A) ##################### 
# For uphyloplot2 determination of phylogeny ----- 
## Set to correct work directory
## (run from repository root)

## STEP I # Set to correct work directory
## (run from repository root)

## STEP II - Delete non-malignant cells 
## sed '/^ImmuneCells/d' <  17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings > trimmed_infercnv.cell_groupings
## Move trimmed_infercnv.cell_groupings to the following directory 

## STEP III - Change work directory to the following where UPhyloplot2 script is stored 
## (run from repository root)

## STEP IV Plot
## python uphyloplot2.py -c 0  

## STEPV Using code below for annotation and cell percentage calculations

# Per-sample Ward.D phylogeny -----
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny"))
library(Seurat)
library(SCP)

# Load data ----- 
OB_cancer <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Malig_OB_for_NMF_sample_w_too_few_cells_removed.rds"))
Idents(OB_cancer) <- "sample"

DefaultAssay(OB_cancer) <- "RNA"
sample_names <- names(table(OB_cancer$sample))
patient_names <- gsub("_AF_", "", sample_names)
patient_names <- gsub("_BF_", "", patient_names) %>% unique()
cluster_group_all <- data.frame()

# To include both DimPlots and phylogeny tree (to include number of cells)
# Clonal structure determination is based on inferCNV heatmaps, UMAP plot structure, and the specific CNV events, amongst other factors


# Run for cells with treatment dynamics ----- 
# HP1 ----- 
H_P1 <- subset(OB_cancer, idents = c("H_BF_P1", "H_AF_P1"))
DimPlot(H_P1)

H_P1[["RNA"]]$scale.data <- NULL
H_P1[["RNA"]]$data <- NULL
H_P1[["RNA"]] <- split(H_P1[["RNA"]], f = H_P1$sample) 
H_P1[["RNA"]]
DefaultAssay(H_P1) <- "RNA"
H_P1 <- NormalizeData(H_P1) %>% FindVariableFeatures() %>% ScaleData()
#H_P1[["percent.ribo"]] <- PercentageFeatureSet(object = H_P1, pattern = "^RP[SL]")
#Variable_features <- H_P1@assays[["RNA"]]@meta.data[["var.features"]][!is.na(H_P1@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(H_P1)[grep("^MT_", rownames(H_P1), invert=F)] 
#RB_genes <- rownames(H_P1)[grep("^RP[SL]", rownames(H_P1), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#H_P1 <- ScaleData(H_P1, vars.to.regress = c("percent.mt", "percent.ribo"))
#H_P1 <- ScaleData(H_P1)
H_P1 <- RunPCA(H_P1)
H_P1 <- FindNeighbors(H_P1, dims = 1:50, reduction = "pca")
H_P1 <- RunUMAP(H_P1, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(H_P1, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "HP1")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
#cnv_table1_m <- as.matrix(cnv_table1)
#cnv_table2_m <- as.matrix(cnv_table2)
#cnv_table3_m <- as.matrix(cnv_table3)
#cnv_table4_m <- as.matrix(cnv_table4)
#cnv_table5_m <- as.matrix(cnv_table5)

DimPlot(H_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(H_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(H_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(H_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(H_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#   "#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#  1       0       0       0       0    5212
#  2    1405       0       0       0       0
#  3       0       0       0     325       0
#  4       0       0     415       0       0
#  5       0    1508       0       0       0

cluster1 <- "C1"
cluster2 <- "C3"
cluster3 <- "C2"
cluster4 <- "C2"
cluster5 <- "C1"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
H_P1 <- AddMetaData(H_P1, cluster_group2)
DimPlot(H_P1, group.by = "Clone", reduction = "umap.unintegrated")
table(H_P1$Clone, H_P1$treatment_status)
length(which(H_P1$Clone == "C1")) / length(H_P1$Clone) # 0.7580372
length(which(H_P1$Clone == "C2")) / length(H_P1$Clone) # 0.08347434
length(which(H_P1$Clone == "C3")) / length(H_P1$Clone) # 0.1584884
saveRDS(H_P1, "H_P1_cloned.rds")
# H_P1 <- readRDS("H_P1_cloned.rds")

Idents(H_P1) <- "Clone"
levels(H_P1) <- c("C1", "C2", "C3")
H_P1$Clone <- Idents(H_P1)
p1 <- CellDimPlot(
  srt = H_P1, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933", "#339999"), 
  raster = T, 
  pt.size = 3, 
  split.by = "treatment_status", label_insitu = T) + theme(legend.position = "none")
p1

pdf("HP1_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

#CellDimPlot3D(srt = H_P1, group.by = "Clone")
# remove object when done 
rm(H_P1)

# HP3----- 
H_P3 <- subset(OB_cancer, idents = c("H_BF_P3", "H_AF_P3"))
DimPlot(H_P3)

H_P3[["RNA"]]$scale.data <- NULL
H_P3[["RNA"]]$data <- NULL
H_P3[["RNA"]] <- split(H_P3[["RNA"]], f = H_P3$sample) 
H_P3[["RNA"]]
DefaultAssay(H_P3) <- "RNA"
H_P3 <- NormalizeData(H_P3) %>% FindVariableFeatures() %>% ScaleData()
#H_P3[["percent.ribo"]] <- PercentageFeatureSet(object = H_P3, pattern = "^RP[SL]")
#Variable_features <- H_P3@assays[["RNA"]]@meta.data[["var.features"]][!is.na(H_P3@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(H_P3)[grep("^MT_", rownames(H_P3), invert=F)] 
#RB_genes <- rownames(H_P3)[grep("^RP[SL]", rownames(H_P3), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#H_P3 <- ScaleData(H_P3, vars.to.regress = c("percent.mt", "percent.ribo"))
#H_P3 <- ScaleData(H_P3)
H_P3 <- RunPCA(H_P3)
H_P3 <- FindNeighbors(H_P3, dims = 1:50, reduction = "pca")
H_P3 <- RunUMAP(H_P3, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
H_P3 <- JoinLayers(H_P3)
DimPlot(H_P3, label = T, reduction = "umap.unintegrated")
DimPlot(H_P3, group.by = "treatment_status", reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "HP3")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)

DimPlot(H_P3, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(H_P3, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(H_P3, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(H_P3, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(H_P3, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#   1       0       0       0       0    5045
#   2       0       0       0    1748       0
#   3       0       0     756       0       0
#   4       0    2523       0       0       0
#   5    2100       0       0       0       0

cluster1 <- "C3"
cluster2 <- "C3"
cluster3 <- "C1"
cluster4 <- "C3"
cluster5 <- "C2"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
H_P3 <- AddMetaData(H_P3, cluster_group2)
DimPlot(H_P3, group.by = "Clone", reduction = "umap.unintegrated")
table(H_P3$Clone, H_P3$treatment_status)
length(which(H_P3$Clone == "C1")) / length(H_P3$Clone) # 0.06210976
length(which(H_P3$Clone == "C2")) / length(H_P3$Clone) # 0.06210976
length(which(H_P3$Clone == "C3")) / length(H_P3$Clone) # 0.06210976

FeaturePlot(H_P3, "MYC", reduction = "umap.unintegrated")
FeaturePlot(H_P3, "IBSP", reduction = "umap.unintegrated")
VlnPlot(H_P3, "MYC", group.by = "Clone")
VlnPlot(H_P3, "IBSP", group.by = "Clone")
Idents(H_P3) <- "Clone"
markers <- FindMarkers(H_P3, ident.1 = "C3", only.pos = F)

Idents(H_P3) <- "Clone"
levels(H_P3) <- c("C1", "C2", "C3")
H_P3$Clone <- Idents(H_P3)
p1 <- CellDimPlot(
  srt = H_P3, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933", "#339999"), 
  raster = T, 
  pt.size = 3, 
  split.by = "treatment_status", label_insitu = T) + 
  theme(legend.position = "none")
p1

pdf("HP3_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

saveRDS(H_P3, "H_P3_cloned.rds")
# remove object when done 
rm(H_P3)

# HP5----- 
H_P5 <- subset(OB_cancer, idents = c("H_BF_P5", "H_AF_P5"))
DimPlot(H_P5)

H_P5[["RNA"]]$scale.data <- NULL
H_P5[["RNA"]]$data <- NULL
H_P5[["RNA"]] <- split(H_P5[["RNA"]], f = H_P5$sample) 
H_P5[["RNA"]]
DefaultAssay(H_P5) <- "RNA"
H_P5 <- NormalizeData(H_P5) %>% FindVariableFeatures() %>% ScaleData()
#H_P5[["percent.ribo"]] <- PercentageFeatureSet(object = H_P5, pattern = "^RP[SL]")
#Variable_features <- H_P5@assays[["RNA"]]@meta.data[["var.features"]][!is.na(H_P5@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(H_P5)[grep("^MT_", rownames(H_P5), invert=F)] 
#RB_genes <- rownames(H_P5)[grep("^RP[SL]", rownames(H_P5), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#H_P5 <- ScaleData(H_P5, vars.to.regress = c("percent.mt", "percent.ribo"))
#H_P5 <- ScaleData(H_P5)
H_P5 <- RunPCA(H_P5)
H_P5 <- FindNeighbors(H_P5, dims = 1:50, reduction = "pca")
H_P5 <- RunUMAP(H_P5, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
H_P5 <- JoinLayers(H_P5)
DimPlot(H_P5, label = T, reduction = "umap.unintegrated")
DimPlot(H_P5, group.by = "treatment_status", reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "HP5")] # Need to manually change for every patient 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)

DimPlot(H_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(H_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(H_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(H_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(H_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#1    8613       0       0       0       0
#2       0       0       0    2983       0
#3       0       0       0       0    2216
#4       0       0     244       0       0
#5       0     324       0       0       0

cluster1 <- "C1"
cluster2 <- "C2"
cluster3 <- "C2"
cluster4 <- "C1"
cluster5 <- "C2"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
H_P5 <- AddMetaData(H_P5, cluster_group2)
DimPlot(H_P5, group.by = "Clone", reduction = "umap.unintegrated")
table(H_P5$Clone, H_P5$treatment_status)
length(which(H_P5$Clone == "C1")) / length(H_P5$Clone) # 0.6159249
length(which(H_P5$Clone == "C2")) / length(H_P5$Clone) # 0.3840751
#length(which(H_P5$Clone == "C3")) / length(H_P5$Clone) # 0.06210976

FeaturePlot(H_P5, "MYC", reduction = "umap.unintegrated")
FeaturePlot(H_P5, "IBSP", reduction = "umap.unintegrated")
VlnPlot(H_P5, "MYC", group.by = "Clone")
VlnPlot(H_P5, "IBSP", group.by = "Clone")
Idents(H_P5) <- "Clone"
markers <- FindMarkers(H_P5, ident.1 = "C1", only.pos = F)
DimPlot(H_P5, group.by = "Clone", reduction = "umap.unintegrated", split.by = "treatment_status")

levels(H_P5$Clone)
H_P5$Clone <- Idents(H_P5) 
levels(H_P5$Clone)
#Idents(H_P5) <- H_P5$Clone
p1 <- CellDimPlot(
  srt = H_P5, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933", "#339999"), 
  raster = T, 
  pt.size = 3, 
  split.by = "treatment_status", label_insitu = T) + theme(legend.position = "none")
p1

pdf("HP5_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

table(H_P5$Clone, H_P5$treatment_status)
#CellDimPlot3D(srt = H_P5, group.by = "Clone")
class(p1)
saveRDS(H_P5, "H_P5_cloned.rds")

# remove object when done 
rm(H_P5)

# LP2----- 
L_P2 <- subset(OB_cancer, idents = c("L_BF_P2", "L_AF_P2"))
DimPlot(L_P2)

L_P2[["RNA"]]$scale.data <- NULL
L_P2[["RNA"]]$data <- NULL
L_P2[["RNA"]] <- split(L_P2[["RNA"]], f = L_P2$sample) 
L_P2[["RNA"]]
DefaultAssay(L_P2) <- "RNA"
L_P2 <- NormalizeData(L_P2) %>% FindVariableFeatures() %>% ScaleData()
#L_P2[["percent.ribo"]] <- PercentageFeatureSet(object = L_P2, pattern = "^RP[SL]")
#Variable_features <- L_P2@assays[["RNA"]]@meta.data[["var.features"]][!is.na(L_P2@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(L_P2)[grep("^MT_", rownames(L_P2), invert=F)] 
#RB_genes <- rownames(L_P2)[grep("^RP[SL]", rownames(L_P2), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#L_P2 <- ScaleData(L_P2, vars.to.regress = c("percent.mt", "percent.ribo"))
#L_P2 <- ScaleData(L_P2)
L_P2 <- RunPCA(L_P2)
L_P2 <- FindNeighbors(L_P2, dims = 1:50, reduction = "pca")
L_P2 <- RunUMAP(L_P2, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
L_P2 <- JoinLayers(L_P2)
DimPlot(L_P2, label = T, reduction = "umap.unintegrated")
DimPlot(L_P2, group.by = "treatment_status", reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "LP2")] # Need to manually change for every patient 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)

DimPlot(L_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(L_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(L_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(L_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(L_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#8DD3C7 #D1C2D2 #D8C965 #EB8E8B #FFED6F
#1       0       0       0     601       0
#2       0       0       0       0     916
#3       0       0     644       0       0
#4     257       0       0       0       0
#5       0     372       0       0       0

cluster1 <- "C1"
cluster2 <- "C1"
cluster3 <- "C1"
cluster4 <- "C1"
cluster5 <- "C1"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
L_P2 <- AddMetaData(L_P2, cluster_group2)
DimPlot(L_P2, group.by = "Clone", reduction = "umap.unintegrated")
table(L_P2$Clone, L_P2$treatment_status)
length(which(L_P2$Clone == "C1")) / length(L_P2$Clone) # 0.6159249
length(which(L_P2$Clone == "C2")) / length(L_P2$Clone) # 0.3840751
#length(which(L_P2$Clone == "C3")) / length(L_P2$Clone) # 0.06210976

FeaturePlot(L_P2, "MYC", reduction = "umap.unintegrated")
FeaturePlot(L_P2, "IBSP", reduction = "umap.unintegrated")
VlnPlot(L_P2, "MYC", group.by = "Clone")
VlnPlot(L_P2, "IBSP", group.by = "Clone")
Idents(L_P2) <- "Clone"
markers <- FindMarkers(L_P2, ident.1 = "C1", only.pos = F)
DimPlot(L_P2, group.by = "Clone", reduction = "umap.unintegrated", split.by = "treatment_status")

p1 <- CellDimPlot(
  srt = L_P2, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933", "#339999"), 
  raster = T, 
  pt.size = 3, 
  split.by = "treatment_status", label_insitu = T) + theme(legend.position = "none")
p1

pdf("LP2_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

table(L_P2$Clone, L_P2$treatment_status)
#CellDimPlot3D(srt = L_P2, group.by = "Clone")
class(p1)
saveRDS(L_P2, "L_P2_cloned.rds")

# remove object when done 
rm(L_P2)

# LP5----- 
L_P5 <- subset(OB_cancer, idents = c("L_BF_P5", "L_AF_P5")) # Need to manually change 
DimPlot(L_P5)

L_P5[["RNA"]]$scale.data <- NULL
L_P5[["RNA"]]$data <- NULL
L_P5[["RNA"]] <- split(L_P5[["RNA"]], f = L_P5$sample) 
L_P5[["RNA"]]
DefaultAssay(L_P5) <- "RNA"
L_P5 <- NormalizeData(L_P5) %>% FindVariableFeatures() %>% ScaleData()
#L_P5[["percent.ribo"]] <- PercentageFeatureSet(object = L_P5, pattern = "^RP[SL]")
#Variable_features <- L_P5@assays[["RNA"]]@meta.data[["var.features"]][!is.na(L_P5@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(L_P5)[grep("^MT_", rownames(L_P5), invert=F)] 
#RB_genes <- rownames(L_P5)[grep("^RP[SL]", rownames(L_P5), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#L_P5 <- ScaleData(L_P5, vars.to.regress = c("percent.mt", "percent.ribo"))
#L_P5 <- ScaleData(L_P5)
L_P5 <- RunPCA(L_P5)
L_P5 <- FindNeighbors(L_P5, dims = 1:50, reduction = "pca")
L_P5 <- RunUMAP(L_P5, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
L_P5 <- JoinLayers(L_P5)
DimPlot(L_P5, label = T, reduction = "umap.unintegrated")
DimPlot(L_P5, group.by = "treatment_status", reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "LP5")] # Need to manually change for every patient 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)

DimPlot(L_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(L_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(L_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(L_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(L_P5, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#8DD3C7 #D1C2D2 #D8C965 #EB8E8B #FFED6F
#1       0       0       0       0     937
#2     839       0       0       0       0
#3       0     895       0       0       0
#4       0       0    1395       0       0
#5       0       0       0     643       0

cluster1 <- "C2"
cluster2 <- "C1"
cluster3 <- "C2"
cluster4 <- "C2"
cluster5 <- "C1"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
L_P5 <- AddMetaData(L_P5, cluster_group2)
DimPlot(L_P5, group.by = "Clone", reduction = "umap.unintegrated")
table(L_P5$Clone, L_P5$treatment_status)
length(which(L_P5$Clone == "C1")) / length(L_P5$Clone) # 0.3147165
length(which(L_P5$Clone == "C2")) / length(L_P5$Clone) # 0.6852835
#length(which(L_P5$Clone == "C3")) / length(L_P5$Clone) # 0.06210976

FeaturePlot(L_P5, "MYC", reduction = "umap.unintegrated")
FeaturePlot(L_P5, "IBSP", reduction = "umap.unintegrated")
VlnPlot(L_P5, "MYC", group.by = "Clone")
VlnPlot(L_P5, "IBSP", group.by = "Clone")
Idents(L_P5) <- "Clone"
markers <- FindMarkers(L_P5, ident.1 = "C1", only.pos = F)
DimPlot(L_P5, group.by = "Clone", reduction = "umap.unintegrated", split.by = "treatment_status")
levels(L_P5) <- c("C1", "C2")
L_P5$Clone <- Idents(L_P5)

p1 <- CellDimPlot(
  srt = L_P5, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933", "#339999"), 
  raster = T, 
  pt.size = 3, 
  split.by = "treatment_status", label_insitu = T) + theme(legend.position = "none")
p1

pdf("LP5_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

table(L_P5$Clone, L_P5$treatment_status)
#CellDimPlot3D(srt = L_P5, group.by = "Clone")
class(p1)
saveRDS(L_P5, "L_P5_cloned.rds")

# remove object when done 
rm(L_P5)


# LP6 ----- 
L_P6 <- subset(OB_cancer, idents = c("L_BF_P6", "L_AF_P6")) # Need to manually change 
DimPlot(L_P6)

L_P6[["RNA"]]$scale.data <- NULL
L_P6[["RNA"]]$data <- NULL
L_P6[["RNA"]] <- split(L_P6[["RNA"]], f = L_P6$sample) 
L_P6[["RNA"]]
DefaultAssay(L_P6) <- "RNA"
L_P6 <- NormalizeData(L_P6) %>% FindVariableFeatures() %>% ScaleData()
#L_P6[["percent.ribo"]] <- PercentageFeatureSet(object = L_P6, pattern = "^RP[SL]")
#Variable_features <- L_P6@assays[["RNA"]]@meta.data[["var.features"]][!is.na(L_P6@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(L_P6)[grep("^MT_", rownames(L_P6), invert=F)] 
#RB_genes <- rownames(L_P6)[grep("^RP[SL]", rownames(L_P6), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#L_P6 <- ScaleData(L_P6, vars.to.regress = c("percent.mt", "percent.ribo"))
#L_P6 <- ScaleData(L_P6)
L_P6 <- RunPCA(L_P6)
L_P6 <- FindNeighbors(L_P6, dims = 1:50, reduction = "pca")
L_P6 <- RunUMAP(L_P6, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
L_P6 <- JoinLayers(L_P6)
DimPlot(L_P6, label = T, reduction = "umap.unintegrated")
DimPlot(L_P6, group.by = "treatment_status", reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "LP6")] # Need to manually change for every patient 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)

DimPlot(L_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(L_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(L_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(L_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(L_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#1       0       0     730       0       0
#2       0       0       0     817       0
#3    2384       0       0       0       0
#4       0       0       0       0    3389
#5       0     291       0       0       0

cluster1 <- "C1"
cluster2 <- "C1"
cluster3 <- "C1"
cluster4 <- "C1"
cluster5 <- "C1"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
L_P6 <- AddMetaData(L_P6, cluster_group2)
DimPlot(L_P6, group.by = "Clone", reduction = "umap.unintegrated")
table(L_P6$Clone, L_P6$treatment_status)
length(which(L_P6$Clone == "C1")) / length(L_P6$Clone) # 0.3147165
length(which(L_P6$Clone == "C2")) / length(L_P6$Clone) # 0.6852835
#length(which(L_P6$Clone == "C3")) / length(L_P6$Clone) # 0.06210976

FeaturePlot(L_P6, "MYC", reduction = "umap.unintegrated")
FeaturePlot(L_P6, "IBSP", reduction = "umap.unintegrated")
VlnPlot(L_P6, "MYC", group.by = "Clone")
VlnPlot(L_P6, "IBSP", group.by = "Clone")
Idents(L_P6) <- "Clone"
markers <- FindMarkers(L_P6, ident.1 = "C1", only.pos = F)
DimPlot(L_P6, group.by = "Clone", reduction = "umap.unintegrated", split.by = "treatment_status")

p1 <- CellDimPlot(
  srt = L_P6, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933", "#339999"), 
  raster = T, 
  pt.size = 3, 
  split.by = "treatment_status", label_insitu = T) + theme(legend.position = "none")
p1

pdf("LP6_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()


table(L_P6$Clone, L_P6$treatment_status)
#CellDimPlot3D(srt = L_P6, group.by = "Clone")
class(p1)
saveRDS(L_P6, "L_P6_cloned.rds")

# remove object when done 
rm(L_P6)

# Baseline (Pre-treatment) only -----
# HP2 ----- 
H_P2 <- subset(OB_cancer, idents = c("H_BF_P2"))
DimPlot(H_P2)

H_P2[["RNA"]]$scale.data <- NULL
H_P2[["RNA"]]$data <- NULL
H_P2[["RNA"]] <- split(H_P2[["RNA"]], f = H_P2$sample) 
H_P2[["RNA"]]
DefaultAssay(H_P2) <- "RNA"
H_P2 <- NormalizeData(H_P2) %>% FindVariableFeatures() %>% ScaleData()
#H_P2[["percent.ribo"]] <- PercentageFeatureSet(object = H_P2, pattern = "^RP[SL]")
#Variable_features <- H_P2@assays[["RNA"]]@meta.data[["var.features"]][!is.na(H_P2@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(H_P2)[grep("^MT_", rownames(H_P2), invert=F)] 
#RB_genes <- rownames(H_P2)[grep("^RP[SL]", rownames(H_P2), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#H_P2 <- ScaleData(H_P2, vars.to.regress = c("percent.mt", "percent.ribo"))
#H_P2 <- ScaleData(H_P2)
H_P2 <- RunPCA(H_P2)
H_P2 <- FindNeighbors(H_P2, dims = 1:50, reduction = "pca")
H_P2 <- RunUMAP(H_P2, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(H_P2, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "HP2")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
#cnv_table1_m <- as.matrix(cnv_table1)
#cnv_table2_m <- as.matrix(cnv_table2)
#cnv_table3_m <- as.matrix(cnv_table3)
#cnv_table4_m <- as.matrix(cnv_table4)
#cnv_table5_m <- as.matrix(cnv_table5)

DimPlot(H_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(H_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(H_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(H_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(H_P2, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#1       0     912       0       0       0
#2       0       0       0       0     447
#3       0       0       0     376       0
#4     900       0       0       0       0
#5       0       0     828       0       0

cluster1 <- "C1"
cluster2 <- "C1"
cluster3 <- "C2"
cluster4 <- "C1"
cluster5 <- "C1"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
H_P2 <- AddMetaData(H_P2, cluster_group2)
DimPlot(H_P2, group.by = "Clone", reduction = "umap.unintegrated")
table(H_P2$Clone, H_P2$treatment_status)
length(which(H_P2$Clone == "C1")) / length(H_P2$Clone) # 0.7580372
length(which(H_P2$Clone == "C2")) / length(H_P2$Clone) # 0.08347434
length(which(H_P2$Clone == "C3")) / length(H_P2$Clone) # 0.1584884
saveRDS(H_P2, "H_P2_cloned.rds")
# H_P2 <- readRDS("H_P2_cloned.rds")

p1 <- CellDimPlot(
  srt = H_P2, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  #split.by = "treatment_status", 
  label_insitu = T) #+ 
#theme(legend.position = "none")
p1

pdf("HP1_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

#CellDimPlot3D(srt = H_P2, group.by = "Clone")
# remove object when done 
rm(H_P2)

# HP4 ----- 
H_P4 <- subset(OB_cancer, idents = c("H_BF_P4"))
DimPlot(H_P4)

H_P4[["RNA"]]$scale.data <- NULL
H_P4[["RNA"]]$data <- NULL
H_P4[["RNA"]] <- split(H_P4[["RNA"]], f = H_P4$sample) 
H_P4[["RNA"]]
DefaultAssay(H_P4) <- "RNA"
H_P4 <- NormalizeData(H_P4) %>% FindVariableFeatures() %>% ScaleData()
#H_P4[["percent.ribo"]] <- PercentageFeatureSet(object = H_P4, pattern = "^RP[SL]")
#Variable_features <- H_P4@assays[["RNA"]]@meta.data[["var.features"]][!is.na(H_P4@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(H_P4)[grep("^MT_", rownames(H_P4), invert=F)] 
#RB_genes <- rownames(H_P4)[grep("^RP[SL]", rownames(H_P4), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#H_P4 <- ScaleData(H_P4, vars.to.regress = c("percent.mt", "percent.ribo"))
#H_P4 <- ScaleData(H_P4)
H_P4 <- RunPCA(H_P4)
H_P4 <- FindNeighbors(H_P4, dims = 1:50, reduction = "pca")
H_P4 <- RunUMAP(H_P4, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(H_P4, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "HP4")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
#cnv_table1_m <- as.matrix(cnv_table1)
#cnv_table2_m <- as.matrix(cnv_table2)
#cnv_table3_m <- as.matrix(cnv_table3)
#cnv_table4_m <- as.matrix(cnv_table4)
#cnv_table5_m <- as.matrix(cnv_table5)

DimPlot(H_P4, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(H_P4, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(H_P4, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(H_P4, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(H_P4, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#8DD3C7 #D1C2D2 #D8C965 #EB8E8B #FFED6F
#1       0       0     499       0       0
#2       0       0       0     295       0
#3       0       0       0       0     392
#4     430       0       0       0       0
#5       0     108       0       0       0

cluster1 <- "C2"
cluster2 <- "C3"
cluster3 <- "C3"
cluster4 <- "C3"
cluster5 <- "C1"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
H_P4 <- AddMetaData(H_P4, cluster_group2)
DimPlot(H_P4, group.by = "Clone", reduction = "umap.unintegrated")
table(H_P4$Clone, H_P4$treatment_status)
length(which(H_P4$Clone == "C1")) / length(H_P4$Clone) # 0.7580372
length(which(H_P4$Clone == "C2")) / length(H_P4$Clone) # 0.08347434
length(which(H_P4$Clone == "C3")) / length(H_P4$Clone) # 0.1584884
saveRDS(H_P4, "H_P4_cloned.rds")
# H_P4 <- readRDS("H_P4_cloned.rds")

p1 <- CellDimPlot(
  srt = H_P6, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  #split.by = "treatment_status", 
  label_insitu = T) #+ 
#theme(legend.position = "none")
p1

pdf("HP1_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

#CellDimPlot3D(srt = H_P4, group.by = "Clone")
# remove object when done 
rm(H_P4)

# HP6 ----- 
H_P6 <- subset(OB_cancer, idents = c("H_BF_P6"))
DimPlot(H_P6)

H_P6[["RNA"]]$scale.data <- NULL
H_P6[["RNA"]]$data <- NULL
H_P6[["RNA"]] <- split(H_P6[["RNA"]], f = H_P6$sample) 
H_P6[["RNA"]]
DefaultAssay(H_P6) <- "RNA"
H_P6 <- NormalizeData(H_P6) %>% FindVariableFeatures() %>% ScaleData()
#H_P6[["percent.ribo"]] <- PercentageFeatureSet(object = H_P6, pattern = "^RP[SL]")
#Variable_features <- H_P6@assays[["RNA"]]@meta.data[["var.features"]][!is.na(H_P6@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(H_P6)[grep("^MT_", rownames(H_P6), invert=F)] 
#RB_genes <- rownames(H_P6)[grep("^RP[SL]", rownames(H_P6), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#H_P6 <- ScaleData(H_P6, vars.to.regress = c("percent.mt", "percent.ribo"))
#H_P6 <- ScaleData(H_P6)
H_P6 <- RunPCA(H_P6)
H_P6 <- FindNeighbors(H_P6, dims = 1:50, reduction = "pca")
H_P6 <- RunUMAP(H_P6, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(H_P6, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "HP6")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
#cnv_table1_m <- as.matrix(cnv_table1)
#cnv_table2_m <- as.matrix(cnv_table2)
#cnv_table3_m <- as.matrix(cnv_table3)
#cnv_table4_m <- as.matrix(cnv_table4)
#cnv_table5_m <- as.matrix(cnv_table5)

DimPlot(H_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(H_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(H_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(H_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(H_P6, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
##8DD3C7 #D1C2D2 #D8C965 #EB8E8B #FFED6F
#1       0       0       0    1374       0
#2       0       0    1385       0       0
#3       0    1224       0       0       0
#4     772       0       0       0       0
#5       0       0       0       0    1408

cluster1 <- "C2"
cluster2 <- "C2"
cluster3 <- "C2"
cluster4 <- "C1"
cluster5 <- "C2"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
H_P6 <- AddMetaData(H_P6, cluster_group2)
DimPlot(H_P6, group.by = "Clone", reduction = "umap.unintegrated")
table(H_P6$Clone, H_P6$treatment_status)
length(which(H_P6$Clone == "C1")) / length(H_P6$Clone) # 0.1252637
length(which(H_P6$Clone == "C2")) / length(H_P6$Clone) # 0.8747363
#length(which(H_P6$Clone == "C3")) / length(H_P6$Clone) # 0.1584884
saveRDS(H_P6, "H_P6_cloned.rds")
# H_P6 <- readRDS("H_P6_cloned.rds")

p1 <- CellDimPlot(
  srt = H_P6, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  #split.by = "treatment_status", 
  label_insitu = T) #+ 
#theme(legend.position = "none")
p1

pdf("HP1_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

#CellDimPlot3D(srt = H_P6, group.by = "Clone")
# remove object when done 
rm(H_P6)

# HP7 ----- 
H_P7 <- subset(OB_cancer, idents = c("H_BF_P7"))
DimPlot(H_P7)

H_P7[["RNA"]]$scale.data <- NULL
H_P7[["RNA"]]$data <- NULL
H_P7[["RNA"]] <- split(H_P7[["RNA"]], f = H_P7$sample) 
H_P7[["RNA"]]
DefaultAssay(H_P7) <- "RNA"
H_P7 <- NormalizeData(H_P7) %>% FindVariableFeatures() %>% ScaleData()
#H_P7[["percent.ribo"]] <- PercentageFeatureSet(object = H_P7, pattern = "^RP[SL]")
#Variable_features <- H_P7@assays[["RNA"]]@meta.data[["var.features"]][!is.na(H_P7@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(H_P7)[grep("^MT_", rownames(H_P7), invert=F)] 
#RB_genes <- rownames(H_P7)[grep("^RP[SL]", rownames(H_P7), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#H_P7 <- ScaleData(H_P7, vars.to.regress = c("percent.mt", "percent.ribo"))
#H_P7 <- ScaleData(H_P7)
H_P7 <- RunPCA(H_P7)
H_P7 <- FindNeighbors(H_P7, dims = 1:50, reduction = "pca")
H_P7 <- RunUMAP(H_P7, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(H_P7, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "HP7")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
#cnv_table1_m <- as.matrix(cnv_table1)
#cnv_table2_m <- as.matrix(cnv_table2)
#cnv_table3_m <- as.matrix(cnv_table3)
#cnv_table4_m <- as.matrix(cnv_table4)
#cnv_table5_m <- as.matrix(cnv_table5)

DimPlot(H_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(H_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(H_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(H_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(H_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
##8DD3C7 #D1C2D2 #D8C965 #EB8E8B #FFED6F
#1       0       0     342       0       0
#2       0       0       0       0      74
#3       0       0       0     250       0
#4       0     117       0       0       0
#5      60       0       0       0       0

cluster1 <- "C1"
cluster2 <- "C2"
cluster3 <- "C1"
cluster4 <- "C3"
cluster5 <- "C1"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
H_P7 <- AddMetaData(H_P7, cluster_group2)
DimPlot(H_P7, group.by = "Clone", reduction = "umap.unintegrated")
table(H_P7$Clone, H_P7$treatment_status)
length(which(H_P7$Clone == "C1")) / length(H_P7$Clone) # 0.1252637
length(which(H_P7$Clone == "C2")) / length(H_P7$Clone) # 0.8747363
length(which(H_P7$Clone == "C3")) / length(H_P7$Clone) # 0.1584884
saveRDS(H_P7, "H_P7_cloned.rds")
# H_P7 <- readRDS("H_P7_cloned.rds")

p1 <- CellDimPlot(
  srt = H_P7, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  #split.by = "treatment_status", 
  label_insitu = T) #+ 
#theme(legend.position = "none")
p1

pdf("HP1_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

#CellDimPlot3D(srt = H_P7, group.by = "Clone")
# remove object when done 
rm(H_P7)

# HP8 ----- 
H_P8 <- subset(OB_cancer, idents = c("H_BF_P8"))
DimPlot(H_P8)

H_P8[["RNA"]]$scale.data <- NULL
H_P8[["RNA"]]$data <- NULL
H_P8[["RNA"]] <- split(H_P8[["RNA"]], f = H_P8$sample) 
H_P8[["RNA"]]
DefaultAssay(H_P8) <- "RNA"
H_P8 <- NormalizeData(H_P8) %>% FindVariableFeatures() %>% ScaleData()
#H_P8[["percent.ribo"]] <- PercentageFeatureSet(object = H_P8, pattern = "^RP[SL]")
#Variable_features <- H_P8@assays[["RNA"]]@meta.data[["var.features"]][!is.na(H_P8@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(H_P8)[grep("^MT_", rownames(H_P8), invert=F)] 
#RB_genes <- rownames(H_P8)[grep("^RP[SL]", rownames(H_P8), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#H_P8 <- ScaleData(H_P8, vars.to.regress = c("percent.mt", "percent.ribo"))
#H_P8 <- ScaleData(H_P8)
H_P8 <- RunPCA(H_P8)
H_P8 <- FindNeighbors(H_P8, dims = 1:50, reduction = "pca")
H_P8 <- RunUMAP(H_P8, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(H_P8, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "HP8")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
#cnv_table1_m <- as.matrix(cnv_table1)
#cnv_table2_m <- as.matrix(cnv_table2)
#cnv_table3_m <- as.matrix(cnv_table3)
#cnv_table4_m <- as.matrix(cnv_table4)
#cnv_table5_m <- as.matrix(cnv_table5)

DimPlot(H_P8, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(H_P8, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(H_P8, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(H_P8, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(H_P8, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
##8DD3C7 #D1C2D2 #D8C965 #EB8E8B #FFED6F
#1       0     649       0       0       0
#2       0       0       0     342       0
#3       0       0       0       0     781
#4     205       0       0       0       0
#5       0       0     558       0       0

cluster1 <- "C1"
cluster2 <- "C3"
cluster3 <- "C2"
cluster4 <- "C4"
cluster5 <- "C2"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
H_P8 <- AddMetaData(H_P8, cluster_group2)
DimPlot(H_P8, group.by = "Clone", reduction = "umap.unintegrated")
table(H_P8$Clone, H_P8$treatment_status)
length(which(H_P8$Clone == "C1")) / length(H_P8$Clone) # 0.2560158
length(which(H_P8$Clone == "C2")) / length(H_P8$Clone) # 0.5282051
length(which(H_P8$Clone == "C3")) / length(H_P8$Clone) # 0.1349112
length(which(H_P8$Clone == "C4")) / length(H_P8$Clone) # 0.08086785
saveRDS(H_P8, "H_P8_cloned.rds")
# H_P8 <- readRDS("H_P8_cloned.rds")

p1 <- CellDimPlot(
  srt = H_P8, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  #split.by = "treatment_status", 
  label_insitu = T) #+ 
#theme(legend.position = "none")
p1

pdf("HP8_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

#CellDimPlot3D(srt = H_P8, group.by = "Clone")
# remove object when done 
rm(H_P8)

# LP1 ----- 
L_P1 <- subset(OB_cancer, idents = c("L_BF_P1"))
DimPlot(L_P1)

L_P1[["RNA"]]$scale.data <- NULL
L_P1[["RNA"]]$data <- NULL
L_P1[["RNA"]] <- split(L_P1[["RNA"]], f = L_P1$sample) 
L_P1[["RNA"]]
DefaultAssay(L_P1) <- "RNA"
L_P1 <- NormalizeData(L_P1) %>% FindVariableFeatures() %>% ScaleData()
#L_P1[["percent.ribo"]] <- PercentageFeatureSet(object = L_P1, pattern = "^RP[SL]")
#Variable_features <- L_P1@assays[["RNA"]]@meta.data[["var.features"]][!is.na(L_P1@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(L_P1)[grep("^MT_", rownames(L_P1), invert=F)] 
#RB_genes <- rownames(L_P1)[grep("^RP[SL]", rownames(L_P1), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#L_P1 <- ScaleData(L_P1, vars.to.regress = c("percent.mt", "percent.ribo"))
#L_P1 <- ScaleData(L_P1)
L_P1 <- RunPCA(L_P1)
L_P1 <- FindNeighbors(L_P1, dims = 1:50, reduction = "pca")
L_P1 <- RunUMAP(L_P1, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(L_P1, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "LP1")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
#cnv_table1_m <- as.matrix(cnv_table1)
#cnv_table2_m <- as.matrix(cnv_table2)
#cnv_table3_m <- as.matrix(cnv_table3)
#cnv_table4_m <- as.matrix(cnv_table4)
#cnv_table5_m <- as.matrix(cnv_table5)

DimPlot(L_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(L_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(L_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(L_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(L_P1, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#8DD3C7 #D1C2D2 #D8C965 #EB8E8B #FFED6F
#1       0       0       0       0     550
#2     300       0       0       0       0
#3       0     318       0       0       0
#4       0       0       0     241       0
#5       0       0     279       0       0
cluster1 <- "C2"
cluster2 <- "C1"
cluster3 <- "C2"
cluster4 <- "C2"
cluster5 <- "C2"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
L_P1 <- AddMetaData(L_P1, cluster_group2)
DimPlot(L_P1, group.by = "Clone", reduction = "umap.unintegrated")
table(L_P1$Clone, L_P1$treatment_status)
length(which(L_P1$Clone == "C1")) / length(L_P1$Clone) # 0.2560158
length(which(L_P1$Clone == "C2")) / length(L_P1$Clone) # 0.5282051
#length(which(L_P1$Clone == "C3")) / length(L_P1$Clone) # 0.1349112
#length(which(L_P1$Clone == "C4")) / length(L_P1$Clone) # 0.08086785
saveRDS(L_P1, "L_P1_cloned.rds")
# L_P1 <- readRDS("L_P1_cloned.rds")

p1 <- CellDimPlot(
  srt = L_P1, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  #split.by = "treatment_status", 
  label_insitu = T) #+ 
#theme(legend.position = "none")
p1

pdf("LP1_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

#CellDimPlot3D(srt = L_P1, group.by = "Clone")
# remove object when done 
rm(L_P1)

# LP3 and LP4 truncated due to monoclonal -----

# LP7 ----- 
L_P7 <- subset(OB_cancer, idents = c("L_BF_P7"))
DimPlot(L_P7)

L_P7[["RNA"]]$scale.data <- NULL
L_P7[["RNA"]]$data <- NULL
L_P7[["RNA"]] <- split(L_P7[["RNA"]], f = L_P7$sample) 
L_P7[["RNA"]]
DefaultAssay(L_P7) <- "RNA"
L_P7 <- NormalizeData(L_P7) %>% FindVariableFeatures() %>% ScaleData()
#L_P7[["percent.ribo"]] <- PercentageFeatureSet(object = L_P7, pattern = "^RP[SL]")
#Variable_features <- L_P7@assays[["RNA"]]@meta.data[["var.features"]][!is.na(L_P7@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(L_P7)[grep("^MT_", rownames(L_P7), invert=F)] 
#RB_genes <- rownames(L_P7)[grep("^RP[SL]", rownames(L_P7), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#L_P7 <- ScaleData(L_P7, vars.to.regress = c("percent.mt", "percent.ribo"))
#L_P7 <- ScaleData(L_P7)
L_P7 <- RunPCA(L_P7)
L_P7 <- FindNeighbors(L_P7, dims = 1:50, reduction = "pca")
L_P7 <- RunUMAP(L_P7, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(L_P7, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- 
patient <- patient_names[which(patient_names == "LP7")] 
dir <- paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
#cnv_table1_m <- as.matrix(cnv_table1)
#cnv_table2_m <- as.matrix(cnv_table2)
#cnv_table3_m <- as.matrix(cnv_table3)
#cnv_table4_m <- as.matrix(cnv_table4)
#cnv_table5_m <- as.matrix(cnv_table5)

DimPlot(L_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(L_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(L_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(L_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(L_P7, cells.highlight = WhichCells(OB_cancer, rownames(cnv_table5)), reduction = "umap.unintegrated")

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#8DD3C7 #D1C2D2 #D8C965 #EB8E8B #FFED6F
#1       0       0       0    2118       0
#2       0     480       0       0       0
#3       0       0     901       0       0
#4     685       0       0       0       0
#5       0       0       0       0     758

cluster1 <- "C3"
cluster2 <- "C2"
cluster3 <- "C1"
cluster4 <- "C4"
cluster5 <- "C2"

colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
table(cluster_group2$cell_type)
colnames(cluster_group2)[2] <- "Clone"

#write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#cluster_group3 <- cbind(cluster_group2, dd)
#cluster_group_all <- rbind(cluster_group_all, cluster_group3)
L_P7 <- AddMetaData(L_P7, cluster_group2)
DimPlot(L_P7, group.by = "Clone", reduction = "umap.unintegrated")
table(L_P7$Clone, L_P7$treatment_status)
length(which(L_P7$Clone == "C1")) / length(L_P7$Clone) # 0.2560158
length(which(L_P7$Clone == "C2")) / length(L_P7$Clone) # 0.5282051
length(which(L_P7$Clone == "C3")) / length(L_P7$Clone) # 0.1349112
length(which(L_P7$Clone == "C4")) / length(L_P7$Clone) # 0.08086785
saveRDS(L_P7, "L_P7_cloned.rds")
# L_P7 <- readRDS("L_P7_cloned.rds")

p1 <- CellDimPlot(
  srt = L_P7, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  #split.by = "treatment_status", 
  label_insitu = T) #+ 
#theme(legend.position = "none")
p1

pdf("LP7_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

#CellDimPlot3D(srt = L_P7, group.by = "Clone")
# remove object when done 
rm(L_P7)

