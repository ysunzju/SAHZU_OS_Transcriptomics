library(here) # project-root-relative paths; run scripts from repository root
##################### InferCNV for ST data (Related to Figures 2G,3J, Figure S3B) ##################### 
##################### InferCNV script ##################### 
rm(list = ls())
library(infercnv)
library(Seurat)

pat <- commandArgs(TRUE)[1]
numcores = 30

# Step1 : PrepareFiles
scObject1 <- readRDS( paste0("/data/s01020/InferCNV_for_ST_tumour_ident_251015/tOS.",pat,".scObject.rds"))
#Idents(scObject1) <- "cluster"
#scObject1 <- subset(scObject1, idents = c("Osteoblasts"))
scObject1$CellType <- "Osteoblasts"
Ref <- readRDS("~/InferCNV_for_ST_tumour_ident_251015/ref.rds")
Ref$CellType <- 'ImmuneCells'
scObject <- merge(x=scObject1, y=Ref)
scObject$anno <- as.character(scObject$CellType)
scObject$anno[scObject$anno %in% c("Osteoblasts")] <- "Tumor Candidate"
write.table(scObject@meta.data[,"anno",drop=FALSE], paste0(pat,"_annotations_file.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t" )
#table(scObject$CellType)
#table(scObject$anno)
ref_group_names <- setdiff(names(table(as.character(scObject$anno))),"Tumor Candidate")
ref_group_names

# Step2 : CreateInfercnvObject
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix= GetAssayData(scObject, layer="counts"),
  annotations_file= paste0(pat,"_annotations_file.txt"),
  delim="\t",
  #gene_order_file= here("data/InferCNV_31_01_25/infercnv_gene_order_nodup.txt"),
  gene_order_file= "~/InferCNV_for_ST_tumour_ident_251015/hg38_gencode_v27.txt", 
  #gene_order_file= here("data/InferCNV_31_01_25/infercnv_gene_order.txt"),
  #gene_order_file= here("data/InferCNV_31_01_25/gencode_v21_gen_pos.complete.txt"), 
  ref_group_names= ref_group_names,
  max_cells_per_group = NULL,
  min_max_counts_per_cell = c(100, +Inf),
  chr_exclude = c("chrX", "chrY", "chrM")
)

new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
  new_gene_order = rbind(new_gene_order, 
                         infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]]
                                                       == chr_name) , , drop=FALSE])
}

names(new_gene_order) <- c("chr", "start", "stop") 
infercnv_obj@gene_order = new_gene_order 
infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]

# Step3 : infercnv
infercnv_obj = infercnv::run( 
  infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10X Genomics
  out_dir= paste0("~/InferCNV_for_ST_tumour_ident_251015/", pat), 
  cluster_by_groups=FALSE, # If observations are defined according to groups (ie. patients), each group of cells will be clustered separately.
  denoise=TRUE,
  HMM=TRUE, # when set to True, runs HMM to predict CNV level (default: FALSE)
  analysis_mode="subclusters", # options(samples|subclusters|cells), Grouping level for image filtering or HMM predictions. default: samples (fastest, but subclusters is ideal)
  tumor_subcluster_partition_method="random_trees",
  scale_data = TRUE,
  leiden_resolution = 0.01,  
  k_obs_groups = 5,     
  output_format="pdf",  
  num_threads = numcores, 
  write_expr_matrix = T
) 

##################### Inferring clonal structure, UMAP, spatial distribution plot, and CNV heatmap (related to Fig 2h, extended data fig 3b) ##################### 
rm(list = ls())
# disabled, run from repo root: setwd(here("data/New_ST_251011/Output"))
#setwd(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny"))
library(Seurat)
library(SCP)
library(ComplexHeatmap)
library(circlize)

# Load data ----- 
ST_merge <- readRDS(here("data/New_ST_251011/Output/ST_merge_251015.rds"))
ref <- readRDS(here("data/InferCNV_rerun_for_phylo_22_06_25/ref.rds"))

Idents(ST_merge) <- "orig.ident"
DefaultAssay(ST_merge) # "SCT"

# SAHST001 ----- 
SAHST001 <- subset(ST_merge, idents = c("SAHST001"))
DimPlot(SAHST001)
DefaultAssay(SAHST001) <- "Spatial"
SAHST001@assays$SCT <- NULL
SAHST001 <- SCTransform(SAHST001, assay = "Spatial", verbose = FALSE)

SAHST001 <- RunPCA(SAHST001)
SAHST001 <- FindNeighbors(SAHST001, dims = 1:50, reduction = "pca")
SAHST001 <- RunUMAP(SAHST001, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(SAHST001, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- ]
patient <- levels(ST_merge)[which(levels(ST_merge) == "SAHST001")] 
dir <- paste0(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/"), patient)
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

DimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table5)), reduction = "umap.unintegrated")

SpatialDimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table1)), image.alpha = 0.5)
SpatialDimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table2)), image.alpha = 0.5)
SpatialDimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table3)), image.alpha = 0.5)
SpatialDimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table4)), image.alpha = 0.5)
SpatialDimPlot(SAHST001, cells.highlight = WhichCells(SAHST001, rownames(cnv_table5)), image.alpha = 0.5)

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#1       0     178       0       0       0
#2       0       0       0       0     689
#3       0       0     536       0       0
#4       0       0       0     213       0
#5     119       0       0       0       0

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
SAHST001 <- AddMetaData(SAHST001, cluster_group2)
DimPlot(SAHST001, group.by = "Clone", reduction = "umap.unintegrated")
SpatialDimPlot(SAHST001, group.by = "Clone")
#table(SAHST001$Clone, SAHST001$treatment_status)
length(which(SAHST001$Clone == "C1")) / length(SAHST001$Clone) # 0.7580372
length(which(SAHST001$Clone == "C2")) / length(SAHST001$Clone) # 0.08347434
length(which(SAHST001$Clone == "C3")) / length(SAHST001$Clone) # 0.1584884
saveRDS(SAHST001, "SAHST001_cloned.rds")
# SAHST001 <- readRDS("SAHST001_cloned.rds")

Idents(SAHST001) <- "Clone"
levels(SAHST001) <- c("C1")
SAHST001$Clone <- Idents(SAHST001)
p1 <- CellDimPlot(
  srt = SAHST001, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933", "#339999"), 
  raster = T, 
  pt.size = 3, 
  #split.by = "treatment_status", 
  label_insitu = T) + theme(legend.position = "none")
p1

pdf("SAHST001_Clonal_DimPlot.pdf", height = 3, width = 3)
p1
dev.off()

### Heatmap
#SAHST001 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P5_cloned.rds"))
expr <- read.delim(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/SAHST001/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(SAHST001)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(SAHST001$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
    colnames(expr)[which(colnames(expr) %in% obs_cell)])]))

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df = kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

annotation_row = data.frame(
  #k_cluster =kmeans_df$k_cluster,
  group = kmeans_df$group
)
row.names(annotation_row) <- kmeans_df$cell_id
#color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
#names(color_cluster)=as.character(1:5)
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(group=c("C1" = "#9999cc",
                                            "C2" = "#669933", 
                                            "C3" = "#339999"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_SAHST001.pdf",width = 6, height = 4)
ht = Heatmap(t(log2(expr))[rownames(annotation_row),],
             col = colorRamp2(c(-0.5,0,0.5), c("#375093","white","#831A21")),
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             show_row_names = F,
             column_split = factor(gene_pos$V2, new_cluster),
             heatmap_legend_param = list(title = "inferCNV",
                                         direction = "vertical",
                                         title_position = "leftcenter-rot",
                                         legend_height = unit(3, "cm")),
             left_annotation = left_anno, 
             row_title = NULL,
             column_title = NULL,
             top_annotation = top_color,
             use_raster = T, 
             border = T)
draw(ht, heatmap_legend_side = "right")
dev.off()

#CellDimPlot3D(srt = SAHST001, group.by = "Clone")
# remove object when done 

### Plot the clones
#color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
#Idents(ST_merge) <- "orig.ident"
#SAHST005 <- subset(ST_merge, idents = "SAHST002")
SAHST001@reductions$spatial = SAHST001@reductions$umap
SAHST001@reductions$spatial@key = 'spatial_'
SAHST001@reductions$spatial@cell.embeddings = as.matrix(SAHST001@images$image@coordinates[,c(3,2)])
SAHST001@reductions$spatial@cell.embeddings[,2] = -SAHST001@reductions$spatial@cell.embeddings[,2]
colnames(SAHST001@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
data_mat <- data.frame(SAHST001$Clone, #SAHST005$S_phase, SAHST002_plot$Os_min, 
                       SAHST001@reductions$spatial@cell.embeddings)

# Clone DimPlot 
# color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST001.Clone)) +
  geom_point(shape = 16, size = 1) +
  scale_color_manual(values = c("#9999cc", "#669933",# "#339999", 
                                "lightgrey")) + 
  #scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
  #                                  "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
  #                                  "#d73027", "#a50026"), 
  #colours = c("lightblue","lightyellow","red"),  
  #                      limits = color_scale_val) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 97/57) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf('SAHST001_Clone_Spatial_mapping.pdf', height = 3, width = 3)
print(p1)
dev.off()

rm(SAHST001)

# SAHST005 ----- 
SAHST005 <- subset(ST_merge, idents = c("SAHST005"))
DimPlot(SAHST005)
DefaultAssay(SAHST005) <- "Spatial"
SAHST005@assays$SCT <- NULL
SAHST005 <- SCTransform(SAHST005, assay = "Spatial", verbose = FALSE)

SAHST005 <- RunPCA(SAHST005)
SAHST005 <- FindNeighbors(SAHST005, dims = 1:50, reduction = "pca")
SAHST005 <- RunUMAP(SAHST005, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(SAHST005, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- ]
patient <- levels(ST_merge)[which(levels(ST_merge) == "SAHST005")] 
dir <- paste0(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/"), patient)
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

DimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table5)), reduction = "umap.unintegrated")

SpatialDimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table1)))
SpatialDimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table2)))
SpatialDimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table3)))
SpatialDimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table4)))
SpatialDimPlot(SAHST005, cells.highlight = WhichCells(SAHST005, rownames(cnv_table5)))

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#1       0       0       0     305       0
#2       0     504       0       0       0
#3       0       0     446       0       0
#4       0       0       0       0     581
#5     124       0       0       0       0

cluster1 <- "Non-tumour-dominant"
cluster2 <- "C1"
cluster3 <- "C2"
cluster4 <- "C1"
cluster5 <- "Non-tumour-dominant"

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
SAHST005 <- AddMetaData(SAHST005, cluster_group2)
DimPlot(SAHST005, group.by = "Clone", reduction = "umap.unintegrated")
SpatialDimPlot(SAHST005, group.by = "Clone")
#table(SAHST005$Clone, SAHST005$treatment_status)
length(which(SAHST005$Clone == "C1")) / length(SAHST005$Clone) # 0.7580372
length(which(SAHST005$Clone == "C2")) / length(SAHST005$Clone) # 0.08347434
length(which(SAHST005$Clone == "C3")) / length(SAHST005$Clone) # 0.1584884
length(which(SAHST005$Clone == "Non-tumour-dominant")) / length(SAHST005$Clone) # 0.1584884
SAHST005$Clone[SAHST005$Clone == "SAHST005"] <- "Non-tumour-dominant"

saveRDS(SAHST005, "SAHST005_cloned.rds")
# SAHST005 <- readRDS("SAHST005_cloned.rds")

Idents(SAHST005) <- "Clone"
levels(SAHST005) <- c("C1", "C2", #"C3", 
                 "Non-tumour-dominant")
SAHST005$Clone <- Idents(SAHST005)
p1 <- CellDimPlot(
  srt = SAHST005, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933",# "#339999", 
               "lightgrey"), 
  raster = T, 
  pt.size = 3, 
  #split.by = "treatment_status", 
  label_insitu = T) + theme(legend.position = "none")
p1

pdf("SAHST005_Clonal_DimPlot.pdf", height = 3, width = 6)
p1
dev.off()

### Plot the clones
#color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
#Idents(ST_merge) <- "orig.ident"
#SAHST005 <- subset(ST_merge, idents = "SAHST002")
SAHST005@reductions$spatial = SAHST005@reductions$umap
SAHST005@reductions$spatial@key = 'spatial_'
SAHST005@reductions$spatial@cell.embeddings = as.matrix(SAHST005@images$image@coordinates[,c(3,2)])
SAHST005@reductions$spatial@cell.embeddings[,2] = -SAHST005@reductions$spatial@cell.embeddings[,2]
colnames(SAHST005@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
data_mat <- data.frame(SAHST005$Clone, #SAHST005$S_phase, SAHST002_plot$Os_min, 
                       SAHST005@reductions$spatial@cell.embeddings)

# Clone DimPlot 
# color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST005.Clone)) +
  geom_point(shape = 16, size = 1) +
  scale_color_manual(values = c("#9999cc", "#669933",# "#339999", 
                                   "lightgrey")) + 
  #scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
  #                                  "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
  #                                  "#d73027", "#a50026"), 
                        #colours = c("lightblue","lightyellow","red"),  
  #                      limits = color_scale_val) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 1.7) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf('SAHST005_Clone_Spatial_mapping.pdf', height = 3, width = 3)
print(p1)
dev.off()

# SAHST005 genes
SAHST005@reductions$spatial = SAHST005@reductions$umap
SAHST005@reductions$spatial@key = 'spatial_'
SAHST005@reductions$spatial@cell.embeddings = as.matrix(SAHST005@images$image@coordinates[,c(3,2)])
SAHST005@reductions$spatial@cell.embeddings[,2] = -SAHST005@reductions$spatial@cell.embeddings[,2]
colnames(SAHST005@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

### Plot
#color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
#data_mat <- data.frame(J, SAHST005_plot$S_phase, SAHST005_plot$Os_min ,SAHST005_plot@reductions$spatial@cell.embeddings)
SAHST005_exp <- as.data.frame(t(GetAssayData(object = SAHST005)[c("ERBB2", "TOP2A", "BIRC5", "MKI67"), ]))
SAHST005_exp$cell <- rownames(SAHST005_exp)

SAHST005_embed <- as.data.frame(SAHST005@reductions$spatial@cell.embeddings)
SAHST005_embed$cell <- rownames(SAHST005_embed)

mat_for_plot <- merge(SAHST005_embed, SAHST005_exp, by = "cell")

# ERBB2 
color_scale_val <- range(mat_for_plot$ERBB2)
color_scale_val <- range(mat_for_plot$TOP2A)
color_scale_val <- range(mat_for_plot$MKI67)
color_scale_val <- range(mat_for_plot$BIRC5)

p1 <- ggplot(mat_for_plot, aes(x = spatial_1, y = spatial_2, colour = TOP2A)) +
  geom_point(shape = 16, size = 1) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026"), 
                        #colours = c("lightblue","lightyellow","red"),  
                        limits = color_scale_val) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 1.7) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf('SAHST005_ECM.pdf', height = 3, width = 3)
print(p1)
dev.off()

#CellDimPlot3D(srt = SAHST005, group.by = "Clone")
# remove object when done 
rm(SAHST005)

# SAHST002 ----- 
SAHST002 <- subset(ST_merge, idents = c("SAHST002"))
DimPlot(SAHST002)
DefaultAssay(SAHST002) <- "Spatial"
SAHST002@assays$SCT <- NULL
SAHST002 <- SCTransform(SAHST002, assay = "Spatial", verbose = FALSE)

SAHST002 <- RunPCA(SAHST002)
SAHST002 <- FindNeighbors(SAHST002, dims = 1:50, reduction = "pca")
SAHST002 <- RunUMAP(SAHST002, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(SAHST002, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- ]
patient <- levels(ST_merge)[which(levels(ST_merge) == "SAHST002")] 
dir <- paste0(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)

DimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table5)), reduction = "umap.unintegrated")

SpatialDimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table1)), image.alpha = 0.5)
SpatialDimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table2)), image.alpha = 0.5)
SpatialDimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table3)), image.alpha = 0.5)
SpatialDimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table4)), image.alpha = 0.5)
SpatialDimPlot(SAHST002, cells.highlight = WhichCells(SAHST002, rownames(cnv_table5)), image.alpha = 0.5)

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#1       0       0       0     199       0
#2       0       0       0       0     269
#3       0       0     145       0       0
#4     147       0       0       0       0
#5       0     128       0       0       0

cluster1 <- "C1"
cluster2 <- "C2"
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
SAHST002 <- AddMetaData(SAHST002, cluster_group2)
DimPlot(SAHST002, group.by = "Clone", reduction = "umap.unintegrated")
SpatialDimPlot(SAHST002, group.by = "Clone")
#table(SAHST002$Clone, SAHST002$treatment_status)
length(which(SAHST002$Clone == "C1")) / length(SAHST002$Clone) # 0.7580372
length(which(SAHST002$Clone == "C2")) / length(SAHST002$Clone) # 0.08347434
length(which(SAHST002$Clone == "C3")) / length(SAHST002$Clone) # 0.1584884
saveRDS(SAHST002, "SAHST002_cloned.rds")
# SAHST002 <- readRDS("SAHST002_cloned.rds")

#Idents(SAHST002) <- "Clone"
SAHST002$Clone <- factor(SAHST002$Clone, levels = c("C1", "C2", "C3"))
#SAHST002$Clone <- Idents(SAHST002)
p1 <- CellDimPlot(
  srt = SAHST002, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = c("#9999cc", "#669933", "#339999"), 
  raster = T, 
  pt.size = 3, 
  #split.by = "treatment_status", 
  label_insitu = T) + theme(legend.position = "none")
p1

pdf("SAHST002_Clonal_DimPlot.pdf", height = 3, width = 3)
p1
dev.off()

### Heatmap
#SAHST002 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P5_cloned.rds"))
expr <- read.delim(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/SAHST002/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(SAHST002)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(SAHST002$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
    colnames(expr)[which(colnames(expr) %in% obs_cell)])]))

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df = kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

annotation_row = data.frame(
  #k_cluster =kmeans_df$k_cluster,
  group = kmeans_df$group
)
row.names(annotation_row) <- kmeans_df$cell_id
#color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
#names(color_cluster)=as.character(1:5)
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(group=c("C1" = "#9999cc",
                                            "C2" = "#669933", 
                                            "C3" = "#339999"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_SAHST002.pdf",width = 6, height = 4)
ht = Heatmap(t(log2(expr))[rownames(annotation_row),],
             col = colorRamp2(c(-0.5,0,0.5), c("#375093","white","#831A21")),
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             show_row_names = F,
             column_split = factor(gene_pos$V2, new_cluster),
             heatmap_legend_param = list(title = "inferCNV",
                                         direction = "vertical",
                                         title_position = "leftcenter-rot",
                                         legend_height = unit(3, "cm")),
             left_annotation = left_anno, 
             row_title = NULL,
             column_title = NULL,
             top_annotation = top_color,
             use_raster = T, 
             border = T)
draw(ht, heatmap_legend_side = "right")
dev.off()

#CellDimPlot3D(srt = SAHST002, group.by = "Clone")
# remove object when done 

### Plot the clones
#color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
#Idents(ST_merge) <- "orig.ident"
#SAHST005 <- subset(ST_merge, idents = "SAHST002")
SAHST002@reductions$spatial = SAHST002@reductions$umap
SAHST002@reductions$spatial@key = 'spatial_'
SAHST002@reductions$spatial@cell.embeddings = as.matrix(SAHST002@images$image@coordinates[,c(3,2)])
SAHST002@reductions$spatial@cell.embeddings[,2] = -SAHST002@reductions$spatial@cell.embeddings[,2]
colnames(SAHST002@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
data_mat <- data.frame(SAHST002$Clone, #SAHST005$S_phase, SAHST002_plot$Os_min, 
                       SAHST002@reductions$spatial@cell.embeddings)

# Clone DimPlot 
# color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST002.Clone)) +
  geom_point(shape = 16, size = 1) +
  scale_color_manual(values = c("#9999cc", "#669933",# "#339999", 
                                "#339999")) + 
  #scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
  #                                  "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
  #                                  "#d73027", "#a50026"), 
  #colours = c("lightblue","lightyellow","red"),  
  #                      limits = color_scale_val) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 1.7) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf('SAHST002_Clone_Spatial_mapping.pdf', height = 3, width = 3)
print(p1)
dev.off()

rm(SAHST002)

# SAHST004 ----- 
SAHST004 <- subset(ST_merge, idents = c("SAHST004"))
DimPlot(SAHST004)
DefaultAssay(SAHST004) <- "Spatial"
SAHST004@assays$SCT <- NULL
SAHST004 <- SCTransform(SAHST004, assay = "Spatial", verbose = FALSE)

SAHST004 <- RunPCA(SAHST004)
SAHST004 <- FindNeighbors(SAHST004, dims = 1:50, reduction = "pca")
SAHST004 <- RunUMAP(SAHST004, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(SAHST004, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- ]
patient <- levels(ST_merge)[which(levels(ST_merge) == "SAHST004")] 
dir <- paste0(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)

DimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table5)), reduction = "umap.unintegrated")

SpatialDimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table1)), image.alpha = 0.5)
SpatialDimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table2)), image.alpha = 0.5)
SpatialDimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table3)), image.alpha = 0.5)
SpatialDimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table4)), image.alpha = 0.5)
SpatialDimPlot(SAHST004, cells.highlight = WhichCells(SAHST004, rownames(cnv_table5)), image.alpha = 0.5)

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#1       0       0     563       0       0
#2       0       0       0    1207       0
#3       0     393       0       0       0
#4       0       0       0       0     442
#5     956       0       0       0       0

cluster1 <- "Peri.tumour"
cluster2 <- "C1"
cluster3 <- "Non.tumour"
cluster4 <- "Peri.tumour"
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
SAHST004 <- AddMetaData(SAHST004, cluster_group2)
DimPlot(SAHST004, group.by = "Clone", reduction = "umap.unintegrated")
SpatialDimPlot(SAHST004, group.by = "Clone")
#table(SAHST004$Clone, SAHST004$treatment_status)
length(which(SAHST004$Clone == "C1")) / length(SAHST004$Clone) # 0.7580372
length(which(SAHST004$Clone == "C2")) / length(SAHST004$Clone) # 0.08347434
length(which(SAHST004$Clone == "C3")) / length(SAHST004$Clone) # 0.1584884
saveRDS(SAHST004, "SAHST004_cloned.rds")
# SAHST004 <- readRDS("SAHST004_cloned.rds")

#Idents(SAHST004) <- "Clone"
SAHST004$Clone <- factor(SAHST004$Clone, levels = c("Non.tumour", "Peri.tumour", "C1"))
#SAHST004$Clone <- Idents(SAHST004)
p1 <- CellDimPlot(
  srt = SAHST004, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = rev(c("#9999cc", "#4d4d4d", "#d3d3d3")), 
  raster = T, 
  pt.size = 3, 
  #split.by = "treatment_status", 
  label_insitu = T) + theme(legend.position = "none")
p1

pdf("SAHST004_Clonal_DimPlot.pdf", height = 3, width = 3)
p1
dev.off()

### Heatmap
#SAHST004 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P5_cloned.rds"))
expr <- read.delim(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/SAHST004/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(SAHST004)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(SAHST004$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
    colnames(expr)[which(colnames(expr) %in% obs_cell)])]))

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df = kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

annotation_row = data.frame(
  #k_cluster =kmeans_df$k_cluster,
  group = kmeans_df$group
)
row.names(annotation_row) <- kmeans_df$cell_id
#color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
#names(color_cluster)=as.character(1:5)
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(group=c("C1" = "#9999cc",
                                            "Non.tumour" = "#d3d3d3", 
                                            "Peri.tumour" = "#4d4d4d"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_SAHST004.pdf",width = 6, height = 4)
ht = Heatmap(t(log2(expr))[rownames(annotation_row),],
             col = colorRamp2(c(-0.5,0,0.5), c("#375093","white","#831A21")),
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             show_row_names = F,
             column_split = factor(gene_pos$V2, new_cluster),
             heatmap_legend_param = list(title = "inferCNV",
                                         direction = "vertical",
                                         title_position = "leftcenter-rot",
                                         legend_height = unit(3, "cm")),
             left_annotation = left_anno, 
             row_title = NULL,
             column_title = NULL,
             top_annotation = top_color,
             use_raster = T, 
             border = T)
draw(ht, heatmap_legend_side = "right")
dev.off()

#CellDimPlot3D(srt = SAHST004, group.by = "Clone")
# remove object when done 

### Plot the clones
#color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
#Idents(ST_merge) <- "orig.ident"
#SAHST005 <- subset(ST_merge, idents = "SAHST004")
SAHST004@reductions$spatial = SAHST004@reductions$umap
SAHST004@reductions$spatial@key = 'spatial_'
SAHST004@reductions$spatial@cell.embeddings = as.matrix(SAHST004@images$image@coordinates[,c(3,2)])
SAHST004@reductions$spatial@cell.embeddings[,2] = -SAHST004@reductions$spatial@cell.embeddings[,2]
colnames(SAHST004@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
data_mat <- data.frame(SAHST004$Clone, #SAHST005$S_phase, SAHST004_plot$Os_min, 
                       SAHST004@reductions$spatial@cell.embeddings)

# Clone DimPlot 
# color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004.Clone)) +
  geom_point(shape = 16, size = 1) +
  scale_color_manual(values = #c("#9999cc", "#669933",# "#339999", 
                       #  "#339999")
                       rev(c("#9999cc", "#4d4d4d", "#d3d3d3"))) + 
  #scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
  #                                  "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
  #                                  "#d73027", "#a50026"), 
  #colours = c("lightblue","lightyellow","red"),  
  #                      limits = color_scale_val) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 127/77) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf('SAHST004_Clone_Spatial_mapping.pdf', height = 3, width = 3)
print(p1)
dev.off()

rm(SAHST004)

# SAHST003 ----- 
SAHST003 <- subset(ST_merge, idents = c("SAHST003"))
DimPlot(SAHST003)
DefaultAssay(SAHST003) <- "Spatial"
SAHST003@assays$SCT <- NULL
SAHST003 <- SCTransform(SAHST003, assay = "Spatial", verbose = FALSE)

SAHST003 <- RunPCA(SAHST003)
SAHST003 <- FindNeighbors(SAHST003, dims = 1:50, reduction = "pca")
SAHST003 <- RunUMAP(SAHST003, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(SAHST003, label = T, reduction = "umap.unintegrated")

## Cluster Ward 2D ----- ]
patient <- levels(ST_merge)[which(levels(ST_merge) == "SAHST003")] 
dir <- paste0(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/"), patient)
dd <- patient
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)

DimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table1)), reduction = "umap.unintegrated")
DimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table2)), reduction = "umap.unintegrated")
DimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table3)), reduction = "umap.unintegrated")
DimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table4)), reduction = "umap.unintegrated")
DimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table5)), reduction = "umap.unintegrated")

SpatialDimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table1)), image.alpha = 0.5)
SpatialDimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table2)), image.alpha = 0.5)
SpatialDimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table3)), image.alpha = 0.5)
SpatialDimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table4)), image.alpha = 0.5)
SpatialDimPlot(SAHST003, cells.highlight = WhichCells(SAHST003, rownames(cnv_table5)), image.alpha = 0.5)

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)

#"#8DD3C7" "#D1C2D2" "#D8C965" "#EB8E8B" "#FFED6F"
#1       0       0       0       0     495
#2       0       0      86       0       0
#3       0       0       0     155       0
#4       0     106       0       0       0
#5     353       0       0       0       0

cluster1 <- "C1"
cluster2 <- "C1"
cluster3 <- "Non.tumour"
cluster4 <- "Peri.tumour"
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
SAHST003 <- AddMetaData(SAHST003, cluster_group2)
DimPlot(SAHST003, group.by = "Clone", reduction = "umap.unintegrated")
SpatialDimPlot(SAHST003, group.by = "Clone")
#table(SAHST003$Clone, SAHST003$treatment_status)
length(which(SAHST003$Clone == "C1")) / length(SAHST003$Clone) # 0.7580372
length(which(SAHST003$Clone == "C2")) / length(SAHST003$Clone) # 0.08347434
length(which(SAHST003$Clone == "C3")) / length(SAHST003$Clone) # 0.1584884
saveRDS(SAHST003, "SAHST003_cloned.rds")
# SAHST003 <- readRDS("SAHST003_cloned.rds")

#Idents(SAHST003) <- "Clone"
SAHST003$Clone <- factor(SAHST003$Clone, levels = c("Non.tumour", "Peri.tumour", "C1", "C2"))
#SAHST003$Clone <- Idents(SAHST003)
p1 <- CellDimPlot(
  srt = SAHST003, group.by = "Clone", #stat.by = "Clone",
  reduction = "umap.unintegrated", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = rev(c("#669933", "#9999cc", "#4d4d4d", "#d3d3d3")), 
  raster = T, 
  pt.size = 3, 
  #split.by = "treatment_status", 
  label_insitu = T) + theme(legend.position = "none")
p1

pdf("SAHST003_Clonal_DimPlot.pdf", height = 3, width = 3)
p1
dev.off()

### Heatmap
#SAHST003 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P5_cloned.rds"))
expr <- read.delim(here("data/New_ST_251011/Output/Spatial_inferCNV_251015/SAHST003/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(SAHST003)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(SAHST003$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
    colnames(expr)[which(colnames(expr) %in% obs_cell)])]))

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df = kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

annotation_row = data.frame(
  #k_cluster =kmeans_df$k_cluster,
  group = kmeans_df$group
)
row.names(annotation_row) <- kmeans_df$cell_id
#color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
#names(color_cluster)=as.character(1:5)
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(group=c("C1" = "#9999cc",
                                            "C2" = "#669933",
                                            "Non.tumour" = "#d3d3d3", 
                                            "Peri.tumour" = "#4d4d4d"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_SAHST003.pdf",width = 6, height = 4)
ht = Heatmap(t(log2(expr))[rownames(annotation_row),],
             col = colorRamp2(c(-0.5,0,0.5), c("#375093","white","#831A21")),
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             show_row_names = F,
             column_split = factor(gene_pos$V2, new_cluster),
             heatmap_legend_param = list(title = "inferCNV",
                                         direction = "vertical",
                                         title_position = "leftcenter-rot",
                                         legend_height = unit(3, "cm")),
             left_annotation = left_anno, 
             row_title = NULL,
             column_title = NULL,
             top_annotation = top_color,
             use_raster = T, 
             border = T)
draw(ht, heatmap_legend_side = "right")
dev.off()

#CellDimPlot3D(srt = SAHST003, group.by = "Clone")
# remove object when done 

### Plot the clones
#color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
#Idents(ST_merge) <- "orig.ident"
#SAHST005 <- subset(ST_merge, idents = "SAHST003")
SAHST003@reductions$spatial = SAHST003@reductions$umap
SAHST003@reductions$spatial@key = 'spatial_'
SAHST003@reductions$spatial@cell.embeddings = as.matrix(SAHST003@images$image@coordinates[,c(3,2)])
SAHST003@reductions$spatial@cell.embeddings[,2] = -SAHST003@reductions$spatial@cell.embeddings[,2]
colnames(SAHST003@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
data_mat <- data.frame(SAHST003$Clone, #SAHST005$S_phase, SAHST003_plot$Os_min, 
                       SAHST003@reductions$spatial@cell.embeddings)

# Clone DimPlot 
# color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST003.Clone)) +
  geom_point(shape = 16, size = 1) +
  scale_color_manual(values = #c("#9999cc", "#669933",# "#339999", 
                       #  "#339999")
                       rev(c("#669933","#9999cc", "#4d4d4d", "#d3d3d3"))) + 
  #scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
  #                                  "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
  #                                  "#d73027", "#a50026"), 
  #colours = c("lightblue","lightyellow","red"),  
  #                      limits = color_scale_val) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 127/77) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf('SAHST003_Clone_Spatial_mapping.pdf', height = 3, width = 3)
print(p1)
dev.off()

rm(SAHST003)
