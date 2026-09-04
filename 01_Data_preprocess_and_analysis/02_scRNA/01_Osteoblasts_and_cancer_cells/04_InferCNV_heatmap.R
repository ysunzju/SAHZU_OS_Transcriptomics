library(here) # project-root-relative paths; run scripts from repository root
##################### InferCNV Heatmap (related to Figures S2D,E, S3A) ##################### 
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps"))
library(ComplexHeatmap)
library(circlize)

ref <- readRDS(here("data/InferCNV_rerun_for_phylo_22_06_25/ref.rds"))

# Demo of cancerous versus non-cancerous cells ----- 
ref_mat <- read.delim(here("data/InferCNV_31_01_25/H_AF_P5_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal/infercnv.references.txt"), sep = " ", check.names = F)
gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(ref_mat),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 
# Extract ref cell positions
ref_cell <- colnames(ref_mat)
#obs_cell <- colnames(H_P1)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(ref_mat)[which(colnames(ref_mat) %in% ref_cell)]),
  group = "Ref"#, 
  #  colnames(expr)[which(colnames(expr) %in% obs_cell)])])
)

# Clustering
set.seed(123)
kmeans.result <- kmeans(t(ref_mat), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df=kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

# Anno
annotation_row = data.frame(
  #k_cluster =kmeans_df$k_cluster,
  group = kmeans_df$group
)
row.names(annotation_row) <- kmeans_df$cell_id
#color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
#names(color_cluster)=as.character(1:5)
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(#group=c("C1" = "#9999cc",
                             #        "C2" = "#669933", 
                             #        "C3" = "#339999")#,
                             #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

# Plot 
pdf("CNV_heatmap_Ref_in_H_AF_P5.pdf",width = 6, height = 3)
ht = ComplexHeatmap::Heatmap(t(log2(ref_mat))[rownames(annotation_row),],
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
             border = T)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# Observations 
H_AF_P5 <- readRDS(here("data/InferCNV_31_01_25/All_OB_fin_withCNVmetadata.rds"))
Idents(H_AF_P5) <- "sample"
H_AF_P5 <- subset(H_AF_P5, idents = "H_AF_P5")
expr <- read.delim(here("data/InferCNV_31_01_25/H_AF_P5_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

# Extract ob cells
obs_cell <- colnames(H_AF_P5)

cell_anno = data.frame(cell_id = c(colnames(expr)[which(colnames(expr) %in% obs_cell)]#, 
                                   #colnames(ref_mat)[which(colnames(ref_mat) %in% ref_cell)]
),
group = c(H_AF_P5$malignancy_status[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)])])
)

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df=kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 
levels(kmeans_df$group) <- c("Malignant", "Non-malignant")

annotation_row = data.frame(
  #k_cluster =kmeans_df$k_cluster,
  group = kmeans_df$group
  
)

row.names(annotation_row) <- kmeans_df$cell_id
#color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
#names(color_cluster)=as.character(1:5)
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(group=c("Malignant" = "#1e2235",
                                            "Non-malignant" = "#8f9fab")#, 
                                    #        "C3" = "#339999")#,
                                    #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_obs_in_H_AF_P5.pdf",width = 6, height = 6)
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
             border = T)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# H_P1 ----- 
H_P1 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P1_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/HP1/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(H_P1)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(H_P1$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
    colnames(expr)[which(colnames(expr) %in% obs_cell)])]))

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df=kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

annotation_row = data.frame(
  #k_cluster =kmeans_df$k_cluster,
  group = kmeans_df$group
  
)
row.names(annotation_row) <- kmeans_df$cell_id
color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
names(color_cluster)=as.character(1:5)
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(group=c("C1" = "#9999cc",
                                            "C2" = "#669933", 
                                            "C3" = "#339999")#,
                                    #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_H_P1.pdf",width = 6, height = 4)
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
             border = T)
draw(ht, heatmap_legend_side = "right")
dev.off()
rm(H_P1)

# H_P3 ----- 
H_P3 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P3_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/HP3/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(H_P3)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(H_P3$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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

pdf("CNV_heatmap_H_P3.pdf",width = 6, height = 4)
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
             border = T)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# H_P5 ----- 
H_P5 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P5_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/HP5/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 
ref_cell <- colnames(ref)
obs_cell <- colnames(H_P5)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(H_P5$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
    colnames(expr)[which(colnames(expr) %in% obs_cell)])]))

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df = kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

#注释
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

pdf("CNV_heatmap_H_P5.pdf",width = 6, height = 4)
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
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# L_P2 ----- 
L_P2 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/L_P2_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/LP2/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(L_P2)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(L_P2$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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

pdf("CNV_heatmap_L_P2.pdf",width = 6, height = 4)
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
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# L_P5----- 
L_P5 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/L_P5_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/LP5/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 
ref_cell <- colnames(ref)
obs_cell <- colnames(L_P5)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(L_P5$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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

pdf("CNV_heatmap_L_P5.pdf",width = 6, height = 4)
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
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# L_P6----- 
L_P6 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/L_P6_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/LP6/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(L_P6)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(L_P6$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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

pdf("CNV_heatmap_L_P6.pdf",width = 6, height = 4)
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

# Baseline only samples 
# H_P2 ----- 
H_P2 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P2_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/HP2/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

ref_cell <- colnames(ref)
obs_cell <- colnames(L_P6)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(H_P2$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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

pdf("CNV_heatmap_H_P2.pdf",width = 6, height = 4)
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

# H_P4 ----- 
H_P4 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps/H_P4_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/HP4/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 
obs_cell <- colnames(H_P4)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(H_P4$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(group=c("C1" = "#9999cc",
                                            "C2" = "#669933", 
                                            "C3" = "#339999"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_H_P4.pdf",width = 6, height = 4)
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

# H_P6 ----- 
H_P6 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps/H_P6_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/HP6/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

obs_cell <- colnames(H_P6)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(H_P6$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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

pdf("CNV_heatmap_H_P6.pdf",width = 6, height = 4)
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
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# H_P7 ----- 
H_P7 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps/H_P7_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/HP7/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

obs_cell <- colnames(H_P7)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(H_P7$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
    colnames(expr)[which(colnames(expr) %in% obs_cell)])]))

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df = kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

#注释
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

#绘制热图
pdf("CNV_heatmap_H_P7.pdf",width = 6, height = 4)
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

# H_P8 ----- 
H_P8 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Ward2D_phylogeny/H_P8_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/HP8/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

obs_cell <- colnames(H_P8)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(H_P8$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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
                                            "C3" = "#339999",
                                            "C4" = "#CCCCCC"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_H_P8.pdf",width = 6, height = 4)
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
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# L_P1 ----- 
L_P1 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps/L_P1_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/LP1/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

obs_cell <- colnames(L_P1)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(L_P1$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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
                                            "C3" = "#339999",
                                            "C4" = "#CCCCCC"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_L_P1.pdf",width = 6, height = 4)
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

# L_P3 ----- 
#L_P3 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps/L_P3_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/LP3/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

OB_cancer <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Malig_OB_for_NMF_sample_w_too_few_cells_removed.rds"))
L_P3 <- subset(OB_cancer, idents = c("L_BF_P3"))
obs_cell <- colnames(L_P3)
L_P3$Clone <- "C1"

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(L_P3$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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
                                            "C3" = "#339999",
                                            "C4" = "#CCCCCC"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_L_P3.pdf",width = 6, height = 4)
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
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# L_P4 ----- 
#L_P3 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps/L_P3_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/LP4/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

OB_cancer <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Malig_OB_for_NMF_sample_w_too_few_cells_removed.rds"))
L_P4 <- subset(OB_cancer, idents = c("L_BF_P4"))
obs_cell <- colnames(L_P4)
L_P4$Clone <- "C1"

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(L_P4$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
    colnames(expr)[which(colnames(expr) %in% obs_cell)])]))

set.seed(123)
kmeans.result <- kmeans(t(expr), 5)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df = kmeans_df%>%inner_join(cell_anno,by="cell_id") %>% arrange(group)
kmeans_df$group=as.factor(kmeans_df$group) 

#注释
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
                                            "C3" = "#339999",
                                            "C4" = "#CCCCCC"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_L_P4.pdf",width = 6, height = 4)
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
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
dev.off()

# L_P7 ----- 
L_P7 <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps/L_P7_cloned.rds"))
expr <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/LP7/infercnv.observations.txt"), sep = " ", check.names = F)

gene_pos <- read.delim(here("data/InferCNV_rerun_for_phylo_22_06_25/hg38_gencode_v27.txt"), header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]
new_cluster <- unique(gene_pos$V2)
top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(color = c("white")),
                                                    labels_gp = gpar(cex = 1, col = c("black" ), fontface = "bold"),
                                                    height = unit(5,"mm"))) 

obs_cell <- colnames(L_P7)

cell_anno = data.frame(cell_id = c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
  colnames(expr)[which(colnames(expr) %in% obs_cell)]),
  group =c(L_P7$Clone[c(#colnames(expr)[which(colnames(expr) %in% ref_cell)], 
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
                                            "C3" = "#339999",
                                            "C4" = "#CCCCCC"
                           )#,
                           #k_cluster=color_cluster
                           ),
                           show_annotation_name = F)

pdf("CNV_heatmap_L_P7.pdf",width = 6, height = 4)
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
