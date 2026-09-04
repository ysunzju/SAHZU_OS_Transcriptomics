library(here) # project-root-relative paths; run scripts from repository root
##################### Figure 1 ##################### 
# All cells plot updated 251107
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/All_cells_F1_main"))

# Load combined file 
library(Seurat)
library(tidyverse)
#seuratObj <- readRDS(here("data/Correlation_and_crosstalk_analysis_18_07_25/All_cells_for_ecoTyper_cell_type_adjusted_250828.rds"))
seuratObj <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))

##################### Dimplot of all cells (Figure 1B) ##################### 
#DimPlot(seuratObjobj, group.by = "cluster")
DimPlot(seuratObj, group.by = "cluster")
dim(seuratObj)

colors <- c('#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',
               '#CC6677', '#882255', '#AA4499', '#662288', '#BB5566',
               '#5D3A9B', '#997700', '#DDDDDD', '#004488',
               '#6699CC', '#66BBAA') 

p1 <- DimPlot(seuratObj, reduction = "umap", group.by = c("cluster"), label = T) + 
  scale_color_manual(values = colors) + coord_fixed(ratio = 0.9) #+ 
  #NoLegend()
p1

pdf(paste0("All_cells", "_UMAP_v251108.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

##################### Annotated per-sample cell count plot (F1c) ##################### 
# Info plot ----- 
## Paired cell proportion plot -----
## Cell prop (malig vs non-malig) plot -----
names(colors) <- levels(seuratObj$cluster)
seuratObj$cluster_level_adjusted <- factor(seuratObj$cluster, 
                                           levels = rev(c("Osteoblasts", "MSCs", "Fibroblasts", "Chondrocytes", 
                                                      "MuralCells", "ECs", "Osteoclasts", "MastCells", "Myeloid", "TandNK", 
                                                      "Pro_B", "BCells", "PlasmaCells", "pDCs", "HSCs", "Erythrocytes"
                                             ))) 
#mat <- data.frame(table(seuratObj$sample, seuratObj$cluster))
mat <- data.frame(table(seuratObj$sample, seuratObj$cluster_level_adjusted))
mat$patient <- gsub("_AF_", "", mat$Var1)
mat$patient <- gsub("_BF_", "", mat$patient)
colnames(mat)[1:3] <- c("celltype","Patient","cell_num")

library(wesanderson)
library(ggpubr)
library(ggalluvial)

levels(mat$celltype)
mat$celltype <- factor(mat$celltype, levels = c("H_BF_P1", "H_AF_P1", "H_BF_P2", "H_AF_P2", "H_BF_P3", "H_AF_P3", 
                                                   "H_BF_P4", "H_AF_P4", "H_BF_P5", "H_AF_P5", "H_BF_P6", "H_AF_P6", 
                                                   "H_BF_P7", "H_AF_P7", "H_BF_P8", "H_AF_P8", "L_BF_P1", "L_AF_P1", 
                                                   "L_BF_P2", "L_AF_P2", "L_BF_P3", "L_AF_P3", "L_BF_P4", "L_AF_P4", 
                                                   "L_BF_P5", "L_AF_P5", "L_BF_P6", "L_AF_P6", "L_BF_P7", "L_AF_P7"))

p2 <- ggbarplot(mat, x = "celltype", y="cell_num", color="black", fill="Patient",
          font.main = c(8,"bold", "black"), font.x = c(8, "bold"), 
          font.y=c(8,"bold")) + 
  theme_bw() + 
  rotate_x_text() + 
  scale_fill_manual(values=colors) + 
  labs(x = "Sample", y = "Cell count") + 
  theme(axis.text.x=element_blank(),
        panel.grid = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold")
        ) 

# Info Plot (barplot w/ anno) -----
meta_data_all_cells <- seuratObj@meta.data
meta_data_all_cells = meta_data_all_cells[!duplicated(meta_data_all_cells$sample),]
meta_data_all_cells$treatment_status <- ifelse(meta_data_all_cells$group %in% c("A", "C"), "Pre", "Post")  

## Rx status ----- 
mat$group_info <- NA
for (i in unique(mat$celltype)) {
  mat$group_info[which(mat$celltype == i)] <- as.character(meta_data_all_cells$treatment_status[which(meta_data_all_cells$sample == i)])
}
mat$group_info <- factor(mat$group_info, levels = c("Post", "Pre"), ordered = T)
Rx_status_Group_info <- ggplot(mat, aes(x = celltype, y = 1)) +  
  geom_point(aes(fill = group_info), size = 5, shape = 22, 
             stroke = NA) +  
  scale_x_discrete(breaks = c(1)) +  
  scale_fill_manual(name = "", values = rev(c("#8f476D", "#F0CFE3"))) +
  theme_void() +  
  theme(legend.position = "top",        
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
Rx_status_Group_info

## MF (Sex) ----- 
Idents(seuratObj) <- "reportname"
seuratObj <- RenameIdents(seuratObj,
                      "SAH_001" = "M", 
                      "SAH_002" = "M", 
                      "SAH_003" = "F",
                      "SAH_004" = "F",
                      "SAH_005" = "M",
                      "SAH_006" = "M", 
                      "SAH_007" = "M",
                      "SAH_008" = "F",
                      "SAH_009" = "M",
                      "SAH_010" = "M",
                      "SAH_011" = "F", 
                      "SAH_012" = "F",
                      "SAH_013" = "M",
                      "SAH_014" = "M",
                      "SAH_015" = "M", 
                      "SAH_016" = "F", 
                      "SAH_017" = "M", 
                      "SAH_018" = "M", 
                      "SAH_019" = "F", 
                      "SAH_020" = "F", 
                      "SAH_021" = "M", 
                      "SAH_022" = "F", 
                      "SAH_023" = "M", 
                      "SAH_024" = "M", 
                      "SAH_025" = "M", 
                      "SAH_026" = "F", 
                      "SAH_027" = "F", 
                      "SAH_028" = "M", 
                      "SAH_029" = "F", 
                      "SAH_030" = "M")
seuratObj$Gender <- Idents(seuratObj)

meta_data_all_cells <- seuratObj@meta.data
meta_data_all_cells = meta_data_all_cells[!duplicated(meta_data_all_cells$sample),]
mat$group_info <- NA
for (i in unique(mat$celltype)) {
  mat$group_info[which(mat$celltype == i)] <- as.character(meta_data_all_cells$Gender[which(meta_data_all_cells$sample == i)])
}
mat$group_info <- factor(mat$group_info, levels = c("M", "F"), ordered = T)
Sex_Group_info <- ggplot(mat, aes(x = celltype, y = 1)) +  
  geom_point(aes(fill = group_info), size = 5, shape = 22, 
             stroke = NA) +  
  scale_x_discrete(breaks = c(1)) +  
  scale_fill_manual(name = "", values = rev(c("#ECE093", "#8F6D5F"))) + 
  theme_void() +  
  theme(legend.position = "top",        
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
Sex_Group_info

## Age -----  
Idents(seuratObj) <- "reportname"
seuratObj <- RenameIdents(seuratObj,
                      "SAH_001" = 20, 
                      "SAH_002" = 17,
                      "SAH_003" = 14,
                      "SAH_004" = 14,
                      "SAH_005" = 28,
                      "SAH_006" = 11, 
                      "SAH_007" = 18,
                      "SAH_008" = 12,
                      "SAH_009" = 20, 
                      "SAH_010" = 17,  
                      "SAH_011" = 14, 
                      "SAH_012" = 14,
                      "SAH_013" = 28,
                      "SAH_014" = 11,
                      "SAH_015" = 18, 
                      "SAH_016" = 12, 
                      "SAH_017" = 13, 
                      "SAH_018" = 16, 
                      "SAH_019" = 9, 
                      "SAH_020" = 12, 
                      "SAH_021" = 7, 
                      "SAH_022" = 16, 
                      "SAH_023" = 11, 
                      "SAH_024" = 13, 
                      "SAH_025" = 16, 
                      "SAH_026" = 9, 
                      "SAH_027" = 12, 
                      "SAH_028" = 7, 
                      "SAH_029" = 16, 
                      "SAH_030" = 11)
seuratObj$Age <- Idents(seuratObj)

meta_data_all_cells <- seuratObj@meta.data
meta_data_all_cells = meta_data_all_cells[!duplicated(meta_data_all_cells$sample),]
mat$group_info <- NA
for (i in unique(mat$celltype)) {
  mat$group_info[which(mat$celltype == i)] <- as.character(meta_data_all_cells$Age[which(meta_data_all_cells$sample == i)])
}
mat$group_info <- as.numeric(mat$group_info)

Age_Group_info <- ggplot(mat, aes(x = celltype, y = 1)) +  
  geom_point(aes(fill = group_info), size = 5, shape = 22, 
             stroke = NA) +  
  scale_x_discrete(breaks = c(1)) +  
  scale_fill_gradientn(colors = c("#eeebe9", "#d7ccc7", "#bbaaa4", "#a0877f", "#8c6d63", 
                                  "#795447", "#6c4c40", "#5d3f37", "#4d332d", "#3e2622")) +
  theme_void() +  
  theme(legend.position = "top",        
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
Age_Group_info

## TNR ----- 
Idents(seuratObj) <- "reportname"
seuratObj <- RenameIdents(seuratObj,
                      "SAH_001" = 99.6,
                      "SAH_002" = 99.8, 
                      "SAH_003" = 100,
                      "SAH_004" = 100,
                      "SAH_005" = 99.7,
                      "SAH_006" = 99.5, 
                      "SAH_007" = 98.9,
                      "SAH_008" = 94.1,
                      "SAH_009" = 99.6,
                      "SAH_010" = 99.8,
                      "SAH_011" = 100, 
                      "SAH_012" = 100,
                      "SAH_013" = 99.7,
                      "SAH_014" = 99.5,
                      "SAH_015" = 98.9, 
                      "SAH_016" = 94.1,
                      "SAH_017" = 70.1, 
                      "SAH_018" = 47.8, 
                      "SAH_019" = 81.4, 
                      "SAH_020" = 85.1,
                      "SAH_021" = 76.2, 
                      "SAH_022" = 68.6, 
                      "SAH_023" = 80.5, 
                      "SAH_024" = 70.1, 
                      "SAH_025" = 47.8, 
                      "SAH_026" = 81.4, 
                      "SAH_027" = 85.1,
                      "SAH_028" = 76.2, 
                      "SAH_029" = 68.6, 
                      "SAH_030" = 80.5)
seuratObj$TNR_val <- Idents(seuratObj)

meta_data_all_cells <- seuratObj@meta.data
meta_data_all_cells = meta_data_all_cells[!duplicated(meta_data_all_cells$sample),]
mat$group_info <- NA
for (i in unique(mat$celltype)) {
  mat$group_info[which(mat$celltype == i)] <- as.character(meta_data_all_cells$TNR_val[which(meta_data_all_cells$sample == i)])
}
mat$group_info <- as.numeric(mat$group_info)

TNR_Group_info <- ggplot(mat, aes(x = celltype, y = 1)) +  
  geom_point(aes(fill = group_info), size = 5, shape = 22, 
             stroke = NA) +  
  scale_x_discrete(breaks = c(1)) +  
  scale_fill_gradientn(colors = c("#eeecdf", "#becdd2", "#6f9ad1", "#44679f", "#3f4f71")) +
  theme_void() +  
  theme(legend.position = "top",        
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
TNR_Group_info

## Anatomical position ----- 
Idents(seuratObj) <- "reportname"
seuratObj <- RenameIdents(seuratObj,
                          "SAH_001" = "Knee", 
                          "SAH_002" = "Tibia",
                          "SAH_003" = "Femur",
                          "SAH_004" = "Femur",
                          "SAH_005" = "Femur",
                          "SAH_006" = "Femur", 
                          "SAH_007" = "Femur",
                          "SAH_008" = "Tibia",
                          "SAH_009" = "Knee", 
                          "SAH_010" = "Tibia", 
                          "SAH_011" = "Femur", 
                          "SAH_012" = "Femur",
                          "SAH_013" = "Femur",
                          "SAH_014" = "Femur",
                          "SAH_015" = "Femur", 
                          "SAH_016" = "Tibia", 
                          "SAH_017" = "Femur", 
                          "SAH_018" = "Fibula", 
                          "SAH_019" = "Femur", 
                          "SAH_020" = "Tibia",
                          "SAH_021" = "Femur", 
                          "SAH_022" = "Femur", 
                          "SAH_023" = "Femur", 
                          "SAH_024" = "Femur", 
                          "SAH_025" = "Fibula", 
                          "SAH_026" = "Femur", 
                          "SAH_027" = "Tibia", 
                          "SAH_028" = "Femur", 
                          "SAH_029" = "Femur", 
                          "SAH_030" = "Femur")
seuratObj$Anatomy <- Idents(seuratObj)

meta_data_all_cells <- seuratObj@meta.data
meta_data_all_cells = meta_data_all_cells[!duplicated(meta_data_all_cells$sample),]
mat$group_info <- NA
for (i in unique(mat$celltype)) {
  mat$group_info[which(mat$celltype == i)] <- as.character(meta_data_all_cells$Anatomy[which(meta_data_all_cells$sample == i)])
}

Anatomy_Group_info <- ggplot(mat, aes(x = celltype, y = 1)) +  
  geom_point(aes(fill = group_info), size = 5, shape = 22, 
             stroke = NA) +  
  scale_x_discrete(breaks = c(1)) + 
  scale_fill_manual(name = "", values = rev(c("#9F6D97", "#FBECF1", "#F9EAE0", "#E4EFF9"))) + 
  theme_void() +  
  theme(legend.position = "top",        
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
Anatomy_Group_info

## +/- Met at diagnosis ----- 
Idents(seuratObj) <- "reportname"
seuratObj <- RenameIdents(seuratObj,
                          "SAH_001" = "No", 
                          "SAH_002" = "No", 
                          "SAH_003" = "No",
                          "SAH_004" = "No",
                          "SAH_005" = "No",
                          "SAH_006" = "No", 
                          "SAH_007" = "No",
                          "SAH_008" = "Yes",
                          "SAH_009" = "No", 
                          "SAH_010" = "No", 
                          "SAH_011" = "No", 
                          "SAH_012" = "No",
                          "SAH_013" = "No",
                          "SAH_014" = "No",
                          "SAH_015" = "No", 
                          "SAH_016" = "Yes", 
                          "SAH_017" = "No", 
                          "SAH_018" = "Yes", 
                          "SAH_019" = "No", 
                          "SAH_020" = "No", 
                          "SAH_021" = "No", 
                          "SAH_022" = "No", 
                          "SAH_023" = "No", 
                          "SAH_024" = "No", 
                          "SAH_025" = "Yes", 
                          "SAH_026" = "No", 
                          "SAH_027" = "No", 
                          "SAH_028" = "No", 
                          "SAH_029" = "No", 
                          "SAH_030" = "No")
seuratObj$Met_status <- Idents(seuratObj)

meta_data_all_cells <- seuratObj@meta.data
meta_data_all_cells = meta_data_all_cells[!duplicated(meta_data_all_cells$sample),]
mat$group_info <- NA
for (i in unique(mat$celltype)) {
  mat$group_info[which(mat$celltype == i)] <- as.character(meta_data_all_cells$Met_status[which(meta_data_all_cells$sample == i)])
}

Met_Group_info <- ggplot(mat, aes(x = celltype, y = 1)) +  
  geom_point(aes(fill = group_info), size = 5, shape = 22, 
             stroke = NA) +  
  scale_x_discrete(breaks = c(1)) + 
  scale_fill_manual(name = "", values = rev(c("black", "grey"))) + 
  theme_void() +  
  theme(legend.position = "top",        
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
Met_Group_info

## Plot combined ----- 
library(patchwork)
p3 <- Age_Group_info + Sex_Group_info + Rx_status_Group_info + Anatomy_Group_info + Met_Group_info + #Rx_eff_Group_info + 
  TNR_Group_info + #p1 + 
  p2 +
  plot_layout(ncol = 1, heights = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.15), guides = "collect"
              ) + 
  theme(legend.position = "none", legend.margin = margin(t = -0.5, unit = "cm"),
        legend.text = element_text(#angle = 0,
          face = "italic",
          size = 5#,
          #hjust= 0,
          #vjust = 0
        ))
p3

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Cell_prop_w_anno_251110.pdf"), 
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches
print(p3, newpage = FALSE)
dev.off()

##################### F1d ##################### 
# Expression heatmap ----- 
library(ComplexHeatmap)
seuratObj <- NormalizeData(seuratObj) %>% FindVariableFeatures()%>% ScaleData()
seuratObj$group_anno <- ifelse(seuratObj$group == "A", "H_Pre", 
                               ifelse(seuratObj$group == "B", "H_Post", 
                                      ifelse(seuratObj$group == "C", "L_Pre", "L_Post")))
seuratObj$group_anno <- factor(seuratObj$group_anno, levels = c("H_Pre", "H_Post", "L_Pre", "L_Post"))
cols_for_group <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
names(cols_for_group) = c("H_Pre", "H_Post", "L_Pre", "L_Post")
seuratObj$cluster_level_adjusted <- factor(seuratObj$cluster, 
                                           levels = c("Osteoblasts", "MSCs", "Fibroblasts", "Chondrocytes", 
                                                          "MuralCells", "ECs", "Osteoclasts", "MastCells", "Myeloid", "TandNK", 
                                                          "Pro_B", "BCells", "PlasmaCells", "pDCs", "HSCs", "Erythrocytes"
                                           )) 
seuratObj$cluster_level_adjusted <- droplevels(seuratObj$cluster_level_adjusted)


# Load function. Credits to KS生信科研课@Wechat
ks_single_marker_heatmap <- function(object, 
                                     assay,
                                     genes,
                                     bar, 
                                     cluster_order=NULL,
                                     border = T,
                                     column_gap=0,
                                     show_rownames = T, 
                                     cols_anno=NULL, 
                                     heat_cols=colorRampPalette(c("dodgerblue3","white","Tomato2"))(n=100),
                                     anno_enrich = F,
                                     cluster_gene_length=NULL,
                                     enrich_terms=NULL,
                                     rownames_anno=NULL,
                                     min = -2.5, max = 2.5, fontsize = 12){
  
  # unique genes
  genes <- intersect(genes, rownames(object))
  
  # exp data
  if(is_empty(GetAssayData(object,assay, layer = 'data'))){
    
    exp <- GetAssayData(object, assay, layer = "counts")[genes,]
    exp <- log10(exp + 1)
    exp <- scale(exp, center = TRUE, scale = TRUE)
    
  }else{
    
    exp <- GetAssayData(object, slot = "data")[genes,]
    
  }
  
  pmed <- function(x) pmin(pmax(x, min), max)
  exp <- apply(exp, 2, pmed)
  
  # Order cells
  if(length(bar) == 1) {
    
    bar_col <- object@meta.data[,c(bar)] %>% as.data.frame()
    rownames(bar_col) <- rownames(object@meta.data)
    colnames(bar_col) <- bar
    
  }else {
    
    bar_col <- object@meta.data[, bar]
    
  }
  
  if(!is.null(cluster_order)){

    bar_col[,1] <- factor(bar_col[,1], levels = cluster_order)
    bar_col <- bar_col[order(bar_col[,1]), ,drop=F]
    
  }
  
  exp <- exp[ , rownames(bar_col)]
  
  if(border == T){
    
    border1  = TRUE
    border2  ='black'
    
  }
  
  # Heatmap annotation
  anno_top <- HeatmapAnnotation(df = bar_col, 
                                col = cols_anno,
                                border = border1,
                                show_annotation_name = F,
                                annotation_name_gp = gpar(fontsize = (fontsize-2)),
                                simple_anno_size = unit(3, "mm"),
                                annotation_legend_param = list(title_gp = gpar(fontsize = (fontsize-2), fontface = "bold"),
                                                               labels_gp = gpar(fontsize = (fontsize-2)),
                                                               grid_height = unit(3.5, "mm"),
                                                               grid_width = unit(3.5, "mm")))

  # legend
  heatmap_legend_param = list(color_bar = "continuous",
                              legend_direction = "vertical",
                              legend_width = unit(0.5, "cm"),
                              legend_height = unit(2, "cm"),
                              title = "Z-Score",
                              title_position="leftcenter-rot",
                              border ="black",
                              labels_gp = gpar(fontsize = 8,col='black',font = 3))
  
  
  
  if(!is.null(rownames_anno)){
    
    gene_indices <- match(rownames_anno, rownames(exp))
    
    ha <- rowAnnotation(
      foo = anno_mark(
        at = gene_indices, 
        labels = rownames_anno, 
        side = "left",
        labels_gp = gpar(fontsize = 10, col = "black") 
      )
    )
    
  }
  
  if(anno_enrich == F){
    
    Heatmap(exp,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_column_names = FALSE,
            show_row_names	= show_rownames,
            column_split = bar_col[, 1],
            top_annotation = anno_top,
            row_names_side = "left",
            column_title = NULL,
            left_annotation  = ha,
            heatmap_legend_param = heatmap_legend_param,
            col = heat_cols,
            border = border2,
            row_names_gp = gpar(fontsize = 8),
            column_gap = unit(column_gap, "pt"))
    
  }else{
    
    
    align_to = split(seq_len(nrow(exp)), rep(cluster_order,cluster_gene_length))
    
    
    Heatmap(exp,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_column_names = FALSE,
            show_row_names	= show_rownames,
            column_split = bar_col[, 1],
            row_names_side = "left",
            top_annotation = anno_top,
            column_title = NULL,
            left_annotation = ha,
            heatmap_legend_param = heatmap_legend_param,
            col = heat_cols,
            border = border2,
            row_names_gp = gpar(fontsize = 8),
            column_gap = unit(column_gap, "pt"),
            right_annotation  = rowAnnotation(textbox = anno_textbox(align_to = align_to,
                                                                     enrich_terms,
                                                                     by="anno_block",
                                                                     just  = 'right',
                                                                     background_gp = gpar(fill = c("white","white"),
                                                                                          col = c("white")),
                                                                     gp=gpar(fontsize = 10,col = "black"),
                                                                     max_width = unit(5, "mm"))))
    
    
  }
  
  
}

p1 <- ks_single_marker_heatmap(object=seuratObj, assay = 'RNA',
                         genes = c("ALPL", "RUNX2", "SFRP2", "THY1", #"LUM", 
                                 "TAGLN", 
                                 "SOX9", "ACAN", "RGS5", "ACTA2", "PECAM1", "PLVAP", 
                                 "ATP6V0D2", "ACP5", "PTPRC", "TPSAB1", "CPA3", "CD14", "LYZ", "CD2","CD3D", 
                                 "NKG7", "VPREB1", "CD79A", "MS4A1", "MZB1", "JCHAIN", 
                                 "IRF7", "LILRA4", "CDK6","SPINK2", "ALAS2", "HBB"
                                 ), bar = c("cluster_level_adjusted", "group_anno"),
                       show_rownames = T,
                       cols_anno = list("cluster_level_adjusted" = colors, "group_anno" = cols_for_group), 
                       heat_cols = colorRampPalette(c("#e9f4f6", "#8f9fab", "#1e2235"))(100),
                       rownames_anno = c("")
                       )

p1
pdf("All_cells_heatmap_251111.pdf", width = 10, height = 3)
print(p1)
dev.off()



##################### F1f and Extended Data Fig 1e ##################### 
# Changes in cellular proportions ----- 
Idents(seuratObj) <- "group"
seuratObj <- RenameIdents(seuratObj, "A" = "H_Pre", 
                      "B" = "H_Post", 
                      "C" = "L_Pre", 
                      "D" = "L_Post")
seuratObj$group_anno <- Idents(seuratObj)
seuratObj$cluster_level_adjusted <- factor(seuratObj$cluster, 
                                           levels = rev(c("Osteoblasts", "MSCs", "Fibroblasts", "Chondrocytes", 
                                                          "MuralCells", "ECs", "Osteoclasts", "MastCells", "Myeloid", "TandNK", 
                                                          "Pro_B", "BCells", "PlasmaCells", "pDCs", "HSCs", "Erythrocytes"
                                           ))) 
seuratObj$cluster_level_adjusted <- droplevels(seuratObj$cluster_level_adjusted)
Cellratio <- prop.table(table(seuratObj$cluster_level_adjusted, seuratObj$sample), margin = 2) 
Cellratio <- data.frame(Cellratio)

library(reshape2)
cellper <- dcast(Cellratio, Var2~Var1, value.var="Freq")
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]

meta <- seuratObj@meta.data
colnames(meta)
meta <- meta[,c(24,4)] # Selecting columns with group and sample information 
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
meta <- unique(meta)
cellper <- merge(cellper, meta, by = "sample")
colnames(cellper)
colnames(cellper)[18] <- "group" # name change - group anno to group 
cellper <- as.data.frame(cellper)

pplist =list()
seuratObj_groups = unique(levels(seuratObj$cluster_level_adjusted))

library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
library(cols4all)
library(tidyverse)
library(ggplot2)
library(gghalves)
library(ggridges)
library(cols4all)
library(ggplot2)
library(patchwork)
library(grid)
library(ggtext)
library(ggsignif)
library(stringr)

mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))

for(group_ in seuratObj_groups){
  cellper_ = cellper %>% dplyr::select(one_of(c('sample', 'group', group_)))
  colnames(cellper_)=c('sample','group','percent')
  cellper_$percent =as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper = quantile(percent,0.75),
                                                      lower = quantile(percent,0.25),
                                                      mean = mean(percent),
                                                      median = median(percent),
                                                      lower_lim = range(percent)[1], 
                                                      upper_lim = range(percent)[2])
  y_limits = c(0, 
               min(1, max(cellper_$percent)*1.2)) 
  y_breaks = round(seq(0, 
                       min(1, max(cellper_$percent)*1.2), 
                       length.out = 5), 2) 
  p_position = min(1, max(cellper_$percent)*1.2) * 0.9
  p_position1 = min(1, max(cellper_$percent)*1.2) * 0.85
  
  print(group_)
  print(cellper_$median)
  mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + 
    geom_jitter(shape =21, aes(fill = group), width = 0.25) + 
    geom_boxplot(alpha = 0.5, size = 1, outlier.shape = NA, aes(fill = group)) + 
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) + 
    theme_cowplot() +
    theme(axis.text = element_text(size =10),axis.title = element_text(size = 10),legend.text = element_text(size =10),
          legend.title = element_text(size =10),plot.title = element_text(size =10,face ='plain'),legend.position ='none')+
    labs(title = group_, y = "Percentage", x = "Group") + 
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) + 
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
          axis.text.y = element_text(color = "black", size = 9),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) + 
    theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid"))
  
  labely = max(cellper_$percent)
  compare_means(percent ~group,  data = cellper_)
  my_comparisons <-list(c("H_Pre","H_Post") 

  )
  pp1 = pp1 + geom_signif(comparisons = list(c("H_Pre","H_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position, 
                          test.args = list(paired = TRUE, 
                                           alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("L_Pre","L_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))}, 
                          y_position = p_position, 
                          test.args = list(paired = TRUE, 
                                           alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("H_Pre","L_Pre")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position1, 
                          test.args   = list(#paired = TRUE, 
                                             #alternative = "two.sided"
                                             ),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black") 
  
  pp1 <- pp1 + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  pplist[[group_]]= pp1
}
names(pplist)
pps <- plot_grid(pplist[[1]], 
                 pplist[[2]],
                 pplist[[3]], 
                 pplist[[4]],
                 pplist[[5]],
                 pplist[[6]],
                 pplist[[7]],
                 pplist[[8]],
                 pplist[[9]],
                 pplist[[10]],
                 pplist[[11]], 
                 pplist[[12]],
                 pplist[[13]],
                 pplist[[14]],
                 pplist[[15]],
                 pplist[[16]], 
                 ncol = 6 
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_all_cell_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 9, 
    height = 6.75) 
print(pps, newpage = FALSE)
dev.off()

##################### Bulk RNA-seq heatmap (F1g) #####################
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/All_cells_F1_main"))
library(tidyverse)

expr_combat <- read.table(here("data/EcoTyper_run_16_06_25/OS_182_post_batch_rem_260731.txt"))
group <- read.csv(here("data/Re-clustering_fine_cell_types_21_05_25/Metadata_260802/updated_bulk_metadta_260802.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")

geneset <- gson::read.gmt(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/h.all.v2025.1.Hs.symbols.gmt.txt"))
table(geneset$term)
geneset$term <- gsub(pattern = "HALLMARK_","", geneset$term)
geneset <- geneset[, c(2, 1)]
colnames(geneset) <- c("Metagene", "Cell type") 

geneSet = split(1:dim(geneset)[1], geneset$`Cell type`) %>%
  lapply(function(x) {
    geneset[unlist(x), 1]
  })
head(geneSet, 4) 

group$TNR <- as.numeric(group$TNR)
group$Sample[!is.na(group$TNR)]
dat <- as.matrix(expr_combat[, -1]) 
group_sel <- group[!is.na(group$TNR), ]
dat <- dat[, colnames(dat) %in% group$Sample[!is.na(group$TNR)]]
geneSet1 <- geneSet[names(geneSet)[c(1, 4, 7, 9, 11, 12, 13, 14, 17, 18, 19, 22, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 
                                     37,  42, 44, 45, 46, 49)]]

library(GSVA) 
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet1,
                       normalize = TRUE) 
ssgsea <- gsva(gsvaPar, verbose = FALSE) 
dt <- scale(t(ssgsea))
dt <- t(dt)

library(ComplexHeatmap)
library(circlize)
dt <- dt[, group_sel$Sample]
identical(group_sel$Sample, colnames(dt))

group1 <- group_sel[, c("Sample","Treatment", "Progression", "TNR")]
group1$Progression <- ifelse(group1$Progression == "1", "Progression", "Non-progression")
group1$TNR <- as.numeric(group1$TNR)
group1 <- group1[order(group1$TNR, decreasing = T), ]
group1 <- group1[order(group1$Progression), ]
group1 <- group1[order(group1$Treatment, decreasing = T), ]

dt <- dt[, group1$Sample]
identical(group_sel$Sample, colnames(dt)) 
identical(group1$Sample, colnames(dt))
group2 <- group1[, c(2, 3, 4)]
identical(group1$Sample, colnames(dt)) 
identical(group2$Treatment, group1$Treatment) 
group2$TNR <- as.numeric(group2$TNR)
group2$TNR_cat <- ifelse(group2$TNR >= 90, "R", "NR")

group3 <- group2
group3$split <- NA
group3$split[group3$TNR_cat == "R" & group3$Treatment == "Pre"] <- "A"
group3$split[group3$TNR_cat == "NR" & group3$Treatment == "Pre"] <- "B"
group3$split[group3$TNR_cat == "R" & group3$Treatment == "Post"] <- "C"
group3$split[group3$TNR_cat == "NR" & group3$Treatment == "Post"] <- "D"
group3$split <- factor(group3$split, levels = c("A", "C", "B", "D"))

fontsize = 8
col_anno <- HeatmapAnnotation(df = group2, 
                              show_annotation_name = T, 
                              gp = gpar(fontsize = 10), 
                              na_col = "black",  
                              col = list(Treatment = c( 
                                'Pre' = '#90486e',  
                                'Post' = '#efcfe3'),  
                                Progression = c( 
                                  "Non-progression" = "#bfbebe",  
                                  "Progression" = "#010101"),  
                                TNR = colorRamp2(c(0, 20, 40, 60, 80, 90, 100),  
                                                 c("#eeecdf","#eeecdf", "#becdd2", "#6f9ad1", "#44679f", "#3f4f71", "#3f4f71")), 
                                TNR_cat = c( 
                                  "R" = "#dca96a",   
                                  "NR" = "#C45C69")),  
                              annotation_legend_param = list(title_gp = gpar(fontsize = (fontsize-2), fontface = "bold"), 
                                                             labels_gp = gpar(fontsize = (fontsize-2)), 
                                                             grid_height = unit(3.5, "mm"), 
                                                             grid_width = unit(3.5, "mm"))) 
range(dt) 
col_fun <- colorRamp2(seq(-2, 2, length.out = 3), c("#50859f", "white", "#d66692"))

# Plot
cellwidth = 0.1
cellheight = 0.225
cn = dim(as.matrix(dt))[2]
rn = dim(as.matrix(dt))[1]
w=cellwidth*cn
h=cellheight*rn

pdf("Heatmap_of_sel_halllmark_geneset_exp_in_TNR_avail_bulk_patients_260131.pdf", width = 9, height = 6)
Heatmap(
  matrix = dt,
  border = TRUE, 
  column_split = group3$split, 
  width = unit(w, "cm"), 
  height = unit(h, "cm"), 
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  clustering_method_rows = 'ward.D', 
  cluster_columns = F, 
  cluster_rows = T, 
  row_title_side = "left", 
  row_names_side = "left", 
  row_dend_side = "left", 
  show_column_names = F, 
  column_dend_height = unit(6, "mm"), 
  row_dend_width = unit(6, "mm"), 
  column_title_gp = gpar(fontsize = 10), 
  column_names_rot = 45, 
  column_names_gp = gpar(fontsize = 8), 
  row_names_gp = gpar(fontsize = 8),
  top_annotation = col_anno, 
  col = col_fun,
)
dev.off()

##################### Grouped boxplot of the key signatures (F1h) #####################
library(cowplot)
library(ggpubr)

dt_for_box <- as.data.frame(t(ssgsea))
dt_for_box$Sample <- rownames(dt_for_box)
dt_for_box <- merge(dt_for_box, group_sel, by = "Sample")
dt_for_box$Strat <- ifelse(dt_for_box$Treatment == "Pre" & dt_for_box$TNR >= 90, "R_Pre", 
                           ifelse(dt_for_box$Treatment == "Pre" & dt_for_box$TNR < 90, "NR_Pre", 
                                  ifelse(dt_for_box$Treatment == "Post" & dt_for_box$TNR >= 90, "R_Post", "NR_Post")))
dt_for_box$Strat <- factor(dt_for_box$Strat, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))
table(dt_for_box$Strat, dt_for_box$Treatment)

#write.table(dt_for_box, file='OS_182_Sigs_260802.txt', quote = F, sep = "\t")

pplist <- list()
for (group_ in names(geneSet1)) {
  cellper_ = dt_for_box %>% dplyr::select(one_of(c('Sample', 'Strat', group_)))
  colnames(cellper_)=c('Sample','Strat','percent')
  cellper_$percent =as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(Strat) %>% mutate(upper = quantile(percent,0.75),
                                                      lower = quantile(percent,0.25),
                                                      mean = mean(percent),
                                                      median = median(percent),
                                                      lower_lim = range(percent)[1], 
                                                      upper_lim = range(percent)[2]) 
  y_limits = c(min(cellper_$percent)*0.8, 
               max(cellper_$percent)*1.3) 
  y_breaks = round(seq(min(cellper_$percent)*0.8, 
                       max(cellper_$percent)*1.3, 
                       length.out = 5), 2) 
  p_position = max(cellper_$percent)*1.3 * 0.9
  p_position1 = max(cellper_$percent)*1.3 * 0.85
  p_position2 = max(cellper_$percent)*1.3 * 0.8
  
  mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
  pp1 = ggplot(cellper_,aes(x=Strat,y=percent)) + 
    geom_jitter(shape =21, aes(fill = Strat), width = 0.25) + 
    geom_boxplot(alpha = 0.5, size = 1, outlier.shape = NA, aes(fill = Strat)) + 
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) + 
    theme_cowplot() +
    theme(axis.text = element_text(size =10),axis.title = element_text(size = 10),legend.text = element_text(size =10),
          legend.title = element_text(size =10),plot.title = element_text(size =10,face ='plain'),legend.position ='none')+
    labs(title = group_, y = "Percentage", x = "Group") + 
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) + 
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
          axis.text.y = element_text(color = "black", size = 9),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) + 
    theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid"))
  
  labely = max(cellper_$percent)
  compare_means(percent ~ Strat,  data = cellper_)
  
  pp1 = pp1 + geom_signif(comparisons = list(c("R_Pre","R_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position, 
                          test.args = list(
                            alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("NR_Pre","NR_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))}, 
                          y_position = p_position, 
                          test.args = list(
                            alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("R_Pre","NR_Pre")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position1, 
                          test.args   = list(
                            alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black") 
  pp1 = pp1 + geom_signif(comparisons = list(c("R_Post","NR_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position2, 
                          test.args = "two.sided", 
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")
  
  pp1 <- pp1 + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  pplist[[group_]]= pp1
}

pps <- plot_grid(pplist[[1]], 
                 pplist[[2]],
                 pplist[[3]], 
                 pplist[[4]],
                 pplist[[5]],
                 pplist[[6]],
                 pplist[[7]],
                 pplist[[8]],
                 pplist[[9]],
                 pplist[[10]],
                 pplist[[11]], 
                 pplist[[12]],
                 pplist[[13]],
                 pplist[[14]],
                 pplist[[15]],
                 pplist[[16]], 
                 pplist[[17]], 
                 pplist[[18]],
                 pplist[[19]], 
                 pplist[[20]],
                 pplist[[21]], 
                 pplist[[22]], 
                 pplist[[23]],
                 pplist[[24]],
                 pplist[[25]],
                 pplist[[26]],
                 pplist[[27]],
                 pplist[[28]],
                 pplist[[29]],
                 pplist[[30]], 
                 ncol = 6 
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_Sigs", "_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 9, 
    height = 11.25) 
print(pps, newpage = FALSE)
dev.off()

##################### Supplementary UMAP plots (Extended Data F1a) #####################
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/All_cells_F1_main"))
library(Seurat)
library(tidyverse)

seuratObj <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))
DimPlot(seuratObj, group.by = "cluster")
names(seuratObj@meta.data) 

# DimPlot by Treatment status ----- 
Idents(seuratObj) <- "group"
seuratObj<- RenameIdents(seuratObj, "A" = "Pre", 
                         "B" = "Post",
                         "C" = "Pre",  
                         "D" = "Post")
seuratObj$treatment_status <- Idents(seuratObj)
colors <- c("#90486e","#efcfe3") 

p1 <- DimPlot(seuratObj, reduction = "umap", group.by = c("treatment_status"), label = T) + 
  scale_color_manual(values = colors) + coord_fixed(ratio = 0.9)
p1

pdf(paste0("All_cells", "_UMAP_by_treatment_stat_v260126.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

# DimPlot by Response ----- 
Idents(seuratObj) <- "group"
seuratObj<- RenameIdents(seuratObj, "A" = "R", 
                         "B" = "R",
                         "C" = "NR",  
                         "D" = "NR")
seuratObj$response <- Idents(seuratObj)
colors <- c("#daa2b3","#6e7b90") 

p1 <- DimPlot(seuratObj, reduction = "umap", group.by = c("response"), label = T) + 
  scale_color_manual(values = colors) + coord_fixed(ratio = 0.9) 
p1

pdf(paste0("All_cells", "_UMAP_by_response_stat_v260126.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

# DimPlot by Sample ----- 
colors <- c("#7E8FB3", "#5CA291", "#E07F59", "#CF7CAE", "#94C54C", "#E6C92B", "#D1B687", "#9E9E9E",
            "#9CAFD4", "#69C0AA", "#F19779", "#F29CD5", "#B5E464", "#FFF76C", "#F4D9B5", "#C7C7C7",
            "#18846A", "#C25401", "#6A659E", "#CB2880", "#5B861B", "#CCA401", "#8F6A1A",
            "#14745A", "#A34A01", "#5B5787", "#AD2067", "#4E7916", "#A88401", "#755517")

p1 <- DimPlot(seuratObj, reduction = "umap", group.by = c("sample"), label = T) + 
  scale_color_manual(values = colors) + coord_fixed(ratio = 0.9) #+ 
p1

pdf(paste0("All_cells", "_UMAP_by_sample_v260126.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

# Unintegrated UMAPs ------- 
seuratObj <- NormalizeData(seuratObj) %>% FindVariableFeatures() %>% ScaleData()
seuratObj <- RunPCA(seuratObj)
seuratObj <- FindNeighbors(seuratObj, dims = 1:50, reduction = "pca")
seuratObj <- RunUMAP(seuratObj, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(seuratObj, reduction = "umap.unintegrated", group.by = c("group"))

## Repeat run but with unintegrated UMAP ----- 
### By group ----- 
Idents(seuratObj) <- "group"
seuratObj<- RenameIdents(seuratObj, "A" = "Pre", 
                         "B" = "Post",
                         "C" = "Pre",  
                         "D" = "Post")
seuratObj$treatment_status <- Idents(seuratObj)
colors <- c("#90486e","#efcfe3")

p1 <- DimPlot(seuratObj, reduction = "umap.unintegrated", group.by = c("treatment_status"), label = T) + 
  scale_color_manual(values = colors) + coord_fixed(ratio = 0.9) #+ 
p1

pdf(paste0("All_cells", "_Uninte_UMAP_by_treatment_stat_v260201.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

# DimPlot by Response ----- 
Idents(seuratObj) <- "group"
seuratObj<- RenameIdents(seuratObj, "A" = "R", 
                         "B" = "R",
                         "C" = "NR",  
                         "D" = "NR")
seuratObj$response <- Idents(seuratObj)
colors <- c("#daa2b3","#6e7b90")  

p1 <- DimPlot(seuratObj, reduction = "umap.unintegrated", group.by = c("response"), label = T) + 
  scale_color_manual(values = colors) + coord_fixed(ratio = 0.9) 
p1

pdf(paste0("All_cells", "_unint_UMAP_by_response_stat_v260201.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

# DimPlot by Sample ----- 
colors <- c("#7E8FB3", "#5CA291", "#E07F59", "#CF7CAE", "#94C54C", "#E6C92B", "#D1B687", "#9E9E9E",
            "#9CAFD4", "#69C0AA", "#F19779", "#F29CD5", "#B5E464", "#FFF76C", "#F4D9B5", "#C7C7C7",
            "#18846A", "#C25401", "#6A659E", "#CB2880", "#5B861B", "#CCA401", "#8F6A1A",
            "#14745A", "#A34A01", "#5B5787", "#AD2067", "#4E7916", "#A88401", "#755517")

p1 <- DimPlot(seuratObj, reduction = "umap.unintegrated", group.by = c("sample"), label = T) + 
  scale_color_manual(values = colors) + coord_fixed(ratio = 0.9) 
p1

pdf(paste0("All_cells", "_unintUMAP_by_sample_v260201.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

##################### nCount and nGenes (Extended Data F1b) #####################
names(seuratObj@meta.data)
mat <- VlnPlot(seuratObj, "nCount_RNA", group.by = "sample", pt.size = 0)
mat_dat <- mat@data 
colors <- c("#7E8FB3", "#5CA291", "#E07F59", "#CF7CAE", "#94C54C", "#E6C92B", "#D1B687", "#9E9E9E",
            "#9CAFD4", "#69C0AA", "#F19779", "#F29CD5", "#B5E464", "#FFF76C", "#F4D9B5", "#C7C7C7",
            "#18846A", "#C25401", "#6A659E", "#CB2880", "#5B861B", "#CCA401", "#8F6A1A",
            "#14745A", "#A34A01", "#5B5787", "#AD2067", "#4E7916", "#A88401", "#755517")

p1 <- ggplot(mat_dat, aes(x=ident, y=nCount_RNA)) +
  geom_boxplot(aes(color=ident), fill=NA, outlier.size = 0) +  
  scale_color_manual(values=colors) + 
  theme_classic() + 
  theme(axis.text = element_text(color = 'black')) + 
  labs(title="", x="", y="nCount_RNA") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p1 

# nFeature 
names(seuratObj@meta.data)
mat <- VlnPlot(seuratObj, "nFeature_RNA", group.by = "sample", pt.size = 0)
mat_dat <- mat@data 

p2 <- ggplot(mat_dat, aes(x=ident, y=nFeature_RNA)) +
  geom_boxplot(aes(color=ident), fill=NA, outlier.size = 0) +  
  scale_color_manual(values=colors) + 
  theme_classic() + 
  theme(axis.text = element_text(color = 'black')) + 
  labs(title="", x="", y="nFeature_RNA") + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p2

p3 <- p1/ p2

pdf('nConut_nFeatures_per_sample.pdf', width = 4.75, height =  8.5, family = "Helvetica")
print(p3)
dev.off()

##################### ROGUE (Extended Data F1c) #####################
# ROGUE by cell type ----- 
library(ROGUE)
rogue_res_mat <- as.data.frame(matrix(
  nrow = length(unique(seuratObj$sample)), ncol = length(unique(seuratObj$cluster)), NA))

colnames(rogue_res_mat) <- unique(seuratObj$cluster)
rownames(rogue_res_mat) <- unique(seuratObj$sample) 
Idents(seuratObj) <- "cluster"

# Loop for subsets 
for (j in levels(seuratObj$cluster)) {
  print(paste0("Run for ", j))
  cells <- WhichCells(seuratObj, idents = j)
  scrna <- subset(seuratObj, cells = cells)
  
  scrna <- NormalizeData(scrna) %>% FindVariableFeatures() %>% ScaleData() 
  expr <- as.matrix(GetAssayData(object = scrna[["RNA"]], slot = "data")) 
  meta <- scrna@meta.data 
  expr <- matr.filter(expr, min.cells = 10, min.genes = 10) 
  ent.res <- SE_fun(expr) 
  head(ent.res) 
  SEplot(ent.res) 
  rogue.value <- CalculateRogue(ent.res, platform = "UMI") 
  rogue.value 
  
  meta$sample <- droplevels(meta$sample) 
  meta$cluster <- droplevels(meta$cluster) 
  rogue.res <- rogue(expr, labels = meta$group, samples = meta$sample, platform = "UMI", span = 0.6)
  saveRDS(rogue.res, paste0('rogue_res_', j,'.rds')) 
  
  val <- vector()
  for (i in unique(seuratObj$sample)) {
    if (FALSE %in% is.na(rogue.res[i, ])) {
      sample_rogue <- rogue.res[i, ][which(!is.na(rogue.res[i, ]))] %>% as.numeric()
    } else {
      sample_rogue <- NA 
    }
    val <- append(val, sample_rogue)
  }
  rogue_res_mat[, j] <- val
}

saveRDS(rogue_res_mat, "Per_sample_mat_for_ROGUE_by_cluster_260126.rds")

rogue_res_mat
mat <- t(rogue_res_mat) %>% as.data.frame()
mat$cluster <- rownames(mat)
vars <- colnames(mat)[!colnames(mat) %in% "cluster"]
rogue_res_mat_long <-melt(mat, id.vars = c("cluster"), 
                          measure.vars = vars,
                          variable.name = c('sample'),
                          value.name = 'value')

mean_mat <- as.data.frame(matrix(ncol = 1, nrow = length(levels(seuratObj))))
rownames(mean_mat) <- levels(seuratObj)
colnames(mean_mat) <- "Mean ROGUE"

for (k in unique(rownames(mean_mat))) {
  temp_mean <- median(as.numeric(mat[k, vars]), na.rm = T)
  mean_mat[k, ] <- temp_mean
}

levels <- rownames(mean_mat)[order(mean_mat$`Mean ROGUE`, decreasing = T)]
colors <- c('#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',
            '#CC6677', '#882255', '#AA4499', '#662288', '#BB5566',
            '#5D3A9B', '#997700', '#DDDDDD', '#004488',
            '#6699CC', '#66BBAA') 
names(colors) <- levels(seuratObj$cluster)

rogue_res_mat_long$cluster <- factor(rogue_res_mat_long$cluster, levels = levels)
p1 <- ggplot(rogue_res_mat_long, aes(x=cluster, y=value)) +
  geom_boxplot(aes(color=cluster), fill=NA, outlier.size = 0) +  
  geom_jitter(aes(color=cluster), width=0.2, size=0.5) +  
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_color_manual(values=colors) + 
  theme_classic() + 
  theme(axis.text = element_text(color = 'black')) + 
  labs(title="", x="", y="ROGUE index") + 
  #scale_y_continuous(limits = c(0.75, 0.95), breaks = seq(0.75, 0.95, by = 0.1)) +  # 设置 y 轴范围和间隔
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
p1 

pdf('Rogue_boxplot_per_cell_group.pdf', width = 4, height = 3, family = "Helvetica")
print(p1)
dev.off()

##################### LISI (Extended Data F1d) #####################
#install.packages("devtools")
#devtools::install_github("immunogenomics/lisi")
library(lisi)

# Unint 
embed <- seuratObj@reductions$umap.unintegrated@cell.embeddings
head(embed)

meta <- data.frame(seuratObj$sample)
meta$patient <- meta$seuratObj.sample
meta$patient <- gsub("_AF_", "", meta$patient)
meta$patient <- gsub("_BF_", "", meta$patient)
head(meta)

res <- compute_lisi(embed, meta, c('seuratObj.sample', 'patient'))
boxplot(res)
mat_for_plot <- res
mat_for_plot$UMAP <- "Unint"

# Int 
embed <- seuratObj@reductions$umap@cell.embeddings
res <- compute_lisi(embed, meta, c('seuratObj.sample', 'patient'))
res$UMAP <- "Int"
mat_for_plot_comb <- rbind(mat_for_plot, res)
head(mat_for_plot_comb)
mat_for_plot_comb$UMAP <- factor(mat_for_plot_comb$UMAP, levels = rev(c("Int", "Unint")))

# Plot
library(ggplot2)

mycol <- c("#6f6f6f", "#6d65a3")
p1 <- ggplot(mat_for_plot_comb,aes(x = UMAP, y = patient, fill = UMAP)) + 
  geom_boxplot(alpha = 0.5, size = 1, outlier.shape = NA, aes(fill = UMAP)) +
  scale_fill_manual(values = mycol) +  
  scale_y_continuous(limits = c(0, 15), breaks = c(0, 5, 10, 15), expand = c(0, 0)) + 
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) + 
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid"))
p1

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_LISI_Unint_int_comparison_260201.pdf"), 
    width = 3, 
    height = 6) 
print(p1, newpage = FALSE)
dev.off()

##################### MuSiC (Extended data fig 1g,h) ##################### 
## MuSiC ------ 
#devtools::install_github('xuranw/MuSiC')
library(MuSiC)
library(SingleCellExperiment)
library(Biobase)
library(tidyverse)

# Bulk expression matrix
## Read the bulk files
### fpkm w/o integration/batch effect removal 
fpkm_merge <- read.csv(here("data/EcoTyper_run_16_06_25/merged_bulk_data_no_batch_correc_260731.csv"), header = T) # Ready for use for MuSiC
fpkm_merge = fpkm_merge[!duplicated(fpkm_merge$GeneSymbol),]
fpkm_merge <- fpkm_merge[, -1]

rownames(fpkm_merge) <- fpkm_merge[, 1]
bulk.mtx = fpkm_merge[, -1]

object <- new("ExpressionSet", exprs=as.matrix(bulk.mtx))
bulk <- exprs(object)

# Deconv w/ major CD45+ and CD45- cell populations for lower RNA-seq exp matrix resolution
seuratObj <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))
Idents(seuratObj) <- "cluster"
seuratObj <- RenameIdents(seuratObj, "BCells" = "CD45+_immune", 
                          "Chondrocytes" = "Osteogenic", 
                          "ECs" = "Endo", 
                          "Erythrocytes" = "Hematopoietic", 
                          "Fibroblasts" = "Stromal", 
                          "HSCs" = "Hematopoietic", 
                          "Myeloid" = "CD45+_immune", 
                          "MSCs" = "Stromal", 
                          "MastCells" = "CD45+_immune", 
                          "MuralCells" = "Stromal", 
                          "Osteoblasts" = "Osteogenic", 
                          "Osteoclasts" = "Osteoclastic", 
                          "PlasmaCells" = "CD45+_immune", 
                          "Pro_B" = "CD45+_immune", 
                          "TandNK" = "CD45+_immune",
                          "pDCs" = "CD45+_immune"
)

seuratObj$cluster <- Idents(seuratObj)
OS.sce <- as.SingleCellExperiment(seuratObj)
rm(seuratObj)
gc()

# Estimate cell type proportions
length(intersect(rownames(bulk.mtx), rownames(OS.sce@assays@data$counts)))
rownames(bulk) <- str_replace(rownames(bulk), "-", ".")
rownames(OS.sce@assays@data$counts) <- str_replace(rownames(OS.sce@assays@data$counts), "-", ".")
length(intersect(rownames(bulk.mtx), rownames(OS.sce@assays@data$counts)))

Est.prop.OS.inhouse = music_prop(bulk.mtx = bulk, 
                                 #sc.eset = OS.sce, 
                                 sc.sce = OS.sce,
                                 clusters = 'cluster',
                                 samples = 'sample', 
                                 #select.ct = c("Osteoblasts", "MSCs", "Fibroblasts", "MuralCells", 
                                 #               "ECs", "Osteoclasts", "Myeloid", "TandNK"), 
                                 select.ct = levels(OS.sce$cluster),
                                 verbose = T)

## Plot ------ 
library(ComplexHeatmap)
library(circlize)

group <- read.csv(here("data/Re-clustering_fine_cell_types_21_05_25/Metadata_260802/updated_bulk_metadta_260802.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
colnames(group)[3] <- "Sample"

group$TNR <- as.numeric(group$TNR)
group_sel <- group[!is.na(group$TNR), ]

results_dt <- Est.prop.OS.inhouse$Est.prop.weighted %>% t()
results_dt <- results_dt[, group_sel$Sample]
identical(group_sel$Sample, colnames(results_dt)) # TRUE

group1 <- group_sel[, c("Sample","Treatment", "Progression", "TNR")]
group1$Progression <- ifelse(group1$Progression == "1", "Progression", "Non-progression")
group1$TNR <- as.numeric(group1$TNR)
group1 <- group1[order(group1$TNR, decreasing = T), ]
group1 <- group1[order(group1$Progression), ]
group1 <- group1[order(group1$Treatment, decreasing = T), ]

results_dt <- results_dt[, group1$Sample]
#identical(group_sel$Sample, colnames(dt)) # FALSE
identical(group1$Sample, colnames(results_dt)) # TRUE

group2 <- group1[, c(2, 3, 4)]
identical(group1$Sample, colnames(results_dt)) # Ensure True 
identical(group2$Treatment, group1$Treatment) # Ensure True 
group2$TNR <- as.numeric(group2$TNR)
group2$TNR_cat <- ifelse(group2$TNR >= 90, "R", "NR")

group3 <- group2
group3$split <- NA
group3$split[group3$TNR_cat == "R" & group3$Treatment == "Pre"] <- "A"
group3$split[group3$TNR_cat == "NR" & group3$Treatment == "Pre"] <- "B"
group3$split[group3$TNR_cat == "R" & group3$Treatment == "Post"] <- "C"
group3$split[group3$TNR_cat == "NR" & group3$Treatment == "Post"] <- "D"
group3$split <- factor(group3$split, levels = c("A", "C", "B", "D"))

fontsize = 8
col_anno <- HeatmapAnnotation(df = group2,
                              show_annotation_name = T,
                              gp = gpar(fontsize = 10),
                              na_col = "black", 
                              col = list(Treatment = c(
                                'Pre' = '#90486e', 
                                'Post' = '#efcfe3'), 
                                Progression = c(
                                  "Non-progression" = "#bfbebe", 
                                  "Progression" = "#010101"), 
                                TNR = colorRamp2(c(0, 20, 40, 60, 80, 90, 100), 
                                                 c("#eeecdf","#eeecdf", "#becdd2", "#6f9ad1", "#44679f", "#3f4f71", "#3f4f71")),
                                TNR_cat = c(
                                  "R" = "#dca96a", 
                                  "NR" = "#c45c69")),
                              annotation_legend_param = list(title_gp = gpar(fontsize = (fontsize-2), fontface = "bold"),
                                                             labels_gp = gpar(fontsize = (fontsize-2)),
                                                             grid_height = unit(3.5, "mm"),
                                                             grid_width = unit(3.5, "mm")))
col_fun <- colorRamp2(seq(0, 1, length.out = 10), c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                                    "#d73027", "#a50026"))


# Plot
cellwidth = 0.1
cellheight = 0.5
cn = dim(as.matrix(results_dt))[2]
rn = dim(as.matrix(results_dt))[1]
w=cellwidth*cn
h=cellheight*rn

levels(OS.sce$cluster)
row.order.for.plot <- c("Osteogenic", "Stromal", "Endo", "Hematopoietic", "Osteoclastic", "CD45+_immune")

pdf("Heatmap_of_meta_celltype_in_TNR_avail_bulk_patients_260204.pdf", width = 10, height = 10)
Heatmap(
  matrix = results_dt,
  border = TRUE, 
  row_order = row.order.for.plot, 
  column_split = group3$split, 
  width = unit(w, "cm"), 
  height = unit(h, "cm"), 
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  clustering_method_rows = 'ward.D', 
  cluster_columns = F, 
  cluster_rows = F, 
  row_title_side = "left", 
  row_names_side = "left", 
  row_dend_side = "left", 
  show_column_names = F, 
  column_dend_height = unit(6, "mm"), 
  row_dend_width = unit(6, "mm"), 
  column_title_gp = gpar(fontsize = 10), 
  column_names_rot = 45, 
  column_names_gp = gpar(fontsize = 8), 
  row_names_gp = gpar(fontsize = 8),
  top_annotation = col_anno, 
  col = col_fun,
)
dev.off()

#### Boxplot ------
library(cowplot)
library(ggpubr)

dt_for_box <- as.data.frame(t(results_dt))
dt_for_box$Sample <- rownames(dt_for_box)
dt_for_box <- merge(dt_for_box, group_sel, by = "Sample")
dt_for_box$Strat <- ifelse(dt_for_box$Treatment == "Pre" & dt_for_box$TNR >= 90, "R_Pre", 
                           ifelse(dt_for_box$Treatment == "Pre" & dt_for_box$TNR < 90, "NR_Pre", 
                                  ifelse(dt_for_box$Treatment == "Post" & dt_for_box$TNR >= 90, "R_Post", "NR_Post")))
dt_for_box$Strat <- factor(dt_for_box$Strat, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))
table(dt_for_box$Strat, dt_for_box$Treatment)

#write.table(dt_for_box, file='OS_MUSIC_Prop_260802.txt', quote = F, sep = "\t")

pplist <- list()
for (group_ in levels(OS.sce$cluster)) {
  cellper_ = dt_for_box %>% dplyr::select(one_of(c('Sample', 'Strat', group_)))
  colnames(cellper_)=c('Sample','Strat','percent')
  cellper_$percent =as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(Strat) %>% mutate(upper = quantile(percent,0.75),
                                                      lower = quantile(percent,0.25),
                                                      mean = mean(percent),
                                                      median = median(percent),
                                                      lower_lim = range(percent)[1], 
                                                      upper_lim = range(percent)[2])
  y_limits = c(min(cellper_$percent)*0.8, 
               max(cellper_$percent)*1.3) 
  y_breaks = round(seq(min(cellper_$percent)*0.8, 
                       max(cellper_$percent)*1.3, 
                       length.out = 5), 2) 
  p_position = max(cellper_$percent)*1.3 * 0.9
  p_position1 = max(cellper_$percent)*1.3 * 0.85
  p_position2 = max(cellper_$percent)*1.3 * 0.8
  
  mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
  pp1 = ggplot(cellper_,aes(x=Strat,y=percent)) + #ggplot作图
    geom_jitter(shape =21, aes(fill = Strat), width = 0.25) + 
    geom_boxplot(alpha = 0.5, size = 1, outlier.shape = NA, aes(fill = Strat)) + 
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) + 
    theme_cowplot() +
    theme(axis.text = element_text(size =10),axis.title = element_text(size = 10),legend.text = element_text(size =10),
          legend.title = element_text(size =10),plot.title = element_text(size =10,face ='plain'),legend.position ='none')+
    labs(title = group_, y = "Percentage", x = "Group") + 
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) + 
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
          axis.text.y = element_text(color = "black", size = 9),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) + 
    theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid"))
  
  labely = max(cellper_$percent)
  compare_means(percent ~ Strat,  data = cellper_)
  
  pp1 = pp1 + geom_signif(comparisons = list(c("R_Pre","R_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position, 
                          test.args = list(
                            alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("NR_Pre","NR_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))}, 
                          y_position = p_position, 
                          test.args = list(
                            alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("R_Pre","NR_Pre")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position1, 
                          test.args   = list(
                            alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black") 
  pp1 = pp1 + geom_signif(comparisons = list(c("R_Post","NR_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position2, 
                          test.args = "two.sided", 
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")
  
  pp1 <- pp1 + 
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  pplist[[group_]]= pp1
}

pps <- plot_grid(pplist[[1]], 
                 pplist[[2]],
                 pplist[[3]], 
                 pplist[[4]],
                 pplist[[5]],
                 pplist[[6]],
                 ncol = 3 
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_Bulk_subsets_all_major_celltype", "_co-occur_plot.pdf"), 
    width = 4.5, 
    height = 4.5) 
print(pps, newpage = FALSE)
dev.off()

##################### Pi and Ti calculations (Extended data fig 1i) #####################
### Pi ----- 
library(reshape2)
Idents(seuratObj) <- "reportname"
seuratObj <- RenameIdents(seuratObj,
                          "SAH_001" = 99.6, 
                          "SAH_002" = 99.8, 
                          "SAH_003" = 100, 
                          "SAH_004" = 100,
                          "SAH_005" = 99.7,
                          "SAH_006" = 99.5, 
                          "SAH_007" = 98.9,
                          "SAH_008" = 94.1,
                          "SAH_009" = 99.6, 
                          "SAH_010" = 99.8, 
                          "SAH_011" = 100, 
                          "SAH_012" = 100,
                          "SAH_013" = 99.7,
                          "SAH_014" = 99.5,
                          "SAH_015" = 98.9, 
                          "SAH_016" = 94.1, 
                          "SAH_017" = 70.1, 
                          "SAH_018" = 47.8, 
                          "SAH_019" = 81.4, 
                          "SAH_020" = 85.1,
                          "SAH_021" = 76.2, 
                          "SAH_022" = 68.6, 
                          "SAH_023" = 80.5, 
                          "SAH_024" = 70.1, 
                          "SAH_025" = 47.8, 
                          "SAH_026" = 81.4, 
                          "SAH_027" = 85.1, 
                          "SAH_028" = 76.2, 
                          "SAH_029" = 68.6, 
                          "SAH_030" = 80.5)
seuratObj$TNR_val <- Idents(seuratObj)

# Pi 
Idents(seuratObj) <- "treatment_status"
scrna_pre <- subset(seuratObj, idents = "Pre")
Cellratio <- prop.table(table(scrna_pre$cluster, scrna_pre$sample), margin = 2) 
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio, Var2~Var1, value.var="Freq")
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]
cellper <- na.omit(cellper)

meta <- scrna_pre@meta.data
colnames(meta)
meta <- meta[,c(4, 20)] # Selecting columns withs sample and TNR information 
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
meta <- unique(meta) 
cellper <- merge(cellper, meta, by = "sample")

colnames(cellper)
cellper <- as.data.frame(cellper)

# Calc 
library(ggplot2)
library(dplyr)
cellper$TNR_val <- as.numeric(as.character(cellper$TNR_val))

# Create the matrix for plot
matrix_for_Plots <- matrix(data = rep(NA, 3*length(unique(seuratObj$cluster))), 
                           ncol = 3, nrow = length(unique(seuratObj$cluster)), dimnames = list(c(as.vector(unique(seuratObj$cluster))), c("Cell_type", "Pi", "pval"))) %>% as.data.frame()
matrix_for_Plots$Cell_type <- rownames(matrix_for_Plots)

for (i in as.vector(unique(seuratObj$cluster))) {
  m <- lm(unlist(cellper["TNR_val"]) ~ unlist(cellper[i]))
  b <- as.numeric(round(unname(coef(m)[2]), digits = 2))
  r2 = as.numeric(round(summary(m)$r.squared, digits = 2))
  Pi <- b/(abs(b)) * r2
  pval <- summary(m)$coefficients
  pval <- as.numeric(round(pval[2, 4], digits = 3))
  matrix_for_Plots[i, "Pi"] <- Pi
  matrix_for_Plots[i, "pval"] <- pval 
  matrix_for_Plots$log10pval <- -log10(matrix_for_Plots$pval)
}

matrix_for_Plots <- matrix_for_Plots %>%
  arrange(Pi) %>%
  mutate(Cell_type = factor(Cell_type,levels = Cell_type)) 
matrix_for_Plots <- na.omit(matrix_for_Plots)

# Plot
mycol <- c('#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',
           '#CC6677', '#882255', '#AA4499', '#662288', '#BB5566',
           '#5D3A9B', '#997700', '#DDDDDD', '#004488',
           '#6699CC', '#66BBAA')  
names(mycol) <- levels(seuratObj$cluster)
mycol <- mycol[levels(matrix_for_Plots$Cell_type)]

p <- ggplot(matrix_for_Plots, aes(x = Cell_type, y = Pi)) +  
  geom_segment(aes(x = Cell_type, xend = Cell_type, y = 0, yend = Pi, color = Cell_type),                 
               linetype = "solid", size = 1, color = mycol) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour="gray80") +  
  geom_point(aes(color = Cell_type, size = log10pval), color = mycol
  ) + 
  scale_size_continuous(range=c(1,9)) +
  geom_text(aes(label = Cell_type, y = 0), 
            hjust = ifelse(matrix_for_Plots$Pi >= 0, 1, 0), 
            angle = 90, 
            fontface = 'italic',
            color = ifelse(matrix_for_Plots$pval >= 0.05, "black", 'red'),             
            size = 4) + 
  theme_classic() +
  theme(axis.text = element_text(size =10),axis.title = element_text(size = 10),legend.text = element_text(size =10),
        legend.title = element_text(size =10),plot.title = element_text(size =10,face ='plain'),legend.position ='none')+
  labs(title = "Predicative Index", y = "Predicative Index") + 
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18, hjust = 0.5),
        legend.position = "none") + 
  NoLegend() + 
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid")) 
p 

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Pi_all_cells_for_fig_T.pdf"), 
    width = 3.5, 
    height = 5) 
print(p, newpage = FALSE)
dev.off()

###### Ti #####
# Needs cell per change 
library(reshape2)
# Calc post
Idents(seuratObj) <- "treatment_status"
scrna_post <- subset(seuratObj, idents = "Post")
Cellratio_post <- prop.table(table(scrna_post$cluster, scrna_post$sample), margin =2)
Cellratio_post <- data.frame(Cellratio_post)

cellper_post <- dcast(Cellratio_post, Var2~Var1, value.var="Freq")
rownames(cellper_post)<- cellper_post[,1]
cellper_post <- cellper_post[,-1]
cellper_post <- na.omit(cellper_post)

# Calc pre
Cellratio_pre <- prop.table(table(scrna_pre$cluster, scrna_pre$sample), margin =2)
Cellratio_pre <- data.frame(Cellratio_pre)

cellper_pre <- dcast(Cellratio_pre, Var2~Var1, value.var="Freq")
rownames(cellper_pre)<- cellper_pre[,1]
cellper_pre <- cellper_pre[,-1]
cellper_pre <- na.omit(cellper_pre)

cellper_delta <- cellper_post - cellper_pre

meta <- scrna_post@meta.data
colnames(meta)
meta <- meta[,c(4,20)] # Selecting columns withs sample and TNR information 
meta <- as.data.frame(meta)
cellper_delta$sample <- rownames(cellper_delta)
meta <- unique(meta)
cellper_delta <- merge(cellper_delta, meta, by = "sample")
colnames(cellper_delta)
cellper_delta <- as.data.frame(cellper_delta)

# Calc 
library(ggplot2)
library(dplyr)
cellper_delta$TNR_val <- as.numeric(as.character(cellper_delta$TNR_val))

# Create the matrix for plot
matrix_for_Plots_delta <- matrix(data = rep(NA, 3*length(unique(seuratObj$cluster))), 
                                 ncol = 3, nrow = length(unique(seuratObj$cluster)), dimnames = list(c(as.vector(unique(seuratObj$cluster))), c("Cell_type", "Ti", "pval"))) %>% as.data.frame()
matrix_for_Plots_delta$Cell_type <- rownames(matrix_for_Plots_delta)

for (i in as.vector(unique(seuratObj$cluster))) {
  m <- lm(unlist(cellper_delta["TNR_val"]) ~ unlist(cellper_delta[i]))
  b <- as.numeric(round(unname(coef(m)[2]), digits = 2))
  r2 = as.numeric(round(summary(m)$r.squared, digits = 2))
  Ti <- b/(abs(b)) * r2
  pval <- summary(m)$coefficients
  pval <- as.numeric(round(pval[2, 4], digits = 3))
  matrix_for_Plots_delta[i, "Ti"] <- Ti
  matrix_for_Plots_delta[i, "pval"] <- pval 
  matrix_for_Plots_delta$log10pval <- -log10(matrix_for_Plots_delta$pval)
}

matrix_for_Plots_delta <- matrix_for_Plots_delta %>%
  arrange(Ti) %>%
  mutate(Cell_type = factor(Cell_type,levels = Cell_type)) 

matrix_for_Plots_delta <- na.omit(matrix_for_Plots_delta)

mycol1 <- as.data.frame(mycol); mycol1$celltype <- rownames(mycol1); 
mycol1 <- mycol1[match(matrix_for_Plots_delta$Cell_type, mycol1$celltype), ]       
mycol <- mycol1$mycol

p1 <- ggplot(matrix_for_Plots_delta, aes(x = Cell_type, y = Ti)) +  
  geom_segment(aes(x = Cell_type, xend = Cell_type, y = 0, yend = Ti, color = Cell_type),                 
               linetype = "solid", size = 1, color = mycol) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour="gray80") +  
  geom_point(aes(color = Cell_type, size = log10pval), color = mycol
  ) + 
  scale_size_continuous(range=c(1,9)) +
  geom_text(aes(label = Cell_type, y = 0), 
            hjust = ifelse(matrix_for_Plots_delta$Ti >= 0, 1, 0), 
            #vjust = -0.5, 
            angle = 90, 
            fontface = 'italic',
            color = ifelse(matrix_for_Plots_delta$pval >= 0.05, "black", 'red'),             
            size = 4) + 
  theme_classic() +
  theme(axis.text = element_text(size =10),axis.title = element_text(size = 10),legend.text = element_text(size =10),
        legend.title = element_text(size =10),plot.title = element_text(size =10,face ='plain'),legend.position ='none')+
  labs(title = "Therapeutic Index", y = "Therapeutic Index") + 
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18, hjust = 0.5),
        legend.position = "none") + 
  NoLegend() + 
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid")) 
scaleFUN <- function(x) sprintf("%.2f", x)
p1 <- p1 + scale_y_continuous(labels=scaleFUN)
p1

library(cowplot)
p2 <-plot_grid(p, p1, ncol = 1)
p2

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Pi_Ti_all_cells_for_fig_1_supp.pdf"), 
    width = 6, 
    height = 6) 
print(p2, newpage = FALSE)
dev.off()
