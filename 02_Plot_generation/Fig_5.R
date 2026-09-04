library(here) # project-root-relative paths; run scripts from repository root
##################### Visualisation of PySCENIC result (related to Fig 5G) #####################
rm(list = ls())
# disabled, run from repo root: setwd(here("data/PySCENIC_10_03_25/run1_myeloid_10_03_25/260125_TNK"))
library(Seurat)
#devtools::install_github("aertslab/SCopeLoomR")
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2)
library(stringr)
library(circlize) 

sce <- readRDS(
  here("data/Correlation_and_crosstalk_analysis_18_07_25/All_cells_for_ecoTyper_cell_type_adjusted_250905.rds"))
Idents(sce) <- "group"
sce <- RenameIdents(sce, "A" = "H_Pre", "B" = "H_Post", "C" = "L_Pre", "D" = "L_Post")
sce$group_anno <- Idents(sce)

Idents(sce) <- "Cell_type_for_ecotyper"
sce <- subset(sce, idents = c("Th", "Tc", "NK", "B"))
Idents(sce) <- "Cell_type_fine_harmony"

loom <- open_loom(here("data/PySCENIC_10_03_25/run1_myeloid_10_03_25/260125_TNK/aucell.loom"))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat) 
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)
embeddings

sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
identical(colnames(sub_regulonAUC), colnames(sce)) # Ensure true (identical)

# pass on cell type information 
cellTypes <- data.frame(row.names = colnames(sce),                         
                        celltype = sce$Cell_type_fine_harmony)
head(cellTypes)
sub_regulonAUC[1:4,1:4] 

# Partition based on own criteria
selectedResolution <- "celltype" # 
cellsPerGroup <- split(rownames(cellTypes),cellTypes[,selectedResolution]) 
regulonsub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 

###### Activity of all regulons #####
regulonActivity_byGroup <- sapply(cellsPerGroup, 
                                  function(x) 
                                    rowMeans(getAUC(sub_regulonAUC)[,x])) 
regulonActivity_byGroup <- regulonActivity_byGroup[, is.finite(colSums(regulonActivity_byGroup))] 
colnames(regulonActivity_byGroup)
#[1] "NK_KIT+"           "NK_FCGR3A+"        "NK_KLRC1+"         "CD8_MAIT"          "Tmem_IL7R+"        "Tmem_GPR183+FOSB+"
#[7] "CD8_Teff_IFNG+"    "CD8_Teff_GZMK+"    "Tex_CXCL13+"       "Tn_CCR7+"          "T_Ki67+"           "Trm_ZNF683+"      
#[13] "T_Mtgenes_hi"      "gdT"               "Treg"              "BCells"            "PlasmaCells"       "Pro_B"   

cellTypes_to_keep <- c()
regulonActivity_byGroup <- regulonActivity_byGroup[, c("Tn_CCR7+", "Tmem_IL7R+", "Tmem_GPR183+FOSB+", "Treg", 
                                                       "CD8_Teff_IFNG+", "CD8_Teff_GZMK+", "Tex_CXCL13+", 
                                                       "NK_KIT+", "NK_KLRC1+", "NK_FCGR3A+"
)]               

regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 

regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

##### DEGs of regulon activity #####
library(dplyr) 
rss=regulonActivity_byGroup_Scaled
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss), 
                 cluster = colnames(rss)[i], 
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2

top5 <- df %>% 
  group_by(cluster) %>% 
  top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 

library(ComplexHeatmap)
heatmap_data_scaled <- n[unique(rownames(n)), ]
cellwidth = 0.5
cellheight = 0.5
cn = dim(as.matrix(heatmap_data_scaled))[2]
rn = dim(as.matrix(heatmap_data_scaled))[1]
w=cellwidth*cn
h=cellheight*rn

pdf(file = "T_NK_cell_type_of_selection_regulon_act_260125.pdf", width = 5, height = 15)
Heatmap(as.matrix(heatmap_data_scaled),
        # Height of the blocks
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        #scale = F,
        #rect_gp = gpar(col = "white", lwd = 1.5),
        #border_g = gpar(col = ,lty = 1,lwd = 1.2),
        # Dent formatting 
        column_dend_height = unit(1, "cm"), 
        row_dend_width = unit(1, "cm"),
        #column_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
        #row_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
        # Gap numbers and style 
        #row_split = 6, column_split = c("1", "1", "2", "2"),
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        # Text formats 
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        show_heatmap_legend = T, 
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = 'left',
        column_names_side = c("top"), 
        show_column_names = T,
        show_row_names = T,
        border = T, 
        col = colorRamp2(breaks = c(-1.5, 0, 1.5), colors = c("#50859f","white","#d66692")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if( heatmap_data_scaled[i, j] > 1) {
            grid.text("+++", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] > 0.5 & heatmap_data_scaled[i, j] <= 1) {
            grid.text("++", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= 0.5 & heatmap_data_scaled[i, j] > 0.2) {
            grid.text("+", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= 0.2 & heatmap_data_scaled[i, j] > -0.2) {
            grid.text("+/-", x, y, gp = gpar(fontsize = 8))
          } #else if( heatmap_data_scaled[i, j] <= -0.2 & heatmap_data_scaled[i, j] > -0.5) {
          #grid.text("-", x, y, gp = gpar(fontsize = 8))
          #} else if( heatmap_data_scaled[i, j] <= -0.5 & heatmap_data_scaled[i, j] > -1) {
          #  grid.text("--", x, y, gp = gpar(fontsize = 8))
          #} else if( heatmap_data_scaled[i, j] <= -1 & heatmap_data_scaled[i, j] > -2) {
          #  grid.text("---", x, y, gp = gpar(fontsize = 8))
          #}
        }
)
dev.off()

