library(here) # project-root-relative paths; run scripts from repository root
##################### cDC2 pipeline #####################
## CD83+/- cDC1 demarcations ----- 
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/DC_251025"))
library(Seurat)
library(tidyverse)
#scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Myeloid/All_myeloid_reident_updated_27_05_25.rds"))
#scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Myeloid/Myeloid_final_as_at_24_07.rds"))
#Wlevels(scrna)

## pipeline ----- 
seuratObj <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))
cDC2 <- subset(seuratObj, idents = c("cDC2_CD1C+"))
cDC2[["RNA"]]$scale.data <- NULL
cDC2[["RNA"]]$data <- NULL
cDC2[["RNA"]] <- split(cDC2[["RNA"]], f = cDC2$sample) 
rm(seuratObj) 

DefaultAssay(cDC2) <- "RNA" 
cDC2 <- NormalizeData(cDC2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
cDC2 <- FindNeighbors(cDC2, dims = 1:50, reduction = "pca") 
cDC2 <- FindClusters(cDC2, resolution = 1, cluster.name = "unintegrated_clusters") 
cDC2 <- RunUMAP(cDC2, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated") 
#DimPlot(cDC2, reduction = "umap.unintegrated", group.by = c("treatment_status", "unintegrated_clusters")) 
#DimPlot(cDC2, group.by = "sample", reduction = "umap.unintegrated") 
cDC2 <- IntegrateLayers(object = cDC2, method = HarmonyIntegration,  
                        orig.reduction = "pca", new.reduction = "integrated.harmony", 
                        verbose = TRUE) 
cDC2[["RNA"]] <- JoinLayers(cDC2[["RNA"]]) 

cDC2 <- FindNeighbors(cDC2, dims = 1:50, reduction = "integrated.harmony")
cDC2 <- RunUMAP(cDC2, dims = 1:50, reduction = "integrated.harmony") #umap as the new reduction
cDC2 <- FindClusters(object = cDC2, resolution = c(0.01, 0.05, 0.08, seq(.1,2.1,.1)))

# Clustree 
#library(clustree)
#clustree_plot <- clustree(cDC2@meta.data, prefix = "RNA_snn_res."); clustree_plot # Manual saves

##### Clustering with SNN res of 0.1 
#cDC2 <- readRDS("D:/02_02_25_OS_LBH_workstream02/Myeloid_subsets_harmony_post_clustree_15_02_25.rds")
DimPlot(cDC2, group.by = "RNA_snn_res.0.1", reduction = "umap", label = T)
Idents(cDC2) <- "RNA_snn_res.0.1"
DimPlot(cDC2, label = T) # Manual saves

### Construct the file for reference ###
markers <- FindAllMarkers(cDC2,
                          logfc.threshold = 0.25, 
                          min.diff.pct = 0.2, 
                          only.pos = T)

markers$gene[which(markers$cluster == 0)][1:50]
cDC2 <- RenameIdents(cDC2, "0" = "cDC2_CD83hi") # 
cDC2 <- RenameIdents(cDC2, "1" = "cDC2_CD83lo") 
cDC2 <- RenameIdents(cDC2, "2" = "cDC2_CD83hi") 
cDC2 <- RenameIdents(cDC2, "3" = "cDC2_CD83hi") 
cDC2 <- RenameIdents(cDC2, "4" = "cDC2_CD83hi") 

## UMAP (Figure 6I) ----- 
DimPlot(cDC2)
FeaturePlot(cDC2, "CD83")
FeaturePlot(cDC2, c("FCER1A", "CD1C", "CLEC10A", "CD83"))

#if (!require("pak", quietly = TRUE)) {
#  install.packages("pak")
#}
#pak::pak("mengxu98/scop")
# Density and DimPlot as per CD83 
library(SCP)
library(scop)
library(viridis)
library(viridisLite)
library(Nebulosa)
cDC2_for_plot <- cDC2
cDC2_for_plot[["RNA"]] <- as(object = cDC2_for_plot[["RNA"]], Class = "Assay")
#p2 <- FeatureDimPlot(cDC2_for_plot, features = "CD83", reduction = "UMAP", theme_use = "theme_blank", add_density = F,
#               palette = "magma", bg_color = "#000004FF")
p2 <- plot_density(cDC2_for_plot, "CD83", reduction = "umap") + 
  #NoLegend() + 
  scale_color_viridis(option = "A") + 
  theme(axis.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank()) + 
  coord_fixed(1) 
p2

pdf("cDC2_CD83exp_251025.pdf", width = 3, height = 3)
print(p2)
dev.off()

cDC2_for_plot$celltype_Min <- Idents(cDC2_for_plot)
cDC2_for_plot$celltype_Min <- droplevels(cDC2_for_plot$celltype_Min)

p1 <- CellDimPlot(
  srt = cDC2_for_plot, group.by = "celltype_Min", stat.by = "group",
  #reduction = "UMAP", 
  theme_use = "theme_blank", 
  palcolor = c("#5B686D", "#C5C7C4"), 
  stat_palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
)
p1

pdf("cDC2_CD83_partition_251025.pdf", width = 5, height = 4)
print(p1)
dev.off()

## Cell proportions (Related to Figure 6J) -----
library(ggplot2)
 
Idents(cDC2_for_plot) <- "group"
cDC2_for_plot_pre <- subset(cDC2_for_plot, ident = c("A", "C"))
cDC2_for_plot_pre$sample <- droplevels(cDC2_for_plot_pre$sample)
cellnum <- table(cDC2_for_plot_pre$celltype_Min, cDC2_for_plot_pre$sample)
cell.prop <- as.data.frame(prop.table(cellnum))

colnames(cell.prop) <- c("Celltype", "Group", "Proportion")

p.bar <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#5B686D", "#C5C7C4")) +  
  #ggtitle("Cell Type Proportion by Group") +
  theme_bw() +
  NoLegend() + 
  coord_fixed(ratio = 5) + 
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  guides(fill = guide_legend(title = NULL))

print(p.bar)
 
Idents(cDC2_for_plot) <- "group"
cDC2_for_plot_pre <- subset(cDC2_for_plot, ident = c("B", "D"))
cDC2_for_plot_pre$sample <- droplevels(cDC2_for_plot_pre$sample)
cellnum <- table(cDC2_for_plot_pre$celltype_Min, cDC2_for_plot_pre$sample)
cell.prop <- as.data.frame(prop.table(cellnum))

colnames(cell.prop) <- c("Celltype", "Group", "Proportion")

p.bar2 <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#5B686D", "#C5C7C4")) +  
  #ggtitle("Cell Type Proportion by Group") +
  theme_bw() +
  NoLegend() + 
  coord_fixed(ratio = 5) + 
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  guides(fill = guide_legend(title = NULL))

print(p.bar2)

library(patchwork)
p3 <- p.bar + p.bar2 + plot_layout(ncol = 1) 

p3
pdf("cDC2_CD83_proportion_251025.pdf", width = 5, height = 4)
print(p3)
dev.off()

## Expression Profile of selection of Co-stimulatory and MHC II molecules (Figure 6L)  -----
seuratObj <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))
levels(seuratObj)
#DC <- subset(seuratObj, idents = c("cDC2_CD1C+", "cDC1_CLEC9A+", "Mature_DC"))
DC <- subset(seuratObj, idents = c("cDC2_CD1C+", "cDC1_CLEC9A+", "Mature_DC", "Macro_APOE+", "Macro_CCL4+", 
                                   "Mono_FCN1+", "Macro_ISGs+", "Macro_Ki67+",  "Macro_MTs+"
))
levels(DC)
DC <- RenameIdents(DC, "Macro_APOE+" = "MonoMac", "Macro_CCL4+" = "MonoMac", "Mono_FCN1+" = "MonoMac", 
                   "Macro_ISGs+" = "MonoMac", "Macro_Ki67+" = "MonoMac", "Macro_MTs+" = "MonoMac"
)
Idents(DC) <- Idents(cDC2)

DC[["RNA"]]$scale.data <- NULL
DC[["RNA"]]$data <- NULL
#cDC2[["RNA"]] <- split(cDC2[["RNA"]], f = cDC2$sample) 
#cDC2[["RNA"]]  
#rm(seuratObj)

DC <- NormalizeData(DC) %>% FindVariableFeatures() %>% ScaleData()
#saveRDS(DC, "Myeloid_w_DC_idents_active_251211.rds")

Features <- c("CD80", "CD86", "CD40", "TNFSF4", "TNFSF9", "ICOSLG", 
              "CCR7", "CD74", 
              rownames(DC)[grepl("HLA-D", rownames(DC))])
Features <- Features[!Features %in% "HLA-DQB1-AS1"]

DC$Cell_ident_for_MHC_plot <- Idents(DC) 
DC$Cell_ident_for_MHC_plot <- factor(DC$Cell_ident_for_MHC_plot, 
                                     levels = rev(c("MonoMac", "Mature_DC", "cDC1_CLEC9A+", 
                                                    "cDC2_CD83lo", "cDC2_CD83hi")))
Idents(DC) <- "Cell_ident_for_MHC_plot"
p <- DotPlot(DC, Features) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1 <- ggplot(p$data, aes(x = features.plot, y = id)) +   
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +   
  #facet_grid(facets = ~feature.groups,  switch = "x", scales = "free_x", space = "free_x") +    
  scale_radius(breaks = c(25, 50, 75, 100), range = c(0,6)) +   
  theme_classic() + 
  scale_color_gradient2(low = "#50859f", mid = "white", high = "#d66692") +  
  theme(axis.text.x = element_text(angle = 90,face = 1, size = 12, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, face = 1, color = "black"), 
        panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid"),
        legend.text = element_text(size = 8, face = 1),         
        legend.title = element_text(size = 10, face = 1),        
        legend.position = 'top', strip.placement = "outside", strip.text.x = element_blank(),        
        axis.title = element_blank()) + guides(colour = guide_colourbar(title.vjust = 0.9, title.hjust = 0)) + 
  labs(size = "Percent Expressed", color = "Average Expression") 
p1

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_", "Costim_and_MHCII.pdf"), 
    width = 9, # The width of the plot in inches
    height = 3) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()

rm(list = ls())
# disabled, run from repo root: setwd(here("data/NicheNet_Mye_and_T_interactions_run_04_03_25_MAC"))

## NicheNET (Related to Figure 6M) ----- 
#devtools::install_github("saeyslab/nichenetr")
library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)

#install.packages("tidyverse")
##### Load scRNA object #####
#seuratObj <- readRDS(here("data/NicheNet_Mye_and_T_interactions_run_04_03_25_MAC/Combined_TIL_and_Myeloid_file_w_fine_cell_type_v1_02_03_25.rds"))
seuratObj <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))
seuratObj <- UpdateSeuratObject(seuratObj)
Idents(seuratObj) <- "sample"
seuratObj <- RenameIdents(seuratObj, "H_BF_P1" = "Pre", "H_BF_P2" = "Pre", "H_BF_P3" = "Pre", "H_BF_P4" = "Pre", "H_BF_P5" = "Pre", "H_BF_P6" = "Pre", "H_BF_P7" = "Pre", "H_BF_P8" = "Pre", 
                          "H_AF_P1" = "Post", "H_AF_P2" = "Post", "H_AF_P3" = "Post", "H_AF_P4" = "Post", "H_AF_P5" = "Post", "H_AF_P6" = "Post", "H_AF_P7" = "Post", "H_AF_P8" = "Post", 
                          "L_BF_P1" = "Pre", "L_BF_P2" = "Pre", "L_BF_P3" = "Pre", "L_BF_P4" = "Pre", "L_BF_P5" = "Pre", "L_BF_P6" = "Pre", "L_BF_P7" = "Pre", 
                          "L_AF_P1" = "Post", "L_AF_P2" = "Post", "L_AF_P3" = "Post", "L_AF_P4" = "Post", "L_AF_P5" = "Post", "L_AF_P6" = "Post", "L_AF_P7" = "Post"
)

seuratObj$treatment_status <- Idents(seuratObj)
Idents(seuratObj) <- "Cell_type_fine_harmony"

seuratObj@meta.data$Cell_type_fine_harmony %>% table() 
lr_network <- readRDS(here("data/NicheNet_Mye_and_T_interactions_run_04_03_25_MAC/lr_network.rds"))
ligand_target_matrix <- readRDS(here("data/NicheNet_Mye_and_T_interactions_run_04_03_25_MAC/ligand_target_matrix.rds"))
weighted_networks <- readRDS(here("data/NicheNet_Mye_and_T_interactions_run_04_03_25_MAC/weighted_networks.rds"))

#organism <- "human"
#if(organism == "human"){
#  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
#  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
#} else if(organism == "mouse"){
#  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
#}

lr_network <- lr_network %>% distinct(from, to)
levels(seuratObj)

# Subsetting to post-treatment and immune cell subsets only 
to_keep <- levels(seuratObj)[!levels(seuratObj) 
                             %in% c("Proliferating", "PC_ISGs","Fibro-like","matPC","SMC_CXCL12",
                                    "SMC_MYH11", "artEC", "capEC", "venEC", "LEC", "EndoMT-I", 
                                    "EndoMT-II", "Cycling", "Chondrocytes", "Erythrocytes", "Fibroblasts", 
                                    "HSCs", "MSCs", "MastCells", "Osteoblasts", "Osteoclasts")]
seuratObj <- subset(seuratObj, idents = to_keep)
levels(seuratObj)
Idents(seuratObj) <- "group"
seuratObj <- subset(seuratObj, idents = c("B", "D"))

## Attaching DC labels ()
Idents(seuratObj) <- "Cell_type_fine_harmony"
DC <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/DC_251025/Myeloid_w_DC_idents_active_251211.rds"))
to_keep <- levels(DC)[!levels(DC) %in% c("MonoMac")]
DC <- subset(DC, idents = to_keep)
Idents(seuratObj) <- Idents(DC)
levels(seuratObj)
rm(DC)

seuratObj$Cell_type_fine_harmony_251212 <- Idents(seuratObj)

#for (i in unique(levels(seuratObj$Cell_type_fine_harmony_251212))) 
#for (i in "cDC2_CD83hi")
for (i in "cDC2_CD1C+") {
  ##### Define a set of potential ligands for both the sender-agnostic and sender-focused approach #####
  receiver = i
  expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.05)
  all_receptors <- unique(lr_network$to)  
  expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
  potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  # Use lapply to get the expressed genes of every sender cell type separately here
  seuratObj_idents <- as.character(unique(Idents(seuratObj)))
  sender_celltypes <- seuratObj_idents[which(seuratObj_idents != i)] 
  
  list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.05)
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  
  potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 
  
  print(paste0("Number of sender expressed genes = ", length(expressed_genes_sender)))
  print(paste0("Number of potential_ligands = ", length(potential_ligands)))
  print(paste0("Number of potential_ligands_focused = ", length(potential_ligands_focused)))
  
  ##### Define the gene set of interest ######
  ## Change as necessary 
  #condition_oi <-  "Post"
  #condition_reference <- "Pre"
  
  condition_oi <-  "B"
  condition_reference <- "D"
  
  seurat_obj_receiver <- subset(seuratObj, idents = receiver)
  #Idents(seurat_obj_receiver) <- "sample"
  Idents(seurat_obj_receiver) <- "group"
  #seurat_obj_receiver <- RenameIdents(seurat_obj_receiver, "H_BF_P1" = "Pre", "H_BF_P2" = "Pre", "H_BF_P3" = "Pre", "H_BF_P4" = "Pre", "H_BF_P5" = "Pre", "H_BF_P6" = "Pre", "H_BF_P7" = "Pre", "H_BF_P8" = "Pre", 
  #                                    "H_AF_P1" = "Post", "H_AF_P2" = "Post", "H_AF_P3" = "Post", "H_AF_P4" = "Post", "H_AF_P5" = "Post", "H_AF_P6" = "Post", "H_AF_P7" = "Post", "H_AF_P8" = "Post", 
  #                                    "L_BF_P1" = "Pre", "L_BF_P2" = "Pre", "L_BF_P3" = "Pre", "L_BF_P4" = "Pre", "L_BF_P5" = "Pre", "L_BF_P6" = "Pre", "L_BF_P7" = "Pre", 
  #                                    "L_AF_P1" = "Post", "L_AF_P2" = "Post", "L_AF_P3" = "Post", "L_AF_P4" = "Post", "L_AF_P5" = "Post", "L_AF_P6" = "Post", "L_AF_P7" = "Post"
  #)
  
  #seurat_obj_receiver$treatment_status <- Idents(seurat_obj_receiver)
  seurat_obj_receiver$group <- Idents(seurat_obj_receiver)
  DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                    ident.1 = condition_oi, ident.2 = condition_reference,
                                    #group.by = "treatment_status",
                                    group.by = "group",
                                    min.pct = 0.05) %>% rownames_to_column("gene")
  
  geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  
  ##### Define the background genes #####
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  print(paste0("Number of background expressed genes = ", length(background_expressed_genes)))
  print(paste0("Number of genes of interest = ", length(geneset_oi)))
  
  ##### Perform NicheNet ligand activity analysis #####
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  
  ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
  
  p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
    geom_histogram(color="black", fill="darkorange")  + 
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
               color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  
  p_hist_lig_activity
  
  best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
  
  vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
    column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)
  
  (make_heatmap_ggplot(vis_ligand_aupr,
                       "Prioritized ligands", "Ligand activity", 
                       legend_title = "AUPR", color = "darkorange") + 
      theme(axis.text.x.top = element_blank()))  
  
  ##### Infer target genes and receptors of top-ranked ligands #####
  ### Active target gene inference ###
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n = 100) %>%
    bind_rows() %>% drop_na()
  
  print(paste0("Number of active ligand target links df = ", nrow(active_ligand_target_links_df)))
  
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0) 
  
  nrow(active_ligand_target_links)
  ## [1] 86
  head(active_ligand_target_links)
  ##        Ifna13 Ifna2 Ifna6 Ifna15 Ifna7 Ifna5 Ifnab Ifna9 Ifna11 Ifna12 Ifna16 Ifna4 Ifna14 Ptprc        Tnf      Il36g       Il10      Il21        Osm
  ## Irf1        0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.27692301 0.07400782 0.07722567 0.1342983 0.16962803
  ## Ddx60       0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.11281871 0.00000000 0.05478472 0.0000000 0.08116101
  ## Parp14      0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.07003101 0.06895448 0.00000000 0.0000000 0.08011593
  ## Ddx58       0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.24433255 0.06891134 0.00000000 0.0000000 0.08862524
  ## Parp12      0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.19298997 0.06687691 0.05621734 0.0000000 0.07252823
  ## Tap1        0     0     0      0     0     0     0     0      0      0      0     0      0     0 0.25076038 0.07099514 0.00000000 0.0280451 0.15842935
  ##             Il27     Ifna1      Ifnb1      Ifng       Ifnk       Ifne      Lrtm2      Ifnl3       Ebi3      Ifnl2
  ## Irf1   0.3393635 0.2493108 0.25825704 0.2864755 0.04430047 0.04063913 0.03483987 0.10021371 0.11317944 0.07408235
  ## Ddx60  0.1596771 0.1218225 0.13569911 0.1171453 0.02629924 0.02771157 0.00000000 0.08217529 0.04882999 0.03746171
  ## Parp14 0.1563348 0.1269487 0.07891102 0.1142710 0.00000000 0.00000000 0.00000000 0.07074981 0.04568410 0.03135479
  ## Ddx58  0.2265024 0.2467722 0.21469112 0.2480807 0.00000000 0.00000000 0.00000000 0.07125327 0.04024212 0.02705719
  ## Parp12 0.1580405 0.1844883 0.14626411 0.1851128 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
  ## Tap1   0.1949126 0.1937607 0.16257249 0.2563000 0.00000000 0.02717588 0.00000000 0.08010402 0.04710621 0.03769675
  
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
  
  vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
  
  make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                      color = "purple", legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  
  ### Receptors of top-ranked ligands ###
  ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
    best_upstream_ligands, expressed_receptors,
    lr_network, weighted_networks$lr_sig) 
  
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
    ligand_receptor_links_df,
    best_upstream_ligands,
    order_hclust = "both") 
  
  (make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                       y_name = "Ligands", x_name = "Receptors",  
                       color = "mediumvioletred", legend_title = "Prior interaction potential"))
  
  ##### Sender-focused approach ####
  ligand_activities_all <- ligand_activities 
  best_upstream_ligands_all <- best_upstream_ligands
  
  ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
  best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
    pull(test_ligand) %>% unique()
  
  ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
    column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
  vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 
  
  p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                       "Prioritized ligands", "Ligand activity", 
                                       legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()) +
    theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid"))
  
  p_ligand_aupr
  
  # Target gene plot
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n = 100) %>%
    bind_rows() %>% drop_na()
  
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0.25) 
  
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
  
  vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
  
  p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                         color = "purple", legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  
  p_ligand_target
  
  # Receptor plot
  ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
    best_upstream_ligands, expressed_receptors,
    lr_network, weighted_networks$lr_sig) 
  
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
    ligand_receptor_links_df,
    best_upstream_ligands,
    order_hclust = "both") 
  
  p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                           y_name = "Ligands", x_name = "Receptors",  
                                           color = "mediumvioletred", legend_title = "Prior interaction potential")
  
  p_ligand_receptor
  
  best_upstream_ligands_all %in% rownames(seuratObj) %>% table()
  
  # Dotplot of sender-focused approach
  p_dotplot <- DotPlot(subset(seuratObj, Cell_type_fine_harmony %in% sender_celltypes), # Change to the celltype metadata name
                       features = rev(best_upstream_ligands), cols = "RdYlBu") + 
    coord_flip() +
    scale_y_discrete(position = "right") +
    theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0))
  
  p_dotplot
  
  celltype_order <- levels(Idents(seuratObj)) 
  # Use this if cell type labels are the identities of your Seurat object
  # if not: indicate the celltype_col properly
  
  DE_table_top_ligands <- lapply(
    celltype_order[celltype_order %in% sender_celltypes],
    get_lfc_celltype, 
    seurat_obj = seuratObj,
    #condition_colname = "treatment_status",
    condition_colname = "group",
    condition_oi = condition_oi,
    condition_reference = condition_reference,
    celltype_col = "Cell_type_fine_harmony",
    #celltype_col = "Cell_type_fine_harmony_251212",
    min.pct = 0, logfc.threshold = 0,
    features = best_upstream_ligands 
  ) 
  
  DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
    column_to_rownames("gene") 
  
  vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])
  
  p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                          "Prioritized ligands", "LFC in Sender",
                                          low_color = "midnightblue", mid_color = "white",
                                          mid = median(vis_ligand_lfc), high_color = "red",
                                          legend_title = "LFC") +
    theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0))
  
  
  p_lfc
  
  # Compare rankings 
  (make_line_plot(ligand_activities = ligand_activities_all,
                  potential_ligands = potential_ligands_focused) +
      theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))
  
  # Combining plots
  figures_without_legend <- cowplot::plot_grid(
    p_ligand_aupr + theme(legend.position = "none"),
    p_dotplot + theme(legend.position = "none",
                      axis.ticks = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_text(size = 12),
                      axis.text.y = element_text(size = 9),
                      axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
      ylab("Expression in Sender"),
    p_lfc + theme(legend.position = "none",
                  axis.title.y = element_blank(), 
                  panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid")),
    p_ligand_target + theme(legend.position = "none",
                            axis.title.y = element_blank(), 
                            panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid")),
    align = "hv",
    nrow = 1,
    rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))
  
  legends <- cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
    nrow = 1,
    align = "h", rel_widths = c(1.5, 1, 1, 1))
  
  combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
  combined_plot
  
  time <- gsub(" ", "_", Sys.time())
  time <- gsub("-", "_", time)
  time <- gsub(":", "_", time)
  pdf(file = paste0(time, i, "_as_receiver_", "NicheNet_w_all_T_Mye_subclusters.pdf"), 
      width = 40, # The width of the plot in inches
      height = 10) # The height of the plot in inches
  print(combined_plot, newpage = FALSE)
  dev.off()
}

# Adapted plot for cDC2 ----
pp1 <- p_ligand_aupr +
  coord_fixed() + 
  scale_fill_gradientn(colors = c('#e9f4f6', "#8f9fab", "#1e2235"))

pp1

pdf("cDC2_as_receiver_in_R_vs_NR_Post.pdf", width = 3, height = 7)
pp1
dev.off()

pp2 <- p_dotplot + 
  coord_fixed(ratio = 0.6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+ 
  scale_color_gradientn(colors = c('#50859f', "#EAEAEA", "#d66692"))

pp2
pdf("EXP_in_immune_celltypes_cDC2_as_receiver_in_R_vs_NR_Post.pdf", width = 10, height = 10)
pp2
dev.off()
