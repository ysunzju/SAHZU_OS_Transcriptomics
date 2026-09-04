library(here) # project-root-relative paths; run scripts from repository root
##################### TNK pipeline ##################### 
# Set up----- 
rm(list = ls())
#setwd()
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/T_and_NK"))
library(Seurat)
library(tidyverse)

scrna <- readRDS(here("data/28_04_25_OS_all_cell_plots/TNK_reident_ordered_01_05_25.rds"))
scrna[["RNA"]]$scale.data <- NULL
scrna[["RNA"]]$data <- NULL
scrna[["RNA"]] <- split(scrna[["RNA"]], f = scrna$sample) 

DefaultAssay(scrna) <- "RNA"
scrna <- NormalizeData(scrna) %>% FindVariableFeatures() %>% ScaleData()

scrna <- RunPCA(scrna)
scrna <- FindNeighbors(scrna, dims = 1:50, reduction = "pca")
scrna <- FindClusters(scrna, resolution = 1, cluster.name = "unintegrated_clusters")
scrna <- RunUMAP(scrna, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(scrna, reduction = "umap.unintegrated", group.by = c("treatment_status", "unintegrated_clusters"))
DimPlot(scrna, group.by = "sample", reduction = "umap.unintegrated", split.by = "treatment_status")
scrna <- IntegrateLayers(object = scrna, method = HarmonyIntegration, 
                         orig.reduction = "pca", new.reduction = "integrated.harmony",
                         verbose = TRUE)
scrna[["RNA"]] <- JoinLayers(scrna[["RNA"]])

scrna <- FindNeighbors(scrna, dims = 1:50, reduction = "integrated.harmony")
scrna <- RunUMAP(scrna, dims = 1:50, reduction = "integrated.harmony") #umap as the new reduction
scrna <- FindClusters(object = scrna, resolution = c(0.01, 0.05, 0.08, seq(.1,2.1,.1)))

##### Clustering with SNN res of 0.6 #####
Idents(scrna) <- "RNA_snn_res.0.6"
DimPlot(scrna, label = T) # Manual saves

markers <- FindAllMarkers(scrna,
                          logfc.threshold = 0.1, 
                          min.diff.pct = 0.2, 
                          only.pos = T)

### Construct the file for refererence ###
library(msigdbr)
library(clusterProfiler)
hs_df = msigdbr(species = "Homo sapiens") %>% as.data.frame()
hs_C8 = msigdbr(species = "Homo sapiens",
                category = "C8",
                subcategory = NULL) %>% as.data.frame() %>% 
  dplyr::select(gs_cat, gs_subcat, gs_name, gene_symbol)

hs_C8 = hs_C8 %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  dplyr::rename("term" = "gs_name", "gene" = "gene_symbol")
rm(hs_df)

library(org.Hs.eg.db)
organism = 'hsa'    #  人类'hsa' 小鼠'mmu'   
OrgDb = 'org.Hs.eg.db'

##### 0 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "0", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 0)][1:50]

scrna <- RenameIdents(scrna, "0" = "HelperT_NR4A1") # IL7R, MYADM, NR4A1, NR4A2/3, PTGER4

##### 1 #####
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "1", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 1)][1:50]

scrna <- RenameIdents(scrna, "1" = "HelperT_IL7R") # IL7R, KLRB1, LTB

##### 2 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "2", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 2)][1:50]

scrna <- RenameIdents(scrna, "2" = "CD8T_GZMK") # 

##### 3 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "3", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 3)][1:50]

scrna <- RenameIdents(scrna, "3" = "CD8T_IFNG") # 

##### 4 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "4", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 4)][1:50]

scrna <- RenameIdents(scrna, "4" = "Trm_ZNF683") # ZNF683, LEF1, CD7

##### 5 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "5", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 5)][1:50]
scrna <- RenameIdents(scrna, "5" = "Likely stromal doublet") # SPARC, COL1A2, LUM 

##### 6 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "T_RBhi_for_del", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 6)][1:50]

scrna <- RenameIdents(scrna, "6" = "Tprolif") # STMN1, MKI67, TOP2A 

##### 7 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "7", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 7)][1:50]

scrna <- RenameIdents(scrna, "7" = "NK_XCL1") # XCL1, XCL2, KLRD1

##### 8 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "8", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 8)][1:50]

scrna <- RenameIdents(scrna, "8" = "MAIT") # KLRB1, SLC4A10, NCR3

##### 9 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "9", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 9)][1:50]

scrna <- RenameIdents(scrna, "9" = "Dying cells") # nCount_low, Mt hi

##### 10 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "10", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 10)][1:50]

scrna <- RenameIdents(scrna, "10" = "T_ISGs") 

##### 11 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "11", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 11)][1:50]

scrna <- RenameIdents(scrna, "11" = "Myeloid_doublets") # Defs doublets given myeloid markers +++, C1QB, CD14, CD1C, SELENOP, APOE

##### 12 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "12", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 12)][1:50]

scrna <- RenameIdents(scrna, "12" = "Treg") # 

##### 13 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "13", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 13)][1:50]

scrna <- RenameIdents(scrna, "13" = "NK_FCGR3A") # FCGR3A (CD16), GNLY, GZMB

##### 14 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "14", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 14)][1:50]

scrna <- RenameIdents(scrna, "14" = "Tn") 

##### 15 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "15", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 15)][1:50]

scrna <- RenameIdents(scrna, "15" = "gdT") #" TRDV2"      "TRGV9"      "TRDC"       "TRGC1" 

##### 16 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "16", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 16)][1:50]

scrna <- RenameIdents(scrna, "16" = "Tex_CXCL13") #" TRDV2"      "TRGV9"      "TRDC"       "TRGC1" 

##### 17 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "17", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 17)][1:50]

scrna <- RenameIdents(scrna, "17" = "iNK") # KIT 

# Renaming -----
unique(Idents(scrna))
#[1] CD8T_GZMK              Myeloid_doublets       Dying cells            Likely stromal doublet CD8T_IFNG              HelperT_IL7R          
#[7] Treg                   HelperT_NR4A1          Tprolif                Tex_CXCL13             gdT                    NK_FCGR3A             
#[13] T_ISGs                 MAIT                   NK_XCL1                Trm_ZNF683             Tn                     iNK    
scrna <- subset(scrna, idents = c("Likely stromal doublet", "Myeloid_doublets"), invert = T)
unique(Idents(scrna))
scrna$Cell_type_fine_harmony <- Idents(scrna)

scrna[["RNA"]]$scale.data <- NULL
scrna[["RNA"]]$data <- NULL
scrna[["RNA"]] <- split(scrna[["RNA"]], f = scrna$sample) 
scrna[["RNA"]] ##查看一下
scrna <- NormalizeData(scrna) %>% FindVariableFeatures() %>% ScaleData()
scrna <- RunPCA(scrna)
scrna <- FindNeighbors(scrna, dims = 1:50, reduction = "pca")
scrna <- FindClusters(scrna, resolution = 1, cluster.name = "unintegrated_clusters")
scrna <- RunUMAP(scrna, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(scrna, reduction = "umap.unintegrated", group.by = c("treatment_status", "unintegrated_clusters"))
DimPlot(scrna, group.by = "sample", reduction = "umap.unintegrated", split.by = "treatment_status")
scrna <- IntegrateLayers(object = scrna, method = HarmonyIntegration, 
                         orig.reduction = "pca", new.reduction = "integrated.harmony",
                         verbose = TRUE)
scrna[["RNA"]] <- JoinLayers(scrna[["RNA"]])

scrna <- FindNeighbors(scrna, dims = 1:50, reduction = "integrated.harmony")
scrna <- RunUMAP(scrna, dims = 1:50, reduction = "integrated.harmony") #umap as the new reduction
DimPlot(scrna, group.by = "Cell_type_fine_harmony")
Idents(scrna) <- "Cell_type_fine_harmony"

updated_markers <- FindAllMarkers(scrna,
                                  logfc.threshold = 0.1, 
                                  min.diff.pct = 0.2, 
                                  only.pos = T)

markers_v1 <- unique(updated_markers$gene[which(updated_markers$cluster %in% c("CD8T_IFNG", "CD8T_GZMK", "HelperT_NR4A1", "HelperT_IL7R", "MAIT"))])
FeaturePlot(scrna, markers_v1, raster = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
spe <- FindMarkers(scrna, ident.1 = "Dying cells")

#write.csv(spe, "Dying_degs.csv", row.names = FALSE, quote = TRUE)
#write.table(spe, row.names = T, col.names = T, "Dying_degs.csv")
saveRDS(scrna, "Temp_TNK_post_del_reint_pending_TSTR_cluster_16_06_25.rds")

## As discussed, Dying cells to be named MT
scrna <- RenameIdents(scrna, "Dying cells" = "T_Mtgenes_hi")
Idents(scrna) <- droplevels(Idents(scrna))
scrna$Cell_type_fine_harmony <- Idents(scrna)
saveRDS(scrna, "T_NK_final_24_06_25.rds")

DimPlot(scrna, label = T)
updated_markers <- FindAllMarkers(scrna,
                                  logfc.threshold = 0.1, 
                                  min.diff.pct = 0.2, 
                                  only.pos = T)

### TNK count (Figure 5A) ----- 
dev.off()
Idents(scrna) <- "Cell_type_fine_harmony"
reduc <- data.frame(Seurat::Embeddings(scrna, reduction = "umap"))
meta <- scrna@meta.data
pc12 <- cbind(reduc, meta)
pc12$idents <- Idents(scrna)

cell_num <- pc12 |>
  dplyr::group_by(idents,Cell_type_fine_harmony) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::arrange(n)

cols <- mycol[cell_num$Cell_type_fine_harmony]
pushViewport(viewport(x = unit(0.61, "npc"), y = unit(0.5, "npc"),
                      width = unit(0.2, "npc"),
                      height = unit(0.7, "npc"),
                      just = "left",
                      # yscale = extendrange(c(0,nrow(cell_num)),f = 0.05),
                      # xscale = extendrange(c(0,max(cell_num$n)),f = 0.05),
                      yscale = c(0,nrow(cell_num) + 0.5),
                      xscale = c(0,max(cell_num$n) + 0.1*max(cell_num$n))))
grid.xaxis()
# grid.yaxis(main = F)
# grid.rect()
grid.rect(x = rep(0,nrow(cell_num)),y = 1:nrow(cell_num),
          width = cell_num$n,height = unit(0.08,"npc"),
          just = "left",
          gp = gpar(fill = cols,col = NA),
          default.units = "native") 
grid.rect(gp = gpar(fill = "transparent"))
grid.text(label = "Number of cells",x = 0.5,y = unit(-2.5,"lines"))
popViewport()

pushViewport(viewport(x = unit(0.81, "npc"), y = unit(0.5, "npc"),
                      width = unit(0.2, "npc"),
                      height = unit(0.7, "npc"),
                      just = "left",
                      yscale = c(0,nrow(cell_num) + 0.5)))
# grid.rect()
grid.points(x = rep(0.1,nrow(cell_num)),y = 1:nrow(cell_num),pch = 19,
            gp = gpar(col = cols),
            size = unit(1.5, "char"))
cell_num$numbering <- c(1:nrow(cell_num))
grid.text(label = cell_num$numbering,x = 0.1,y = 1:nrow(cell_num),
          default.units = "native")
grid.text(label = cell_num$idents,x = 0.2,y = 1:nrow(cell_num),
          just = "left",
          default.units = "native")

popViewport()

## Cell counts (Fig 5B) ----- 
TNK <-readRDS("T_NK_final_24_06_25.rds")
Idents(TNK) <- "group"
TNK <- RenameIdents(TNK, "A" = "R", 
                    "B" = "R", 
                    "C" = "NR", 
                    "D" = "NR")
TNK$response <- Idents(TNK)

TNK$Cell_type_fine_harmony <- droplevels(TNK$Cell_type_fine_harmony)
Cellratio <- table(TNK$Cell_type_fine_harmony, TNK$treatment_status, TNK$response)
Cellratio <- data.frame(Cellratio)

colnames(Cellratio) <- c("CellType", "Treatment", "Response", "Count")

mycol <- c(
  "#A7B3C6", "#4AAE79", "#D1A054", "#AC3E4E", "#8ED9E5", "#E89E2A", "#CED350", "#D94C43", 
  "#9A7FC1", "#A66C5B", "#60C8AE", "#B48A78", "#F6E27F", "#708CC2", "#F4A6B8", "#3A6CA8"
)

p1 <- ggplot(Cellratio, aes(x=Treatment, y=Count, fill=CellType, alluvium = CellType, stratum = CellType)) +
  #geom_bar(stat="identity", width = 0.5) +
  #geom_col(position = "fill") +
  facet_wrap(~ Response, ncol=3) +
  theme(strip.background = element_rect(fill="grey90")) +
  geom_flow(aes(fill = CellType), alpha = 0.5, width = 0.5) +
  scale_x_discrete(expand = expansion(add = 0)) + 
  geom_stratum(aes(fill = CellType), width = 0.5) +
  scale_fill_manual(values = mycol) +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.line.y = element_line())
p1

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "TNK_counts_per_group.pdf"), 
    width = 6,  
    height = 8)  
print(p1, newpage = FALSE)
dev.off()

## Pi (Figure 5D) -----
library(reshape2)
Idents(scrna) <- "reportname"
scrna <- RenameIdents(scrna,
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
scrna$TNR_val <- Idents(scrna)

Idents(scrna) <- "treatment_status"
scrna_pre <- subset(scrna, idents = "Pre")
Cellratio <- prop.table(table(scrna_pre$Cell_type_fine_harmony, scrna_pre$sample), margin = 2)  
Cellratio <- data.frame(Cellratio)
cellper <- reshape2::dcast(Cellratio, Var2~Var1, value.var="Freq") 
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]
cellper <- na.omit(cellper)

meta <- scrna_pre@meta.data
colnames(meta)
meta <- meta[,c(4, 60)] # TNR and sample columns 
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
#cellper$group <- NA
meta <- unique(meta) 
cellper <- merge(cellper, meta, by = "sample")
colnames(cellper)
#colnames(cellper)[18] <- "group"
cellper <- as.data.frame(cellper)
#cellper$patient <- str_replace(cellper$sample, "_BF_", "")

# Calc 
library(ggplot2)
library(dplyr)
cellper$TNR_val <- as.numeric(as.character(cellper$TNR_val))

# Create the matrix for plot
matrix_for_Plots <- matrix(data = rep(NA, 3*length(unique(scrna$Cell_type_fine_harmony))), 
                           ncol = 3, nrow = length(unique(scrna$Cell_type_fine_harmony)), dimnames = list(c(as.vector(unique(scrna$Cell_type_fine_harmony))), c("Cell_type", "Pi", "pval"))) %>% as.data.frame()
matrix_for_Plots$Cell_type <- rownames(matrix_for_Plots)

for (i in as.vector(unique(scrna$Cell_type_fine_harmony))) {
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
#mycol <- c("#92699e", "#a5a9b0", "#696a6c", "#d4c2db", "#b95055", "#6d6fa0", "#c49abc", "#d8a0c0",
#"#63a3b8", "#9a70a8", "#4490c4", "#a44e89", "#927c9a", "#de8b36", "#9f8d89", "#408444")
zhongguose_palette_23 <- rev(c(
  "#A7B3C6", "#4AAE79", "#D1A054", "#AC3E4E", "#8ED9E5", "#E89E2A", "#CED350", "#D94C43", 
  "#9A7FC1", "#A66C5B", "#60C8AE", "#B48A78", "#F6E27F", "#708CC2", "#F4A6B8", "#3A6CA8"
))

mycol <- rev(zhongguose_palette_23)
names(mycol) <- levels(scrna$Cell_type_fine_harmony)
mycol <- mycol[levels(matrix_for_Plots$Cell_type)]

#bg_colors = c("white", "#F6F9E4")
#gradient_grob <- rasterGrob(colorRampPalette(bg_colors)(256), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
p <- ggplot(matrix_for_Plots, aes(x = Cell_type, y = Pi)) +  
  #annotation_custom(gradient_grob,xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_segment(aes(x = Cell_type, xend = Cell_type, y = 0, yend = Pi, color = Cell_type),                 
               linetype = "solid", size = 1, color = mycol) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour="gray80") +  
  geom_point(aes(color = Cell_type, size = log10pval), color = mycol
  ) + 
  scale_size_continuous(range=c(3,15)) +
  #geom_text(aes(label = Ti), size = 3) +  
  geom_text(aes(label = Cell_type, y = 0), 
            hjust = ifelse(matrix_for_Plots$Pi >= 0, 1, 0), 
            #vjust = -0.5, 
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

## Cell prop plots (Fig 5J,M,P) -----
scrna <- readRDS("T_NK_final_24_06_25.rds")
Idents(scrna) <- "group"
scrna <- RenameIdents(scrna, "A" = "H_Pre", 
                      "B" = "H_Post", 
                      "C" = "L_Pre", 
                      "D" = "L_Post")
scrna$group_anno <- Idents(scrna)
scrna$Cell_type_fine_harmony <- droplevels(scrna$Cell_type_fine_harmony)
Cellratio <- prop.table(table(scrna$Cell_type_fine_harmony, scrna$sample), margin = 2) 
Cellratio <- data.frame(Cellratio)

cellper <- reshape2::dcast(Cellratio, Var2~Var1, value.var="Freq") 
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]

#write.csv()
#to_save <- t(cellper)
write.csv(to_save, file='TNK_cells_prop_260806.csv', quote = F)

meta <- scrna@meta.data
colnames(meta)
meta <- meta[,c(53,4)] # Group and sample
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
#cellper$group <- NA
meta <- unique(meta)
cellper <- merge(cellper, meta, by = "sample")
colnames(cellper)
colnames(cellper)[18] <- "group" # Group anno to group 
cellper <- as.data.frame(cellper)

pplist =list()
scrna_groups = unique(levels(scrna$Cell_type_fine_harmony))

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
#install.packages("ggtext")
library(ggsignif)
library(stringr)

mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
#fill_color <- zzm60colors
#names(fill_color) <- unique(levels(scrna$cluster))

for(group_ in scrna_groups){
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
  #p_position2 = min(max(cellper_$percent)*1.2) * 0.8
  
  print(group_)
  print(cellper_$median)
  mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape =21, aes(fill = group), width = 0.25) + 
    geom_boxplot(alpha = 0.5, size = 1, outlier.shape = NA, aes(fill = group)) + 
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) + 
    #stat_summary(fun = mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size =10),axis.title = element_text(size = 10),legend.text = element_text(size =10),
          legend.title = element_text(size =10),plot.title = element_text(size =10,face ='plain'),legend.position ='none')+
    labs(title = group_, y = "Percentage", x = "Group") + 
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = c(0, 0)) + 
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
          #plot.title = element_textbox_simple(size = 10, color = "black", halign = 0.5,
          #                                    fill = fill_color[group_], width = 1.2, 
          #                                    padding = margin(3, 0, 3, 0),
          #                                    margin = margin(0, 0, 10, 0)),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
          axis.text.y = element_text(color = "black", size = 9),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "none") + 
    theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid"))
  
  labely = max(cellper_$percent)
  compare_means(percent ~group,  data = cellper_)
  my_comparisons <-list(c("H_Pre","H_Post") 
                        #c("H_Pre","L_Pre"), 
                        #c("L_Pre","L_Post")
  )
  pp1 = pp1 + geom_signif(comparisons = list(c("H_Pre","H_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position, 
                          paired = T, 
                          test = "t.test",
                          test.args = "two.sided", 
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  #my_comparisons1 <-list(c("L_Pre","L_Post"))
  pp1 = pp1 + geom_signif(comparisons = list(c("L_Pre","L_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))}, 
                          y_position = p_position, 
                          test = "t.test",
                          test.args = "two.sided", 
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("H_Pre","L_Pre")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position1, 
                          test = "t.test", 
                          test.args = "two.sided", 
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black") 
  #pp1 = pp1 + geom_signif(comparisons = list(c("H_Post","L_Post")), 
  #                        map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
  #                        y_position = p_position2, 
  #                        test = "t.test", 
  #                        test.args = "two.sided", 
  #                        textsize = 4, tip_length = 0,
  #                        parse = TRUE, 
  #                        color = "black")
  
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
                 ncol = 8#, nrow = 2
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Before_integration_TNK_cell_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 16,  
    height = 6)  
print(pps, newpage = FALSE)
dev.off()

pps <- plot_grid(pplist[[1]], 
                 pplist[[2]],
                 pplist[[3]],
                 pplist[[3]],
                 pplist[[3]],
                 ncol = 4#, nrow = 2
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Th_cell_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 6,  
    height = 4.5)  
print(pps, newpage = FALSE)
dev.off()

pps <- plot_grid(pplist[[6]], 
                 pplist[[7]],
                 pplist[[8]],
                 pplist[[8]],
                 pplist[[8]],
                 ncol = 4#, nrow = 2
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Tc_cell_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 6,  
    height = 4.5)  
print(pps, newpage = FALSE)
dev.off()

## VlnPlots of all TNK cells for selected genes (Figure S5B) -----
library(ggplot2)
library(reshape2)
dt <- scrna@meta.data[which(names(scrna@meta.data) %in% c("Cell_type_fine_harmony"))]
dt$Cell_names <- rownames(dt)
#mat <- GetAssayData(scrna, slot = "data")[genes, ] %>% t() %>% as.data.frame()
mat <- data.frame(scrna@meta.data[, c(54:56)])
mat$Cell_names <- rownames(mat)
dt <- merge(dt, mat, by = "Cell_names")
dt <- dt[, -1]

dt <- reshape2::melt(dt,
                     id.vars = c("Cell_type_fine_harmony"), # cols to keep 
                     measure.vars = names(scrna@meta.data)[54:56], # cols to combine 
                     variable.name = "signature", # the colname of the combined items
                     value.name = "expressions") # the colname of the values 
dt$Cell_type_fine_harmony <- factor(dt$Cell_type_fine_harmony, levels = rev(
  c("Tn_CCR7+", "Tmem_IL7R+", "Tmem_GPR183+FOSB+", "Treg", 
    "Trm_ZNF683+", "CD8_Teff_IFNG+", "CD8_Teff_GZMK+", "Tex_CXCL13+", "CD8_MAIT", 
    "gdT", "T_ISGs", "T_Mtgenes_hi", "T_Ki67+", "NK_KIT+", "NK_KLRC1+", "NK_FCGR3A+"
  )
)
)
dt$signature <- factor(dt$signature, levels = unique(dt$signature))

zhongguose_palette_23 <- rev(c(
  "#A7B3C6", "#4AAE79", "#D1A054", "#AC3E4E", "#8ED9E5", "#E89E2A", "#CED350", "#D94C43", 
  "#9A7FC1", "#A66C5B", "#60C8AE", "#B48A78", "#F6E27F", "#708CC2", "#F4A6B8", "#3A6CA8"
))

p3 <- ggplot(data = dt,
             aes(x = expressions, y = Cell_type_fine_harmony, fill = Cell_type_fine_harmony)) +
  geom_violin(scale = 'width',
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = 'black',
              size = 0.5, 
              alpha = 0.6) + 
  facet_grid(cols = vars(signature), scales = 'free_x') +
  scale_fill_manual(values = zhongguose_palette_23) + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    axis.text.x = element_text(size = 8), 
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size = 8), 
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8), 
    legend.position = "none"
  ) +
  labs(x = 'Log Normalized Expression') 
p3

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Vlnplots_for_TNK_sigs.pdf"), 
    width = 4,  
    height = 6)  
print(p3, newpage = FALSE)
dev.off() 

## Feature density plot (Figure S5C)----- 
library(Nebulosa)
markers_for_plot <- vector()
for (i in unique(updated_markers$cluster)) {
  m <- updated_markers$gene[which(updated_markers$cluster == i)][1:3]
  markers_for_plot <- union(markers_for_plot, m)
}
markers_for_plot <- unique(markers_for_plot)

p <- plot_density(scrna, markers_for_plot, reduction = "umap", combine = F, raster = T
                  #raster = F, max.cutoff = 0.5, 
                  #cols = c("lightgrey", "#A6559D"), 
                  #reduction = "umap.unintegrated" 
)

for (i in 1:length(p)){
  p[[i]] = p[[i]] + 
    scale_color_viridis_c(option = "magma") + 
    theme(#panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid"),
      plot.title = element_text(margin = margin(t = 0, b = 0)),
      axis.ticks = element_blank(),
      axis.line.x = element_blank(), 
      axis.line.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(), 
      plot.subtitle = element_blank(), 
      legend.position = "none", 
      #panel.background = element_rect(fill = "white"), 
      panel.grid = element_blank()
    ) + 
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank()) + 
    annotate("text", x = range(p[[i]]$data$umap_1)[2]*0.8, y = range(p[[i]]$data$umap_2)[2], size = 5,
             label = p[[i]][["labels"]][["title"]],
             fontface="italic",
             colour="black") +
    theme(plot.title = element_blank())
}

f = CombinePlots(plots = list(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], 
                              p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], 
                              p[[11]], p[[12]], p[[13]], p[[14]], p[[15]], 
                              p[[16]], p[[17]], p[[18]], p[[19]], p[[20]],
                              p[[21]], p[[22]], p[[23]], p[[24]], p[[25]], 
                              p[[26]], p[[27]], p[[28]], p[[29]], p[[30]], 
                              p[[31]], p[[32]], p[[33]], p[[34]], p[[35]],
                              p[[36]], p[[37]], p[[38]], p[[39]], p[[40]],
                              p[[41]], p[[42]], p[[43]]), ncol = 5)
f

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "FeaturePlot_TNK_cluster_cluster_marker_Genes_v3_20_60_25.pdf"), 
    width = 15,  
    height = 27)  
print(f, newpage = FALSE)
dev.off()

### Scoring of CD8, T, and NK cells ----- 
mat <- readxl::read_xlsx(here("data/Re-clustering_fine_cell_types_21_05_25/T_and_NK/NK_CD8_gdT_sigs.xlsx"))
genelist <- list(CD8 = c("CD8B", "CD8A", "CD3D", "LTB", "CD27", "COTL1", "LDHB", "TRAC"), NK = mat$NK, gdT = mat$gdT[!is.na(mat$gdT)])

scrna <- AddModuleScore(scrna, features = genelist, name = names(genelist))
#scrna@meta.data[, 53:56] <- NULL
names(scrna@meta.data)[54:56] <- names(genelist)
FeaturePlot(scrna, names(genelist), keep.scale = "all", ncol =1)
VlnPlot(scrna, names(genelist), pt.size = 0, ncol = 1)
#VlnPlot(scrna, "CD8A")

## Monocle 3 trajectory (Fig 5H,I,K,L,N,O)----- 
### Tmem ----- 
library(monocle3)
Idents(scrna) <- "Cell_type_fine_harmony"
unique(Idents(scrna)) # Selected subsets only 
scrna_subset <- subset(scrna, idents = c("Tn_CCR7+", "Tmem_IL7R+", "Tmem_GPR183+FOSB+"))

data <- GetAssayData(scrna_subset, assay = 'RNA', layer = 'counts')
cell_metadata <- scrna_subset@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(expression_dat = data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds) # Preprocessing data 
cds <- align_cds(cds, alignment_group = "sample") # Remove batch effects
cds <- reduce_dimension(cds) # Dimension reduction
cds <- cluster_cells(cds) # Cluster cells 
cds <- learn_graph(cds) # Predict trajectory 
cds <- learn_graph(cds, use_partition = T, 
                   close_loop = T) 

get_earliest_principal_node <- function(cds, time_bin="Tn_CCR7+"){
  cell_ids <- which(colData(cds)[, "Cell_type_fine_harmony"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

mycol <- c(
  "#A7B3C6", "#4AAE79", "#D1A054", "#AC3E4E", "#8ED9E5", "#E89E2A", "#CED350", "#D94C43", 
  "#9A7FC1", "#A66C5B", "#60C8AE", "#B48A78", "#F6E27F", "#708CC2", "#F4A6B8", "#3A6CA8"
)
#mycol <- celltype_colors
names(mycol) <- levels(scrna$Cell_type_fine_harmony)
mycol <- mycol[levels(scrna$Cell_type_fine_harmony)]

p1 <- plot_cells(cds,
                 color_cells_by = "Cell_type_fine_harmony",
                 label_cell_groups=F,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size = 1.5, 
                 group_label_size = 2,
                 trajectory_graph_color = "#EAEAEA", 
                 trajectory_graph_segment_size = 0.5, 
                 cell_size = 0.5) + 
  scale_color_manual(values = mycol[c("Tn_CCR7+", "Tmem_IL7R+", "Tmem_GPR183+FOSB+")]) + 
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.line.y = element_blank(), 
        panel.border = element_blank())
p1 <- p1 + coord_equal()

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_cell_type_plots_for_Th.pdf"), 
    width = 5,  
    height = 5)  
print(p1, newpage = FALSE)
dev.off()

p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5, 
                 trajectory_graph_color = "#CCCCCC", 
                 trajectory_graph_segment_size = 0.5)
p2 <- p1 + theme(legend.position = "none",
                 axis.ticks.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.title.x = element_blank(), 
                 axis.line.x = element_blank(), 
                 axis.line.y = element_blank(), 
                 panel.border = element_blank()) +
  scale_color_gradientn(colours = c(
    "#3A9ED9", "#4FA1DC", "#64A4DE", "#79A7E0",
    "#8EABE2", "#A3AFE4", "#B7AEE2", "#C4A4DD",
    "#D19BD8", "#DE92D3", "#EB89CE", "#F080C9"
  ))
p2 <- p2 + coord_equal()

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_Th.pdf"), 
    width = 5,  
    height = 5)  
print(p2, newpage = FALSE)
dev.off()

##### Graph test  ----- 
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.1))

##### Heatmap plot of key genes (w complext heatmap) ------ 
library(dplyr)
library(tidyr)
library(viridis)
library(ClusterGVis)
pseudotime <- pseudotime(cds) %>% as.data.frame()
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "pseudotime"

celltype <- cds@colData$Cell_type_fine_harmony %>% as.data.frame()
celltype$cell <- colnames(cds)
colnames(celltype)[1] <- "clusters"

Treatment <- cds@colData$group_anno %>% as.data.frame()
Treatment$cell <- colnames(cds)
colnames(Treatment)[1] <- "treatment"

merge <- merge(pseudotime, celltype, by = 'cell')
merge <- merge(merge, Treatment, by = 'cell')

#devtools::install_github("junjunlab/ClusterGVis", force = T)
exp <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

exp <- data.frame(t(exprs(cds)))
exp <- exp[, c("LEF1", "SELL", "TCF7", "IL7R", "KLRB1", "EGR1", "ID2", "CD69", "FOSB",
               "GPR183", #"CD27", 
               "CD28", 
               #"CD40LG", 
               "BIRC3"#, 
               #"BCL2L1"
)]

library(dplyr)
heatmap_data_scaled <- data.frame(scale(exp)) 

library(ComplexHeatmap)
library(circlize)
heatmap_data_scaled$cell <- rownames(heatmap_data_scaled)
merge_mat <- merge(merge, heatmap_data_scaled, by = "cell")
merge_mat <- merge_mat[order(merge_mat$pseudotime), ]
merge_mat_for_plot <- as.data.frame(t(merge_mat)[c(5:16), ]) # The Gene Exp Matrix
for (i in 1:ncol(merge_mat_for_plot)) {
  merge_mat_for_plot[, i] <- as.numeric(merge_mat_for_plot[, i])
}

library(circlize)
range(merge_mat$pseudotime)[2]
col_fun = colorRamp2(seq(0, range(merge_mat$pseudotime)[2], 
                         length.out = 12), c(
                           "#3A9ED9", "#4FA1DC", "#64A4DE", "#79A7E0",
                           "#8EABE2", "#A3AFE4", "#B7AEE2", "#C4A4DD",
                           "#D19BD8", "#DE92D3", "#EB89CE", "#F080C9"
                         ))

ht_list = HeatmapAnnotation(
  Pseudotime = anno_barplot(merge_mat$pseudotime, 
                            gp = gpar(fill = col_fun(merge_mat$pseudotime), col = col_fun(merge_mat$pseudotime)), 
                            height = unit(1.5, "cm")), 
  Celltype = merge_mat$clusters, 
  Treatment = merge_mat$treatment,
  col = list(Celltype = c("Tn_CCR7+" = "#A7B3C6", "Tmem_IL7R+" = "#4AAE79", "Tmem_GPR183+FOSB+" = "#D1A054"), 
             Treatment = c("H_Pre" = "#54426D", "H_Post" = "#D9A0B3", "L_Pre" = "#0f5688", "L_Post" = "#6B798E"))
)

cellheight = 0.5
rn = dim(as.matrix(merge_mat_for_plot))[1]
h=cellheight*rn
pdf(file = "Th_pseudotime_gene_HM.pdf", 
    width = 6,  
    height = 6)  
Heatmap(as.matrix(#merge_mat[c("ATF3", "BCL2A1", "CD83", "FOS", "IL1B", "JUN", "NR4A1", "NR4A2", "NR4A3"), ]
  merge_mat_for_plot),
  # Height of the blocks
  #width = unit(w, "cm"),
  height = unit(h, "cm"),
  #scale = F,
  #rect_gp = gpar(col = "white", lwd = 1.5),
  #border_g = gpar(col = ,lty = 1,lwd = 1.2),
  # Dent formatting 
  #column_dend_height = unit(1.5, "cm"), 
  #row_dend_width = unit(1.5, "cm"),
  #column_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
  #row_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
  # Gap numbers and style 
  #row_split = 6, column_split = c("1", "1", "2", "2"),
  #row_gap = unit(2, "mm"),
  #column_gap = unit(2, "mm"),
  # Text formats 
  row_title = NULL,column_title = NULL,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  show_heatmap_legend = T, 
  cluster_rows = F,
  cluster_columns = F,
  row_names_side = 'right',
  show_column_names = F,
  show_row_names = T,
  border = T, 
  top_annotation = ht_list, 
  col = colorRamp2(breaks = seq(range(merge_mat_for_plot)[1], 2.5, length.out = 10), 
                   colors = viridis(10, alpha = 1, begin = 0, end = 1, direction = -1, option = "F"))
)
dev.off()

### Tc ----- 
library(monocle3)
Idents(scrna) <- "Cell_type_fine_harmony"
unique(Idents(scrna)) # Selected subsets only 
scrna_subset <- subset(scrna, idents = c(#"Trm_ZNF683+", 
  "CD8_Teff_IFNG+", "CD8_Teff_GZMK+", "Tex_CXCL13+"))
#"Mono_FCN1+", "Macro_CCL4+", "Macro_ISGs+", "Macro_MTs+", #"Macro_LYVE1+", 
#"Macro_APOE+", 
#"Macro_Ki67+", "Promyelocytes", "Neu", "cDC1_CLEC9A+", "cDC2_CD1C+", "Mature_DC"

data <- GetAssayData(scrna_subset, assay = 'RNA', layer = 'counts')
cell_metadata <- scrna_subset@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(expression_dat = data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds) # Preprocessing data 
cds <- align_cds(cds, alignment_group = "sample") # Remove batch effects
cds <- reduce_dimension(cds) # Dimension reduction
cds <- cluster_cells(cds) # Cluster cells 
#cds <- learn_graph(cds) # Predict trajectory 
cds <- learn_graph(cds, use_partition = T, 
                   close_loop = T
) 

get_earliest_principal_node <- function(cds, time_bin="CD8_Teff_IFNG+"){
  cell_ids <- which(colData(cds)[, "Cell_type_fine_harmony"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

get_earliest_principal_node(cds)
#cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
cds <- order_cells(cds, root_pr_nodes = "Y_70")


#mycol[c("Tn_CCR7+", "Tmem_IL7R+", "Tmem_GPR183+FOSB+")]
#plot_cells_3d(cds_3d, color_cells_by = "Cell_type_fine_harmony", color_scale = mycol[c("Tn_CCR7+", "Tmem_IL7R+", "Tmem_GPR183+FOSB+")]) 
#plot_cells_3d(cds_3d, color_cells_by = "sample")
#plot_cells_3d(cds_3d, color_cells_by = "pseudotime")

mycol <- c(
  "#A7B3C6", "#4AAE79", "#D1A054", "#AC3E4E", "#8ED9E5", "#E89E2A", "#CED350", "#D94C43", 
  "#9A7FC1", "#A66C5B", "#60C8AE", "#B48A78", "#F6E27F", "#708CC2", "#F4A6B8", "#3A6CA8"
)
#mycol <- celltype_colors
names(mycol) <- levels(scrna$Cell_type_fine_harmony)
mycol <- mycol[levels(scrna$Cell_type_fine_harmony)]

p1 <- plot_cells(cds,
                 color_cells_by = "Cell_type_fine_harmony",
                 label_cell_groups=F,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size = 1.5, 
                 group_label_size = 2,
                 trajectory_graph_color = "#EAEAEA", 
                 trajectory_graph_segment_size = 0.5, 
                 cell_size = 0.5) + 
  scale_color_manual(values = mycol[c(#"Trm_ZNF683+", 
    "CD8_Teff_IFNG+", "CD8_Teff_GZMK+", "Tex_CXCL13+")]) + 
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.line.y = element_blank(), 
        panel.border = element_blank())
p1 <- p1 + coord_equal()

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_cell_type_plots_for_Tc.pdf"), 
    width = 3,  
    height = 3)  
print(p1, newpage = FALSE)
dev.off()

p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5, 
                 trajectory_graph_color = "#CCCCCC", 
                 trajectory_graph_segment_size = 0.5, cell_size = 0.5)
p2 <- p1 + theme(legend.position = "none",
                 axis.ticks.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.title.x = element_blank(), 
                 axis.line.x = element_blank(), 
                 axis.line.y = element_blank(), 
                 panel.border = element_blank()) +
  scale_color_gradientn(colours = c(
    "#EF88AD",
    "#A53860",
    "#670D2F",
    "#3A0519"
  ))

p2 <- p2 + coord_equal()
p2

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_Tc.pdf"), 
    width = 5,  
    height = 5)  
print(p2, newpage = FALSE)
dev.off()

##### Graph test  ----- 
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
#[1] "AC011462.1" "ACTB"       "CCL5"       "CREM"       "CXCL13"     "DNAJA1"     "DNAJB1"     "DUSP2"      "HSP90AA1"  
#[10] "HSP90AB1"   "HSPA1A"     "HSPA1B"     "HSPA8"      "HSPE1"      "HSPH1"      "IFITM1"     "IGFL2"      "JUN"       
#[19] "MALAT1"     "MT-CO1"     "NR4A1"      "RGCC"       "RPL41"      "YPEL5"     

##### Heatmap plot of key genes (w complext heatmap) ------ 
library(dplyr)
library(tidyr)
library(viridis)
pseudotime <- pseudotime(cds) %>% as.data.frame()
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "pseudotime"

celltype <- cds@colData$Cell_type_fine_harmony %>% as.data.frame()
celltype$cell <- colnames(cds)
colnames(celltype)[1] <- "clusters"

Treatment <- cds@colData$group_anno %>% as.data.frame()
Treatment$cell <- colnames(cds)
colnames(Treatment)[1] <- "treatment"

merge <- merge(pseudotime, celltype, by = 'cell')
merge <- merge(merge, Treatment, by = 'cell')

#devtools::install_github("junjunlab/ClusterGVis", force = T)
exp <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

exp <- data.frame(t(exprs(cds)))
exp <- exp[, c("IFNG", 
               #"CCL5",
               "CD69", 
               "DUSP2", 
               "FOSB",# "ID2", 
               "JUND",
               "CST7", #"CD27", 
               #"TNF", 
               "NKG7", 
               "GZMA", 
               "GZMK", 
               #"CD40LG", 
               "TIGIT", 
               "TOX",
               "CXCL13"
               #, 
               #"BCL2L1"
)]

heatmap_data_scaled <- data.frame(scale(exp)) 
heatmap_data_scaled$cell <- rownames(heatmap_data_scaled)
merge_mat <- merge(merge, heatmap_data_scaled, by = "cell")
merge_mat <- merge_mat[order(merge_mat$pseudotime), ]
merge_mat_for_plot <- as.data.frame(t(merge_mat)[c(5:16), ]) # The Gene Exp Matrix
for (i in 1:ncol(merge_mat_for_plot)) {
  merge_mat_for_plot[, i] <- as.numeric(merge_mat_for_plot[, i])
}

library(circlize)
range(merge_mat$pseudotime)[2]
col_fun = colorRamp2(seq(0, range(merge_mat$pseudotime)[2], 
                         length.out = 4), c(
                           "#EF88AD",
                           "#A53860",
                           "#670D2F",
                           "#3A0519"))

ht_list = HeatmapAnnotation(
  Pseudotime = anno_barplot(merge_mat$pseudotime, 
                            gp = gpar(fill = col_fun(merge_mat$pseudotime), col = col_fun(merge_mat$pseudotime)), 
                            height = unit(1.5, "cm")), 
  Celltype = merge_mat$clusters, 
  Treatment = merge_mat$treatment,
  col = list(Celltype = c(#"Trm_ZNF683+" = "#8ED9E5", 
    "Tmem_GPR183+FOSB+" = "#D1A054",
    "CD8_Teff_IFNG+" = "#E89E2A", "CD8_Teff_GZMK+" = "#CED350", 
    "Tex_CXCL13+" = "#D94C43"), 
    Treatment = c("H_Pre" = "#54426D", "H_Post" = "#D9A0B3", "L_Pre" = "#0f5688", "L_Post" = "#6B798E"))
)

cellheight = 0.5
rn = dim(as.matrix(merge_mat_for_plot))[1]
h=cellheight*rn
pdf(file = "Tc_pseudotime_gene_HM.pdf", 
    width = 6,  
    height = 6)  
Heatmap(as.matrix(#merge_mat[c("ATF3", "BCL2A1", "CD83", "FOS", "IL1B", "JUN", "NR4A1", "NR4A2", "NR4A3"), ]
  merge_mat_for_plot),
  # Height of the blocks
  #width = unit(w, "cm"),
  height = unit(h, "cm"),
  #scale = F,
  #rect_gp = gpar(col = "white", lwd = 1.5),
  #border_g = gpar(col = ,lty = 1,lwd = 1.2),
  # Dent formatting 
  #column_dend_height = unit(1.5, "cm"), 
  #row_dend_width = unit(1.5, "cm"),
  #column_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
  #row_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
  # Gap numbers and style 
  #row_split = 6, column_split = c("1", "1", "2", "2"),
  #row_gap = unit(2, "mm"),
  #column_gap = unit(2, "mm"),
  # Text formats 
  row_title = NULL,column_title = NULL,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  show_heatmap_legend = T, 
  cluster_rows = F,
  cluster_columns = F,
  row_names_side = 'right',
  show_column_names = F,
  show_row_names = T,
  border = T, 
  top_annotation = ht_list, 
  col = colorRamp2(breaks = seq(-0.75, 
                                3, length.out = 10), 
                   colors = viridis(10, alpha = 1, begin = 0, end = 1, direction = -1, option = "F"))
)
dev.off()

### NK ----- 
library(monocle3)
Idents(scrna) <- "Cell_type_fine_harmony"
unique(Idents(scrna)) # Selected subsets only 
scrna_subset <- subset(scrna, idents = c(#"Trm_ZNF683+", 
  "NK_KIT+", "NK_KLRC1+", "NK_FCGR3A+"))

data <- GetAssayData(scrna_subset, assay = 'RNA', layer = 'counts')
cell_metadata <- scrna_subset@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(expression_dat = data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds) # Preprocessing data 
cds <- align_cds(cds, alignment_group = "sample") # Remove batch effects
cds <- reduce_dimension(cds) # Dimension reduction
cds <- cluster_cells(cds) # Cluster cells 
#cds <- learn_graph(cds) # Predict trajectory 
cds <- learn_graph(cds, use_partition = T, 
                   close_loop = T
) 
#CD8_Teff_IFNG
get_earliest_principal_node <- function(cds, time_bin="NK_KIT+"){
  cell_ids <- which(colData(cds)[, "Cell_type_fine_harmony"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}


get_earliest_principal_node(cds)
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

mycol <- c(
  "#A7B3C6", "#4AAE79", "#D1A054", "#AC3E4E", "#8ED9E5", "#E89E2A", "#CED350", "#D94C43", 
  "#9A7FC1", "#A66C5B", "#60C8AE", "#B48A78", "#F6E27F", "#708CC2", "#F4A6B8", "#3A6CA8"
)
#mycol <- celltype_colors
names(mycol) <- levels(scrna$Cell_type_fine_harmony)
mycol <- mycol[levels(scrna$Cell_type_fine_harmony)]

p1 <- plot_cells(cds,
                 color_cells_by = "Cell_type_fine_harmony",
                 label_cell_groups=F,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size = 1.5, 
                 group_label_size = 2,
                 trajectory_graph_color = "#EAEAEA", 
                 trajectory_graph_segment_size = 0.5, 
                 cell_size = 0.5) + 
  scale_color_manual(values = mycol[c(#"Trm_ZNF683+", 
    "NK_KIT+", "NK_KLRC1+", "NK_FCGR3A+")]) + 
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.line.y = element_blank(), 
        panel.border = element_blank())
p1 <- p1 + coord_equal()
p1

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_cell_type_plots_for_NK.pdf"), 
    width = 3,  
    height = 3)  
print(p1, newpage = FALSE)
dev.off()

p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5, 
                 trajectory_graph_color = "#CCCCCC", 
                 trajectory_graph_segment_size = 0.5, cell_size = 0.5)
p2 <- p1 + theme(legend.position = "none",
                 axis.ticks.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.title.x = element_blank(), 
                 axis.line.x = element_blank(), 
                 axis.line.y = element_blank(), 
                 panel.border = element_blank()) +
  scale_color_gradientn(colours =  c("#4B3F72", "#EDEBE6", "#F2E86D"))

p2 <- p2 + coord_equal()
p2

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_NK.pdf"), 
    width = 3,  
    height = 3)  
print(p2, newpage = FALSE)
dev.off()

##### Graph test  ----- 
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
#[1] "AC011462.1" "ACTB"       "CCL5"       "CREM"       "CXCL13"     "DNAJA1"     "DNAJB1"     "DUSP2"      "HSP90AA1"  
#[10] "HSP90AB1"   "HSPA1A"     "HSPA1B"     "HSPA8"      "HSPE1"      "HSPH1"      "IFITM1"     "IGFL2"      "JUN"       
#[19] "MALAT1"     "MT-CO1"     "NR4A1"      "RGCC"       "RPL41"      "YPEL5"     

##### Heatmap plot of key genes (w complext heatmap) ------ 
pseudotime <- pseudotime(cds) %>% as.data.frame()
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "pseudotime"

celltype <- cds@colData$Cell_type_fine_harmony %>% as.data.frame()
celltype$cell <- colnames(cds)
colnames(celltype)[1] <- "clusters"

Treatment <- cds@colData$group_anno %>% as.data.frame()
Treatment$cell <- colnames(cds)
colnames(Treatment)[1] <- "treatment"

merge <- merge(pseudotime, celltype, by = 'cell')
merge <- merge(merge, Treatment, by = 'cell')

exp <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

exp <- data.frame(t(exprs(cds)))
exp <- exp[, c("KIT", 
               #"CCL5",
               #"NCAM1",
               "FXYD5", 
               "XCL1",
               "KLRC1",
               "FCER1G", 
               "GZMK",
               "CST7", 
               "NKG7", 
               "GZMA", 
               "FCGR3A",
               "FGFBP2",
               "PRF1"
               #, 
               #"BCL2L1"
)]

heatmap_data_scaled <- data.frame(scale(exp)) 
heatmap_data_scaled$cell <- rownames(heatmap_data_scaled)
merge_mat <- merge(merge, heatmap_data_scaled, by = "cell")
merge_mat <- merge_mat[order(merge_mat$pseudotime), ]
merge_mat_for_plot <- as.data.frame(t(merge_mat)[c(5:16), ]) # The Gene Exp Matrix
for (i in 1:ncol(merge_mat_for_plot)) {
  merge_mat_for_plot[, i] <- as.numeric(merge_mat_for_plot[, i])
}

library(circlize)
range(merge_mat$pseudotime)[2]
col_fun = colorRamp2(seq(0, range(merge_mat$pseudotime)[2], 
                         length.out = 3), c("#4B3F72", "#EDEBE6", "#F2E86D"))

ht_list = HeatmapAnnotation(
  Pseudotime = anno_barplot(merge_mat$pseudotime, 
                            gp = gpar(fill = col_fun(merge_mat$pseudotime), col = col_fun(merge_mat$pseudotime)), 
                            height = unit(1.5, "cm")), 
  Celltype = merge_mat$clusters, 
  Treatment = merge_mat$treatment,
  col = list(Celltype = c(#"Trm_ZNF683+" = "#8ED9E5", 
    "Tmem_GPR183+FOSB+" = "#D1A054",
    "NK_KIT+" = "#708CC2", "NK_KLRC1+" = "#F4A6B8", 
    "NK_FCGR3A+" = "#3A6CA8"), 
    Treatment = c("H_Pre" = "#54426D", "H_Post" = "#D9A0B3", "L_Pre" = "#0f5688", "L_Post" = "#6B798E"))
)


cellheight = 0.5
rn = dim(as.matrix(merge_mat_for_plot))[1]
h=cellheight*rn
pdf(file = "NK_pseudotime_gene_HM.pdf", 
    width = 6,  
    height = 6)  
Heatmap(as.matrix(#merge_mat[c("ATF3", "BCL2A1", "CD83", "FOS", "IL1B", "JUN", "NR4A1", "NR4A2", "NR4A3"), ]
  merge_mat_for_plot),
  # Height of the blocks
  #width = unit(w, "cm"),
  height = unit(h, "cm"),
  #scale = F,
  #rect_gp = gpar(col = "white", lwd = 1.5),
  #border_g = gpar(col = ,lty = 1,lwd = 1.2),
  # Dent formatting 
  #column_dend_height = unit(1.5, "cm"), 
  #row_dend_width = unit(1.5, "cm"),
  #column_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
  #row_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
  # Gap numbers and style 
  #row_split = 6, column_split = c("1", "1", "2", "2"),
  #row_gap = unit(2, "mm"),
  #column_gap = unit(2, "mm"),
  # Text formats 
  row_title = NULL,column_title = NULL,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  show_heatmap_legend = T, 
  cluster_rows = F,
  cluster_columns = F,
  row_names_side = 'right',
  show_column_names = F,
  show_row_names = T,
  border = T, 
  top_annotation = ht_list, 
  col = colorRamp2(breaks = seq(range(merge_mat_for_plot)[1]#-0.75
                                , 
                                3, length.out = 10), 
                   colors = viridis(10, alpha = 1, begin = 0, end = 1, direction = -1, option = "F"))
)
dev.off()

## MiloR (Figure S5F)----- 
### Use for all cells partitioned by pre and post status ----- 
library(miloR)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(ggbeeswarm)
library(scater)
library(scales)
library(forcats)
library(data.table)
library(stringr)
library(dplyr)

## R and NR with treatment statuses ----- 
## R ----- 
# Load data 
scrna <- readRDS("T_NK_final_24_06_25.rds")
Idents(scrna) <- "group_anno" 
scrna <- subset(scrna, idents = c("H_Pre", "H_Post"))
scrna$sample <- droplevels(scrna$sample)

scrna_sce <- as.SingleCellExperiment(scrna)
scrna_sce_milo <- miloR::Milo(scrna_sce) 
scrna_sce_milo <- miloR::buildGraph(scrna_sce_milo, k = 30, d = 50)
scrna_sce_milo <- makeNhoods(scrna_sce_milo, 
                             prop = 0.2,  
                             k = 30,  
                             d = 50,  
                             refined = TRUE)
scrna_sce_milo <- countCells(scrna_sce_milo, meta.data = data.frame(colData(scrna_sce_milo)), 
                             sample="sample")  
 
traj_design <- data.frame(colData(scrna_sce_milo))[,c("sample", "treatment_status")] 
traj_design$sample <- as.factor(traj_design$sample)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample
 
scrna_sce_milo <- calcNhoodDistance(scrna_sce_milo, d = 50)

da_results <- testNhoods(scrna_sce_milo, 
                         design = ~ treatment_status, 
                         design.df = traj_design) 
scrna_sce_milo <- buildNhoodGraph(scrna_sce_milo) 
da_results <- annotateNhoods(scrna_sce_milo, da_results, coldata_col = "Cell_type_fine_harmony")
da_results$Cell_type_fine_harmony <- factor(da_results$Cell_type_fine_harmony, levels = levels(scrna$Cell_type_fine_harmony))
da_results$Cell_type_fine_harmony <- factor(da_results$Cell_type_fine_harmony, 
                                            levels = rev(levels(da_results$Cell_type_fine_harmony)))
p4 <- plotDAbeeswarm(da_results, group.by = "Cell_type_fine_harmony") +
  scale_color_gradient2(low="#50859f", 
                        mid="darkgrey",
                        high="#d66692",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black')) 
p5 <- p4 + #+ scale_y_continuous(limits = c(-6, 6),                      
  #                  breaks = c(-6, -4, -2, 0, 2, 4, 6),                     
  #                  expand = c(0, 0)) +  
  labs(y = "") +  #coord_flip() + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme(legend.position = "none",        
        axis.title.x = element_text(size = 12),        
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),        
        axis.title.y = element_blank(),        
        axis.ticks.y = element_blank(),        
        axis.text.y = element_blank(),        
        axis.line.y = element_blank(),        
        plot.margin = unit(c(0, 1, 0, 0), "cm")) +  
  ylab("R")
#ylab(bquote(Log[2] ~ italic('Fold Change')))
p5

saveRDS(da_results, "da_results_post_vs_pre_R_TNK_260125.rds")

### Use for all cells partitioned by efficacy status ----- 
scrna <- readRDS("T_NK_final_24_06_25.rds")
Idents(scrna) <- "group_anno" 
scrna <- subset(scrna, idents = c("L_Pre", "L_Post"))
scrna$sample <- droplevels(scrna$sample)

scrna_sce <- as.SingleCellExperiment(scrna)
scrna_sce_milo <- miloR::Milo(scrna_sce)#milo object构建
scrna_sce_milo <- miloR::buildGraph(scrna_sce_milo, k = 30, d = 50)
scrna_sce_milo <- makeNhoods(scrna_sce_milo, 
                             prop = 0.2, 
                             k = 30, 
                             d = 50,
                             refined = TRUE)
scrna_sce_milo <- countCells(scrna_sce_milo, meta.data = data.frame(colData(scrna_sce_milo)), 
                             sample="sample") 

traj_design <- data.frame(colData(scrna_sce_milo))[,c("sample", "treatment_status")]
traj_design$sample <- as.factor(traj_design$sample)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample

scrna_sce_milo <- calcNhoodDistance(scrna_sce_milo, d = 50)

da_results <- testNhoods(scrna_sce_milo, 
                         design = ~ treatment_status, 
                         design.df = traj_design)
scrna_sce_milo <- buildNhoodGraph(scrna_sce_milo)
da_results <- annotateNhoods(scrna_sce_milo, da_results, coldata_col = "Cell_type_fine_harmony")
da_results$Cell_type_fine_harmony <- factor(da_results$Cell_type_fine_harmony, levels = rev(levels(da_results$Cell_type_fine_harmony)))
p6 <- plotDAbeeswarm(da_results, group.by = "Cell_type_fine_harmony") +
  scale_color_gradient2(low="#50859f", 
                        mid="darkgrey",
                        high="#d66692",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black')) 
p7 <- p6 + #+ scale_y_continuous(limits = c(-6, 6),                      
  #                  breaks = c(-6, -4, -2, 0, 2, 4, 6),                     
  #                  expand = c(0, 0)) +  
  labs(y = "") +  #coord_flip() + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme(legend.position = "none",        
        axis.title.x = element_text(size = 12),        
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),        
        axis.title.y = element_blank(),        
        axis.ticks.y = element_blank(),        
        axis.text.y = element_blank(),        
        axis.line.y = element_blank(),        
        plot.margin = unit(c(0, 1, 0, 0), "cm")) +  
  ylab("NR")
#ylab(bquote(Log[2] ~ italic('Fold Change')))
p7

saveRDS(da_results, "da_results_post_vs_pre_NR_TNK_260125.rds")

comb <- p5 + p7
comb

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Combined_TNK_MiloR_260125", "_v1.pdf"), 
    width = 5,  
    height = 5)  
print(comb, newpage = FALSE)
dev.off()

