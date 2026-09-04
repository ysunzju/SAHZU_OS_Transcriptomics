library(here) # project-root-relative paths; run scripts from repository root
##################### Endothelial cell pipeline ##################### 
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/Endo"))
library(Seurat)
library(tidyverse)

scrna <- readRDS(here("data/28_04_25_OS_all_cell_plots/Endo_reident_ordered_01_05_25.rds"))
scrna[["RNA"]]$scale.data <- NULL
scrna[["RNA"]]$data <- NULL
scrna[["RNA"]] <- split(scrna[["RNA"]], f = scrna$sample) 
scrna[["RNA"]] 

# Endo ----- 
# Pipeline with regressing out ------ 
DefaultAssay(scrna) <- "RNA"
scrna <- NormalizeData(scrna) %>% FindVariableFeatures() 
#scrna[["percent.ribo"]] <- PercentageFeatureSet(object = scrna, pattern = "^RP[SL]")
#Variable_features <- scrna@assays[["RNA"]]@meta.data[["var.features"]][!is.na(scrna@assays[["RNA"]]@meta.data[["var.features"]])]
#MT_genes <- rownames(scrna)[grep("^MT_", rownames(scrna), invert=F)] 
#RB_genes <- rownames(scrna)[grep("^RP[SL]", rownames(scrna), invert = F)]
#genes_to_regress_out <- c(MT_genes, RB_genes)
#Variable_features[Variable_features %in% genes_to_regress_out]
#Variable_features_to_keep <- Variable_features[!Variable_features %in% genes_to_regress_out]
#scrna <- ScaleData(scrna, vars.to.regress = c("percent.mt", "percent.ribo"))

scrna <- ScaleData(scrna)
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

# Clustree 
library(clustree)
clustree_plot <- clustree(scrna@meta.data, prefix = "RNA_snn_res."); clustree_plot # Manual saves

##### Clustering with SNN res of 1.5 #####
DimPlot(scrna, group.by = "RNA_snn_res.1.5", reduction = "umap", split.by = "treatment_status", label = T)
#Idents(scrna) <- "RNA_snn_res.2.1"
Idents(scrna) <- "RNA_snn_res.1.5"
DimPlot(scrna, label = T) # Manual saves
FeaturePlot(scrna, "RGCC")
# PECAM1 and VWF as general markers for endothelial cells 
# FLT1 and EMCN and more of a specific marker for same?
# Arteries - GJA4, GJA5, FBLN5 
# Tip-cell - COL4A1, KDR (VEGFR2), ESM1 
# Veins - ACKR1, SELP, CLU 
# Cap - CA4, CD36, RGCC
# Art - EFNB2, SEMA3G
# Art - IGFBP3 maybe also 
# If RGCC expression, mark as cap even if vein marker expressions (ie ACKR1)

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
organism = 'hsa'     
OrgDb = 'org.Hs.eg.db'

markers <- FindAllMarkers(scrna,
                          logfc.threshold = 0.25, 
                          min.diff.pct = 0.2, 
                          only.pos = T)

# Adding tip and cell signatures 
Tip_sig <- list(c("ADM", "ANGPT2", "ANKRD37", "APLN", "C1QTNF6", "CD93", "CLDN5", "COL4A1", "COL4A2", 
                  "COTL1", "CXCR4", "DLL4", "EDNRB", "ESM1", "FSCN1", "GPIHBP1", "HSPG2", "IGFBP3", 
                  "INHBB", "ITGA5", "JUP", "KCNE3", "KCNJ8", "KDR", "LAMA4", "LAMB1", "LAMC1", "LXN", 
                  "MARCKS", "MARCKSL1", "MCAM", "MEST", "MYH9", "MYO1B", "N4BP3", "NID2", "NOTCH4", 
                  "PDGFB", "PGF", "PLOD1", "PLXND1", "PMEPA1", "PTN", "RAMP3", "RBP1", "RGCC", "RHOC", 
                  "SMAD1", "SOX17", "SOX4", "SPARC", "TCF4", "UNC5B", "VIM"))
scrna <- AddModuleScore(scrna, features = Tip_sig, name = "Tip", assay = "RNA") 
names(scrna@meta.data)[47] <- "Tip"
Stalk_sig <- list(c("ACKR1", "AQP1", "C1QTNF9", "CD36", "CSRP2", "EHD4", "FBLN5", "HSPB1", "LIGP1", 
                    "IL6ST", "JAM2", "LGALS3", "LRG1", "MEOX2", "PLSCR2", "CAVIN2", "SELP", "SPINT2", 
                    "TGFBI", "TGM2", "TMEM176A", "TMEM176B", "TMEM252", "TSPAN7", "FLT1", "VWF"))
scrna <- AddModuleScore(scrna, features = Stalk_sig, name = "Stalk", assay = "RNA") 
names(scrna@meta.data)[48] <- "Stalk"
#scrna@meta.data[38] <- NULL

scrna$Tip_like <- scrna$Tip - scrna$Stalk
scrna$Stalk_like <- scrna$Stalk - scrna$Tip

##### 0 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "0", only.pos = F, pct.min = 0.25)
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

scrna <- RenameIdents(scrna, "0" = "capEC") # 

##### 1 #####
#devtools::install_github('immunogenomics/presto')
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

scrna <- RenameIdents(scrna, "1" = "capEC") 

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

scrna <- RenameIdents(scrna, "2" = "venEC") # 

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

scrna <- RenameIdents(scrna, "3" = "venEC") # 

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

scrna <- RenameIdents(scrna, "4" = "capEC") # 

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

scrna <- RenameIdents(scrna, "5" = "capEC") # 

##### 6 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "6", only.pos = F)
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

scrna <- RenameIdents(scrna, "6" = "capEC") 

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

scrna <- RenameIdents(scrna, "7" = "capEC") 

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

scrna <- RenameIdents(scrna, "8" = "capEC")

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

scrna <- RenameIdents(scrna, "9" = "artEC") 

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

scrna <- RenameIdents(scrna, "10" = "capEC") 

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

scrna <- RenameIdents(scrna, "11" = "artEC") #

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

scrna <- RenameIdents(scrna, "12" = "Likely_mye_doub") 

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

scrna <- RenameIdents(scrna, "13" = "capEC") # 

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

scrna <- RenameIdents(scrna, "14" = "artEC") 

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

scrna <- RenameIdents(scrna, "15" = "capEC") 

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

scrna <- RenameIdents(scrna, "16" = "capEC") # 

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

scrna <- RenameIdents(scrna, "17" = "EndoMT-II") 
#scrna <- RenameIdents(scrna, "EndoMT-II" = "EndoMT-II_C01") 

##### 18 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "18", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 18)][1:50]

scrna <- RenameIdents(scrna, "18" = "Cycling_ECs") #  

##### 19 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "19", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 19)][1:50]

scrna <- RenameIdents(scrna, "19" = "artEC") # 

##### 20 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "20", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 20)][1:50]

scrna <- RenameIdents(scrna, "20" = "Cycling_ECs") 

##### 21 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "21", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 21)][1:50]

scrna <- RenameIdents(scrna, "21" = "capEC") # Cycling

##### 22 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "22", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 22)][1:50]

scrna <- RenameIdents(scrna, "22" = "capEC")

##### 23 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "23", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 23)][1:50]

scrna <- RenameIdents(scrna, "23" = "EndoMT-I")

##### 24 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "24", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 24)][1:50]

scrna <- RenameIdents(scrna, "24" = "venEC") 

##### 25 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "25", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 25)][1:50]

scrna <- RenameIdents(scrna, "25" = "Likely_T_doub") # 

##### 26 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "26", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 26)][1:50]

scrna <- RenameIdents(scrna, "26" = "LEC") 

##### 27 #####
#devtools::install_github('immunogenomics/presto')
DefaultAssay(scrna) <- "RNA"
need_DEG <- FindMarkers(scrna, ident.1 = "27", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C8)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 27)][1:50]

scrna <- RenameIdents(scrna, "27" = "Likely_multiplets") 

## Subset and Attach labels ----- 
levels(scrna)
#[1] "Likely_multiplets" "capEC"             "LEC"               "Likely_T_doub"    
#[5] "venEC"             "EndoMT-I"          "Cycling_ECs"       "artEC"            
#[9] "EndoMT-II"         "Likely_mye_doub"  
to_keep <- levels(scrna)[levels(scrna) %in% c("capEC", "LEC", "venEC", 
                                              "EndoMT-I", "Cycling_ECs", "artEC", 
                                              "EndoMT-II")]
scrna <- subset(scrna, idents = to_keep)
scrna$Cell_type_fine_harmony <- Idents(scrna)
#saveRDS(scrna, "All_myeloid_reident_updated_27_05_25.rds")

saveRDS(scrna, "All_endo_reident_updated_251227.rds") 
FeaturePlot(scrna, c("Tip", "Stalk", "Tip_like", "Stalk_like"))

#### Cell prop plots -----
# FLT1 and EMCN and more of a specific marker for same?
# Arteries - GJA4, GJA5, FBLN5 
# Tip-cell - COL4A1, KDR (VEGFR2), ESM1 
# Veins - ACKR1, SELP, CLU 
# Cap - CA4, CD36, RGCC
# Art - EFNB2, SEMA3G
# Art - IGFBP3 maybe also 

# Feautureplots ----- 
FeaturePlot(scrna, c("GJA4", "GJA5", "FBLN5", 
                     "COL4A1", "KDR", "ESM1", 
                     "ACKR1", "SELP", "CLU", 
                     "CA4", "CD36", "RGCC"), ncol = 3, raster = F)

scrna$Cell_type_fine_harmony <- 
  factor(scrna$Cell_type_fine_harmony, 
         levels = c("artEC", "capEC", "venEC", "EndoMT-I", "EndoMT-II", "LEC", "Cycling_ECs"))

dim(scrna_sub)
dim(scrna)
Idents(scrna) <- scrna$Cell_type_fine_harmony
DimPlot(scrna, label = T)
#saveRDS(scrna_sub, "Endo_reident_updated_doublets_removed_31_05_25.rds") 

# Re UMAP
scrna[["RNA"]]$scale.data <- NULL
scrna[["RNA"]]$data <- NULL
scrna[["RNA"]] <- split(scrna[["RNA"]], f = scrna$sample) 
scrna[["RNA"]] 

DefaultAssay(scrna) <- "RNA"
scrna <- NormalizeData(scrna) %>% FindVariableFeatures() %>% ScaleData()
scrna <- RunPCA(scrna)
scrna <- FindNeighbors(scrna, dims = 1:50, reduction = "pca")
scrna <- FindClusters(scrna, resolution = 1, cluster.name = "unintegrated_clusters")
scrna <- RunUMAP(scrna, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
#DimPlot(scrna, reduction = "umap.unintegrated", group.by = c("treatment_status", "unintegrated_clusters"))
#DimPlot(scrna, group.by = "sample", reduction = "umap.unintegrated", split.by = "treatment_status")
scrna <- IntegrateLayers(object = scrna, method = HarmonyIntegration, 
                         orig.reduction = "pca", new.reduction = "integrated.harmony",
                         verbose = TRUE)
scrna[["RNA"]] <- JoinLayers(scrna[["RNA"]])

scrna <- FindNeighbors(scrna, dims = 1:50, reduction = "integrated.harmony")
scrna <- RunUMAP(scrna, dims = 1:50, reduction = "integrated.harmony") #umap as the new reduction
Idents(scrna) <- "Cell_type_fine_harmony"
DimPlot(scrna, label = T)
saveRDS(scrna, "Final_endo_251227.rds")

##################### Endothelial cell figure panels (related to Figure 4 and Figure S4) ##################### 
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/Endo"))

## Endo UMAP (Figure 4A) ----- 
library(Seurat)
Endo <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Endo/Final_endo_251227.rds"))
#DimPlot(Endo, group.by = "Cell_type_fine_harmony")

colors <- c("#C43C39", "#bb82b1", "#9cd2ed","#ea9994", "#bc9a7f",# d4de9c
            "#94c58f", '#e0cfda') 
mycol <- colors
names(mycol) <- levels(Endo$Cell_type_fine_harmony)

p1 <- DimPlot(Endo, reduction = "umap", group.by = c("Cell_type_fine_harmony"), 
              label = T, raster = T, pt.size = 5,raster.dpi = c(2048, 2048)) + 
  scale_color_manual(values = colors) + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        axis.title = element_blank(), 
        plot.title = element_blank()) + 
  coord_fixed(ratio = 1) + 
  NoLegend()
p1

pdf(paste0("Endo_cells", "_UMAP_v251224.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

## Marker dot plot (Figure S4A) ----- 
mar <- FindAllMarkers(Endo, only.pos = T)
levels(Endo$Cell_type_fine_harmony)
mar$gene[mar$cluster == "artEC"][1:50]
mar$gene[mar$cluster == "capEC"][1:50]
mar$gene[mar$cluster == "venEC"][1:50]
mar$gene[mar$cluster == "EndoMT-I"][1:50]
mar$gene[mar$cluster == "EndoMT-II"][1:50]
mar$gene[mar$cluster == "LEC"][1:50]

Features <- c("GJA4", "GJA5", "EFNB2", "RGCC", "ACKR1", "SELP", "RGS5", "PDGFRB", "NOTCH3", 
              "PCOLCE", "DCN", "LUM", "CCL21", "PDPN", "PROX1", "TOP2A", "MKI67")

#Idents(DC) <- "Cell_ident_for_MHC_plot"
p <- DotPlot(Endo, Features) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1 <- ggplot(p$data, aes(x = features.plot, y = id)) +   
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +   
  #facet_grid(facets = ~feature.groups,  switch = "x", scales = "free_x", space = "free_x") +    
  scale_radius(breaks = c(25, 50, 75, 100), range = c(0,6)) +   
  theme_classic() + 
  coord_equal() + 
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
pdf(file = paste0(time, "_", "All_EC_subs_sig_genes.pdf"), 
    width = 9, # The width of the plot in inches
    height = 3) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()

## Endothelial cell proportion (Figure 4B) ------
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

Idents(Endo) <- "group"
Endo <- RenameIdents(Endo, "A" = "H_Pre", 
                     "B" = "H_Post", 
                     "C" = "L_Pre", 
                     "D" = "L_Post")
Endo$group_anno <- Idents(Endo)
Endo$Cell_type_fine_harmony <- droplevels(Endo$Cell_type_fine_harmony)
Cellratio <- prop.table(table(Endo$Cell_type_fine_harmony, Endo$sample), margin = 2) 
Cellratio <- data.frame(Cellratio)

cellper <- reshape2::dcast(Cellratio, Var2~Var1, value.var="Freq")
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]

#write.csv()
#to_save <- t(cellper)
#write.csv(to_save, file='all_endo_cells_prop_260803.csv', quote = F)

meta <- Endo@meta.data
colnames(meta)
meta <- meta[,c(32,4)] # group and sample information
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
#cellper$group <- NA
meta <- unique(meta)
cellper <- merge(cellper, meta, by = "sample")
colnames(cellper)
colnames(cellper)[9] <- "group" # Group anno to group 
cellper <- as.data.frame(cellper)

pplist =list()
seuratObj_groups = unique(levels(Endo$Cell_type_fine_harmony))

mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
#fill_color <- zzm60colors
#names(fill_color) <- unique(levels(seuratObj$cluster))

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
  #p_position2 = min(max(cellper_$percent)*1.2) * 0.8
  
  print(group_)
  print(cellper_$median)
  mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + 
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
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) + 
    theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid"))
  
  ###组间t检验分析/wilcox.test
  labely = max(cellper_$percent)
  compare_means(percent ~group,  data = cellper_)
  my_comparisons <-list(c("H_Pre","H_Post") 
                        #c("H_Pre","L_Pre"), 
                        #c("L_Pre","L_Post")
  )
  pp1 = pp1 + geom_signif(comparisons = list(c("H_Pre","H_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position, 
                          #paired = T, 
                          #test = "t.test",
                          test.args = list(paired = TRUE, 
                                           alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  #my_comparisons1 <-list(c("L_Pre","L_Post"))
  pp1 = pp1 + geom_signif(comparisons = list(c("L_Pre","L_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))}, 
                          y_position = p_position, 
                          #paired = T, 
                          #test = "t.test",
                          test.args = list(paired = TRUE, 
                                           alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("H_Pre","L_Pre")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position1, 
                          #paired = T, 
                          #test = "t.test", 
                          test.args   = list(#paired = TRUE, 
                            alternative = "two.sided"),
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
                 #pplist[[8]],
                 #pplist[[9]],
                 #pplist[[10]],
                 #pplist[[11]], 
                 #pplist[[12]],
                 #pplist[[13]],
                 #pplist[[14]],
                 #pplist[[15]],
                 #pplist[[16]], 
                 ncol = 4 #, nrow = 2
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_Endo_cell_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 6, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches
print(pps, newpage = FALSE)
dev.off()

## VEC trajectory (Figure 4C,D) ----- 
library(monocle3)
Idents(Endo) <- "Cell_type_fine_harmony"
unique(Idents(Endo)) # Selected subsets only 
scrna_subset <- subset(Endo, idents = c("artEC", "capEC", "venEC"))
#Idents(scrna_subset) <- "sample"
#scrna_subset <- subset(scrna_subset, ident = c("H_BF_P1", "H_AF_P1"))
scrna_subset$Cell_type_fine_harmony <-droplevels(scrna_subset$Cell_type_fine_harmony)
scrna_subset$group <-droplevels(scrna_subset$group)
scrna_subset$sample <-droplevels(scrna_subset$sample)

unique(scrna_subset$Cell_type_fine_harmony); unique(scrna_subset$group)

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
cds <- learn_graph(cds, use_partition = F, 
                   close_loop = F) 
get_earliest_principal_node <- function(cds, time_bin="artEC"){
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
colors <- c("#C43C39", "#bb82b1", "#9cd2ed","#ea9994", "#bc9a7f",# d4de9c
            "#94c58f", '#e0cfda') 
mycol <- colors
names(mycol) <- levels(scrna_subset$Cell_type_fine_harmony)

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
  scale_color_manual(values = mycol[c("artEC", "capEC", "venEC")]) + 
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
p1 + coord_equal()

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_cell_type_plots_for_VEC.pdf"), 
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
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
  scale_color_gradientn(colours = c("#B36A6A",  "#C99191",  "#D9B0A0",  "#E8C7B6",  "#D2C0CF",  "#B8A8C6",  "#9B8AB3", "#7A6E9C"))
#                          c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124"))
p2

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_VEC.pdf"), 
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(p2, newpage = FALSE)
dev.off()

##### Graph test  ----- 
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
saveRDS(cds, "VEC_cds_obj_251230.rds") 

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

Treatment <- cds@colData$treatment_status %>% as.data.frame()
Treatment$cell <- colnames(cds)
colnames(Treatment)[1] <- "treatment"

merge <- merge(pseudotime, celltype, by = 'cell')
merge <- merge(merge, Treatment, by = 'cell')

#devtools::install_github("junjunlab/ClusterGVis", force = T)
library(ClusterGVis)
exp <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

exp <- data.frame(t(exprs(cds)))
exp <- exp[, c("GJA4", "RGCC", "ACKR1")]

library(dplyr)
heatmap_data_scaled <- data.frame(scale(exp)) 

library(ComplexHeatmap)
library(circlize)
heatmap_data_scaled$cell <- rownames(heatmap_data_scaled)
merge_mat <- merge(merge, heatmap_data_scaled, by = "cell")
merge_mat <- merge_mat[order(merge_mat$pseudotime), ]
merge_mat_for_plot <- as.data.frame(t(merge_mat)[c(5:7), ]) # The Gene Exp Matrix
for (i in 1:ncol(merge_mat_for_plot)) {
  merge_mat_for_plot[, i] <- as.numeric(merge_mat_for_plot[, i])
}

library(circlize)
range(merge_mat$pseudotime)
col_fun = colorRamp2(seq(0, 30, length.out = 8), c("#B36A6A",  "#C99191",  "#D9B0A0",  "#E8C7B6",  "#D2C0CF",  "#B8A8C6",  "#9B8AB3", "#7A6E9C"))

ht_list = HeatmapAnnotation(
  Pseudotime = anno_barplot(merge_mat$pseudotime, 
                            gp = gpar(fill = col_fun(merge_mat$pseudotime), col = col_fun(merge_mat$pseudotime)), 
                            height = unit(2.1, "cm")), 
  Celltype = merge_mat$clusters, 
  Treatment = merge_mat$treatment,
  col = list(Celltype = c("artEC" = "#C43C39", "capEC" = "#bb82b1", "venEC" = "#9cd2ed"), 
             Treatment = c("Pre" = "#90486e", "Post" = "#efcfe3"))
)

cellheight = 0.5
rn = dim(as.matrix(merge_mat_for_plot))[1]
h=cellheight*rn
pdf(file = "VEC_pseudotime_gene_HM.pdf", 
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
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
  col = colorRamp2(breaks = seq(-0.5, 2.5, length.out = 10), 
                   colors = viridis(10, alpha = 1, begin = 0, end = 1, direction = -1, option = "F") )
)
dev.off()

# Cap EC pipeline ----- 
scrna_cap <- subset(Endo, idents = "capEC")
scrna_cap[["RNA"]]$scale.data <- NULL
scrna_cap[["RNA"]]$data <- NULL
scrna_cap[["RNA"]] <- split(scrna_cap[["RNA"]], f = scrna_cap$sample) 

DefaultAssay(scrna_cap) <- "RNA"
scrna_cap <- NormalizeData(scrna_cap) %>% FindVariableFeatures() 
scrna_cap <- ScaleData(scrna_cap)
scrna_cap <- RunPCA(scrna_cap)
scrna_cap <- FindNeighbors(scrna_cap, dims = 1:50, reduction = "pca")
scrna_cap <- FindClusters(scrna_cap, resolution = 1, cluster.name = "unintegrated_clusters")
scrna_cap <- RunUMAP(scrna_cap, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(scrna_cap, reduction = "umap.unintegrated", group.by = c("treatment_status", "unintegrated_clusters"))
DimPlot(scrna_cap, group.by = "sample", reduction = "umap.unintegrated", split.by = "treatment_status")
scrna_cap <- IntegrateLayers(object = scrna_cap, method = HarmonyIntegration, 
                             orig.reduction = "pca", new.reduction = "integrated.harmony",
                             verbose = TRUE)
scrna_cap[["RNA"]] <- JoinLayers(scrna_cap[["RNA"]])

scrna_cap <- FindNeighbors(scrna_cap, dims = 1:50, reduction = "integrated.harmony")
scrna_cap <- RunUMAP(scrna_cap, dims = 1:50, reduction = "integrated.harmony") #umap as the new reduction
scrna_cap <- FindClusters(object = scrna_cap, resolution = c(0.01, 0.05, 0.08, seq(.1,2.1,.1)))

Tip_sig <- list(c("ADM", "ANGPT2", "ANKRD37", "APLN", "C1QTNF6", "CD93", "CLDN5", "COL4A1", "COL4A2", 
                  "COTL1", "CXCR4", "DLL4", "EDNRB", "ESM1", "FSCN1", "GPIHBP1", "HSPG2", "IGFBP3", 
                  "INHBB", "ITGA5", "JUP", "KCNE3", "KCNJ8", "KDR", "LAMA4", "LAMB1", "LAMC1", "LXN", 
                  "MARCKS", "MARCKSL1", "MCAM", "MEST", "MYH9", "MYO1B", "N4BP3", "NID2", "NOTCH4", 
                  "PDGFB", "PGF", "PLOD1", "PLXND1", "PMEPA1", "PTN", "RAMP3", "RBP1", "RGCC", "RHOC", 
                  "SMAD1", "SOX17", "SOX4", "SPARC", "TCF4", "UNC5B", "VIM"))
scrna_cap <- AddModuleScore(scrna_cap, features = Tip_sig, name = "Tip_cap_only", assay = "RNA") 
names(scrna_cap@meta.data)[51] <- "Tip_cap_only"
Stalk_sig <- list(c("ACKR1", "AQP1", "C1QTNF9", "CD36", "CSRP2", "EHD4", "FBLN5", "HSPB1", "LIGP1", 
                    "IL6ST", "JAM2", "LGALS3", "LRG1", "MEOX2", "PLSCR2", "CAVIN2", "SELP", "SPINT2", 
                    "TGFBI", "TGM2", "TMEM176A", "TMEM176B", "TMEM252", "TSPAN7", "FLT1", "VWF"))
scrna_cap <- AddModuleScore(scrna_cap, features = Stalk_sig, name = "Stalk_cap_only", assay = "RNA") 
names(scrna_cap@meta.data)[52] <- "Stalk_cap_only"

Cap_only_mar <- FindAllMarkers(scrna_cap, group.by = "RNA_snn_res.0.2", only.pos = T)
DimPlot(scrna_cap, group.by = "sample") 
FeaturePlot(scrna_cap, "Tip") 
FeaturePlot(scrna_cap, "Tip_cap_only") 
FeaturePlot(scrna_cap, "Stalk") 
FeaturePlot(scrna_cap, "Stalk_cap_only")

## Tip vs stalk-like cap EC signature heatmap (Figure 4E) -----
library(tidyverse)
library(dplyr)
#dplyr::select
genelistnames <- c("Tip", "Stalk")
heatmap_data <- scrna_cap@meta.data %>%
  dplyr::select(group_anno, all_of(genelistnames)) %>%
  dplyr::group_by(group_anno) %>%
  dplyr::summarise(across(all_of(genelistnames), mean, na.rm = TRUE)) %>%
  column_to_rownames("group_anno") %>%
  as.matrix() %>%
  t()  

heatmap_data_scaled <- t(scale(t(heatmap_data)))

p_value_matrix <- scrna_cap@meta.data %>%
  dplyr::group_by(group_anno) %>%
  dplyr::summarise(across(all_of(genelistnames), ~t.test(.x)$p.value)) %>%
  column_to_rownames("group_anno") %>%
  as.matrix() %>%
  t()

library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

cellwidth = 0.7
cellheight = 0.7
cn = dim(as.matrix(heatmap_data_scaled))[2]
rn = dim(as.matrix(heatmap_data_scaled))[1]
w=cellwidth*cn
h=cellheight*rn

pdf("Tip_Stalk_sig_by_group.pdf", width = 3, height = 3)
Heatmap(as.matrix(heatmap_data_scaled),
        # Height of the blocks
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        #scale = F,
        #rect_gp = gpar(col = "white", lwd = 1.5),
        #border_g = gpar(col = ,lty = 1,lwd = 1.2),
        # Dent formatting 
        column_dend_height = unit(1.5, "cm"), 
        row_dend_width = unit(1.5, "cm"),
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
        row_names_side = 'right',
        show_column_names = T,
        show_row_names = T,
        border = T, 
        col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#50859f","white","#d66692")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if( heatmap_data_scaled[i, j] > 1) {
            grid.text("+++", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] > 0.5 & heatmap_data_scaled[i, j] <= 1) {
            grid.text("++", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= 0.5 & heatmap_data_scaled[i, j] > 0.2) {
            grid.text("+", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= 0.2 & heatmap_data_scaled[i, j] > -0.2) {
            grid.text("+/-", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= -0.2 & heatmap_data_scaled[i, j] > -0.5) {
            grid.text("-", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= -0.5 & heatmap_data_scaled[i, j] > -1) {
            grid.text("--", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= -1 & heatmap_data_scaled[i, j] > -2) {
            grid.text("---", x, y, gp = gpar(fontsize = 8))
          }
        }
)
dev.off()

## Tip and stalk feature density plot (related to Figure 4F) ----- 
library(Nebulosa)
p <- plot_density(scrna_cap, c("Tip", "Stalk"), reduction = "umap", combine = F, raster = T
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

f = p[[1]] + p[[2]]
f

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "FeaturePlot_Endo_cluster_marker_sig.pdf"), 
    width = 6, # The width of the plot in inches
    height = 3) # The height of the plot in inches
print(f, newpage = FALSE)
dev.off()

## CapEC DimPlot by treatment and response (Figure 4H) ----- 
library(SCP)
p1 <- CellDimPlot(
  srt = scrna_cap, group.by = "group_anno", 
  reduction = "umap", theme_use = "theme_blank", label = T, label_repel = T, 
  palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")), 
  raster = T, 
  pt.size = 3, 
  split.by = "group_anno", label_insitu = T) + theme(legend.position = "none")
p1

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "DimPlot_by_Rx_and_Res.pdf"), 
    width = 6, 
    height = 6) 
print(p1, newpage = FALSE)
dev.off()
          
## CapEC subset annotation -----
#Idents(scrna_cap) <- "RNA_snn_res.0.3"
Idents(scrna_cap) <- "RNA_snn_res.0.2"
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

hs_C5 = msigdbr(species = "Homo sapiens",
                category = "C5",
                subcategory = NULL) %>% as.data.frame() %>% 
  dplyr::select(gs_cat, gs_subcat, gs_name, gene_symbol)

hs_C5 = hs_C5 %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  dplyr::rename("term" = "gs_name", "gene" = "gene_symbol")

rm(hs_df)

library(org.Hs.eg.db)
organism = 'hsa'   
OrgDb = 'org.Hs.eg.db'

markers <- FindAllMarkers(scrna_cap,
                          #logfc.threshold = 0.25, 
                          min.diff.pct = 0.2, 
                          only.pos = T)

###### 0 ------
DefaultAssay(scrna_cap) <- "RNA"
need_DEG <- FindMarkers(scrna_cap, ident.1 = "0", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C5)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 0)][1:50]

scrna_cap <- RenameIdents(scrna_cap, "0" = "Transition_FLT1hi") 

###### 1 ------
DefaultAssay(scrna_cap) <- "RNA"
need_DEG <- FindMarkers(scrna_cap, ident.1 = "1", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C5)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 1)][1:50] 

scrna_cap <- RenameIdents(scrna_cap, "1" = "Transition_FLT1lo") 

###### 2 ------
DefaultAssay(scrna_cap) <- "RNA"
need_DEG <- FindMarkers(scrna_cap, ident.1 = "2", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C5)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2
need_DEG$Symbol <- rownames(need_DEG)
markers$gene[which(markers$cluster == 2)][1:50]

scrna_cap <- RenameIdents(scrna_cap, "2" = "Stalk-like") 

###### 3 ------ 
DefaultAssay(scrna_cap) <- "RNA"
need_DEG <- FindMarkers(scrna_cap, ident.1 = "3", only.pos = F)
geneList <- need_DEG$avg_log2FC
names(geneList) <- rownames(need_DEG)
geneList <- sort(geneList, decreasing = T)    
egmt <- GSEA(geneList, TERM2GENE = hs_C5)
gsea_results <- egmt@result
gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]
View(gsea_results) # 
need_DEG$d_pct <- need_DEG$pct.1 - need_DEG$pct.2 
need_DEG$Symbol <- rownames(need_DEG) 
markers$gene[which(markers$cluster == 3)][1:50] 

scrna_cap <- RenameIdents(scrna_cap, "3" = "Tip-like") 

###### 4 ------
DefaultAssay(scrna_cap) <- "RNA"
need_DEG <- FindMarkers(scrna_cap, ident.1 = "4", only.pos = F)
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

scrna_cap <- RenameIdents(scrna_cap, "4" = "Stalk-like") 

scrna_cap$Cell_type_fine_harmony <- Idents(scrna_cap)
saveRDS(scrna_cap, "CapEC_updated_251231.rds")
#Idents(scrna_cap) <- "Cell_type_fine_harmony"
#scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Endo/All_endo_reident_updated_251227.rds"))
#Idents(scrna) <- "Cell_type_fine_harmony"
#Idents(scrna) <- Idents(scrna_cap)
#scrna$Cell_type_fine_harmony <- Idents(scrna)
#saveRDS(scrna, "EC_w_Cap_updated_251231.rds")

#### CapEC cellular proportions (related to Figure 4G)----- 
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

Idents(scrna_cap) <- "group"
scrna_cap <- RenameIdents(scrna_cap, "A" = "H_Pre", 
                          "B" = "H_Post", 
                          "C" = "L_Pre", 
                          "D" = "L_Post")
scrna_cap$group_anno <- Idents(scrna_cap)
scrna_cap$Cell_type_fine_harmony <- droplevels(scrna_cap$Cell_type_fine_harmony)
scrna_cap$Cell_type_fine_harmony <- factor(scrna_cap$Cell_type_fine_harmony, levels = c("Tip-like", "Transition_FLT1hi", "Transition_FLT1lo", "Stalk-like"))
Cellratio <- prop.table(table(scrna_cap$Cell_type_fine_harmony, scrna_cap$sample), margin = 2) #计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)

cellper <- reshape2::dcast(Cellratio, Var2~Var1, value.var="Freq") 
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]

#write.csv()
#to_save <- t(cellper)
write.csv(to_save, file='matPC_cells_prop_260803.csv', quote = F)

meta <- scrna_cap@meta.data
colnames(meta) 
meta <- meta[,c(32,4)] # group and sample information 
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
#cellper$group <- NA
meta <- unique(meta)
cellper <- merge(cellper, meta, by = "sample")
colnames(cellper)
colnames(cellper)[6] <- "group" # Group anno to group 
cellper <- as.data.frame(cellper)

pplist =list()
seuratObj_groups = unique(levels(scrna_cap$Cell_type_fine_harmony))

mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
#fill_color <- zzm60colors
#names(fill_color) <- unique(levels(seuratObj$cluster))

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
  #p_position2 = min(max(cellper_$percent)*1.2) * 0.8
  
  print(group_)
  print(cellper_$median)
  mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + 
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
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) + 
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
                          #paired = T, 
                          #test = "t.test",
                          test.args = list(paired = TRUE, 
                                           alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  #my_comparisons1 <-list(c("L_Pre","L_Post"))
  pp1 = pp1 + geom_signif(comparisons = list(c("L_Pre","L_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))}, 
                          y_position = p_position, 
                          #paired = T, 
                          #test = "t.test",
                          test.args = list(paired = TRUE, 
                                           alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("H_Pre","L_Pre")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position1, 
                          #paired = T, 
                          #test = "t.test", 
                          test.args   = list(#paired = TRUE, 
                            alternative = "two.sided"),
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
                 pplist[[4]],pplist[[4]],pplist[[4]],pplist[[4]],
                 #pplist[[5]],
                 #pplist[[6]],
                 #pplist[[7]],
                 #pplist[[8]],
                 #pplist[[9]],
                 ncol = 4 #, nrow = 2
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_Endo_subclus.pdf"), 
    width = 6, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches
print(pps, newpage = FALSE)
dev.off()

## UMAP of CapEC subsets (related to Figure 4F) -----
colors <- c("#9A9A9A", "#6B4E3D", "#B08D57","#7A1E2D") 
mycol <- colors
names(mycol) <- levels(scrna_cap$Cell_type_fine_harmony)

p1 <- DimPlot(scrna_cap, reduction = "umap", group.by = c("Cell_type_fine_harmony"), 
              label = T, raster = T, pt.size = 8,raster.dpi = c(2048, 2048)) + 
  scale_color_manual(values = colors) + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        axis.title = element_blank(), 
        plot.title = element_blank()) + 
  coord_fixed(ratio = 1) + 
  NoLegend()
p1

pdf(paste0("CapEC_cells", "_UMAP_v251231.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

#### Dot plot of markers and Expression across groups (Figure 4I) ------ 
scrna_cap <- NormalizeData(scrna_cap) %>% FindVariableFeatures() %>% ScaleData()
#mar <- FindAllMarkers(scrna_cap, group.by = "Cell_type_fine_harmony", only.pos = T)
#mar$gene[mar$cluster == "Tip-like"][1:200]
#mar$gene[mar$cluster == "Transition_FLT1hi"][1:200]
#mar$gene[mar$cluster == "Stalk-like"][1:200]
Features <- c("PGF", "PXDN", "CD99", "PDGFB", "APLN", "CD276", "VSIR", "FLT1", "KDR", "ICAM1", "DNASE1L3", 
              "SELE", "SELP")

Idents(scrna_cap) <- "Cell_type_fine_harmony"
p <- DotPlot(scrna_cap, Features) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1 <- ggplot(p$data, aes(x = features.plot, y = id)) +   
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +   
  #facet_grid(facets = ~feature.groups,  switch = "x", scales = "free_x", space = "free_x") +    
  scale_radius(breaks = c(25, 50, 75, 100), range = c(0,6)) +   
  theme_classic() + 
  coord_equal() + 
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
pdf(file = paste0(time, "_", "EC_sig_genes.pdf"), 
    width = 9, # The width of the plot in inches
    height = 3) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()

