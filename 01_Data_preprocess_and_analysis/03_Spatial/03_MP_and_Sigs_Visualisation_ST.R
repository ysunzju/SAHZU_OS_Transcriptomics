library(here) # project-root-relative paths; run scripts from repository root

##################### ST data signature visualisation (Related to Figures 3J, Figure S3B) ##################### 
rm(list = ls())
# disabled, run from repo root: setwd(here("data/New_ST_251011/Output"))
library(Seurat)

ST_comb <- readRDS(here("data/TUmorST_comb_with_hallmarkscores_and_meta.rds"))
unique(ST_comb$orig.ident) #  "SAHST002" "SAHST005"  "SAHST003" "SAHST004" "SAHST001"

# Seurat pipeline for ST ----- 
DefaultAssay(ST_comb) <- "Spatial"
Idents(ST_comb) <- "orig.ident"

## SAHST002 
SAHST002 <- subset(ST_comb, idents = "SAHST002")
SAHST002 <- SCTransform(SAHST002, assay = "Spatial", verbose = FALSE)

## SAHST005 
SAHST005 <- subset(ST_comb, idents = "SAHST005")
SAHST005 <- SCTransform(SAHST005, assay = "Spatial", verbose = FALSE)

## SAHST003 
SAHST003 <- subset(ST_comb, idents = "SAHST003")
SAHST003 <- SCTransform(SAHST003, assay = "Spatial", verbose = FALSE)

## SAHST004 
SAHST004 <- subset(ST_comb, idents = "SAHST004")
SAHST004 <- SCTransform(SAHST004, assay = "Spatial", verbose = FALSE)

## SAHST001 
SAHST001 <- subset(ST_comb, idents = "SAHST001")
SAHST001 <- SCTransform(SAHST001, assay = "Spatial", verbose = FALSE)

## Merge 
ST_merge <- merge(x = SAHST002, y = c(SAHST005, SAHST003, SAHST004, SAHST001))

## UMAPing 
DefaultAssay(ST_merge) <- "SCT"
VariableFeatures(ST_merge) <- c(VariableFeatures(SAHST002), VariableFeatures(SAHST005), VariableFeatures(SAHST003), 
                                VariableFeatures(SAHST004), VariableFeatures(SAHST001))
ST_merge <- RunPCA(ST_merge, verbose = FALSE)
ST_merge <- FindNeighbors(ST_merge, dims = 1:50)
ST_merge <- FindClusters(ST_merge, verbose = FALSE)
ST_merge <- RunUMAP(ST_merge, dims = 1:50)

DimPlot(ST_merge, reduction = "umap", group.by = c("ident", "orig.ident"))
saveRDS(ST_merge, "ST_merge_251015.rds")

# Scoring of MPs ----- 
MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
list_for_sigs <- MP_list[c(2, 3, 4)]

ST_merge <- AddModuleScore(ST_merge, list_for_sigs, name = names(list_for_sigs))
names(ST_merge@meta.data)[c(79, 80, 81)] <- c("S_phase", "Os_min", "ECM")
SpatialFeaturePlot(ST_merge, c("Os_min"), crop = T)
SpatialFeaturePlot(ST_merge, c("ECM"), crop = T)
SpatialFeaturePlot(ST_merge, c("S_phase"), crop = T, image.alpha = 0)

# Spatial feature plot (adapted) ----- 
## For F3 - demo of ECM, S phase, and osteo features 
## SAHST005 ----- 
Idents(ST_merge) <- "orig.ident" 
SAHST005_plot <- subset(ST_merge, idents = "SAHST005") 
SAHST005_plot@reductions$spatial = SAHST005_plot@reductions$umap 
SAHST005_plot@reductions$spatial@key = 'spatial_' 
SAHST005_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST005_plot@images$image@coordinates[,c(3,2)]) 
SAHST005_plot@reductions$spatial@cell.embeddings[,2] = -SAHST005_plot@reductions$spatial@cell.embeddings[,2] 
colnames(SAHST005_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2') 

### Plot
color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
data_mat <- data.frame(SAHST005_plot$ECM, SAHST005_plot$S_phase, SAHST005_plot$Os_min ,SAHST005_plot@reductions$spatial@cell.embeddings)
#write.csv(data_mat, file = "SAHST005_MP_sigs.csv")

# ECM 
color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST005_plot.ECM)) +
  geom_point(shape = 16) +
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

# Osteo 
color_scale_val <- range(ST_merge$Os_min)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST005_plot.Os_min)) +
  geom_point(shape = 16) +
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

pdf('SAHST005_Osteo.pdf', height = 3, width = 3)
print(p1)
dev.off()

# S phase 
color_scale_val <- range(ST_merge$S_phase)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST005_plot.S_phase)) +
  geom_point(shape = 16) +
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

pdf('SAHST005_S_phase.pdf', height = 3, width = 3)
print(p1)
dev.off()


## SAHST001 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST001_plot <- subset(ST_merge, idents = "SAHST001")
SAHST001_plot@reductions$spatial = SAHST001_plot@reductions$umap
SAHST001_plot@reductions$spatial@key = 'spatial_'
SAHST001_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST001_plot@images$image@coordinates[,c(3,2)])
SAHST001_plot@reductions$spatial@cell.embeddings[,2] = -SAHST001_plot@reductions$spatial@cell.embeddings[,2]
colnames(SAHST001_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

### Plot
color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
data_mat <- data.frame(SAHST001_plot$ECM, SAHST001_plot$S_phase, SAHST001_plot$Os_min ,SAHST001_plot@reductions$spatial@cell.embeddings)
#write.csv(data_mat, file = "SAHST001_MP_sigs.csv")

# ECM 
color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST001_plot.ECM)) +
  geom_point(shape = 16) +
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

pdf('SAHST001_ECM.pdf', height = 3, width = 3)
print(p1)
dev.off()

# Osteo 
color_scale_val <- range(ST_merge$Os_min)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST001_plot.Os_min)) +
  geom_point(shape = 16) +
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

pdf('SAHST001_Osteo.pdf', height = 3, width = 3)
print(p1)
dev.off()

# S phase 
color_scale_val <- range(ST_merge$S_phase)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST001_plot.S_phase)) +
  geom_point(shape = 16) +
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

pdf('SAHST001_S_phase.pdf', height = 3, width = 3)
print(p1)
dev.off()


## SAHST002 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST002_plot <- subset(ST_merge, idents = "SAHST002")
SAHST002_plot@reductions$spatial = SAHST002_plot@reductions$umap
SAHST002_plot@reductions$spatial@key = 'spatial_'
SAHST002_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST002_plot@images$image@coordinates[,c(3,2)])
SAHST002_plot@reductions$spatial@cell.embeddings[,2] = -SAHST002_plot@reductions$spatial@cell.embeddings[,2]
colnames(SAHST002_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

### Plot
color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
data_mat <- data.frame(SAHST002_plot$ECM, SAHST002_plot$S_phase, SAHST002_plot$Os_min ,SAHST002_plot@reductions$spatial@cell.embeddings)

# ECM 
color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST002_plot.ECM)) +
  geom_point(shape = 16) +
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

pdf('SAHST002_ECM.pdf', height = 3, width = 3)
print(p1)
dev.off()

# Osteo 
color_scale_val <- range(ST_merge$Os_min)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST002_plot.Os_min)) +
  geom_point(shape = 16) +
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

pdf('SAHST002_Osteo.pdf', height = 3, width = 3)
print(p1)
dev.off()

# S phase 
color_scale_val <- range(ST_merge$S_phase)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST002_plot.S_phase)) +
  geom_point(shape = 16) +
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

pdf('SAHST002_S_phase.pdf', height = 3, width = 3)
print(p1)
dev.off()

## SAHST003 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST003_plot <- subset(ST_merge, idents = "SAHST003")
SAHST003_plot@reductions$spatial = SAHST003_plot@reductions$umap
SAHST003_plot@reductions$spatial@key = 'spatial_'
SAHST003_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST003_plot@images$image@coordinates[,c(3,2)])
SAHST003_plot@reductions$spatial@cell.embeddings[,2] = -SAHST003_plot@reductions$spatial@cell.embeddings[,2]
colnames(SAHST003_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

### Plot
color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
data_mat <- data.frame(SAHST003_plot$ECM, SAHST003_plot$S_phase, SAHST003_plot$Os_min ,SAHST003_plot@reductions$spatial@cell.embeddings)

# ECM 
color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST003_plot.ECM)) +
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

pdf('SAHST003_ECM.pdf', height = 3, width = 3)
print(p1)
dev.off()

# Osteo 
color_scale_val <- range(ST_merge$Os_min)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST003_plot.Os_min)) +
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

pdf('SAHST003_Osteo.pdf', height = 3, width = 3)
print(p1)
dev.off()

# S phase 
color_scale_val <- range(ST_merge$S_phase)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST003_plot.S_phase)) +
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

pdf('SAHST003_S_phase.pdf', height = 3, width = 3)
print(p1)
dev.off()

## SAHST004 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST004_plot <- subset(ST_merge, idents = "SAHST004")
SAHST004_plot@reductions$spatial = SAHST004_plot@reductions$umap
SAHST004_plot@reductions$spatial@key = 'spatial_'
SAHST004_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST004_plot@images$image@coordinates[,c(3,2)])
SAHST004_plot@reductions$spatial@cell.embeddings[,2] = -SAHST004_plot@reductions$spatial@cell.embeddings[,2]
colnames(SAHST004_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

### Plot
color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
data_mat <- data.frame(SAHST004_plot$ECM, SAHST004_plot$S_phase, SAHST004_plot$Os_min ,SAHST004_plot@reductions$spatial@cell.embeddings)

# ECM 
color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot.ECM)) +
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

pdf('SAHST004_ECM.pdf', height = 3, width = 3)
print(p1)
dev.off()

# Osteo 
color_scale_val <- range(ST_merge$Os_min)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot.Os_min)) +
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

pdf('SAHST004_Osteo.pdf', height = 3, width = 3)
print(p1)
dev.off()

# S phase 
color_scale_val <- range(ST_merge$S_phase)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot.S_phase)) +
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

pdf('SAHST004_S_phase.pdf', height = 3, width = 3)
print(p1)
dev.off()

# Visualise (violin) of key MPs -----
library(SCP)
library(scop)
library(ggplot2)
Idents(ST_merge) <- "orig.ident"
ST_for_plot <- subset(ST_merge, idents = c("SAHST005", "SAHST001"))

my_comparisons <- list(c("SAHST005", "SAHST001"))
FeatureStatPlot(srt = ST_merge, group.by = "orig.ident",  
                stat.by = "CD1C", 
                add_box = TRUE,  
                #palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")), 
                #bg_palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")),
                comparisons = my_comparisons,  
                add_trend = T, 
                bg.by = 'orig.ident') + theme(legend.position = "none") + 
  theme(legend.position = "none", 
        #axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.y = element_blank())

pps <- plot_grid(pplist[[1]], 
                 pplist[[2]],
                 pplist[[3]], 
                 pplist[[4]],
                 pplist[[5]],
                 pplist[[6]],
                 pplist[[7]],
                 pplist[[8]],
                 pplist[[9]])
pps



# Visualisation of colocalisation of signatures ----
## Signatures ---- 
# General inspection
ST_merge <- readRDS("ST_merge_251015.rds")

genelist_down_manual <- c("FCER1A", "CD1C", "CLEC10A")
#gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "cDC2_CD83hi")
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "cDC2")
#gene.set_v2 <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v1 <- rbind(gene.set_v1, gene.set_v2)
#rm(gene.set_v2)
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

#ST_merge <- readRDS("ST_merge_251015.rds")
ST_merge <- AddModuleScore(ST_merge, geneSet, name = names(geneSet))

names(ST_merge@meta.data)[c(79)] <- c("cDC2")
SpatialFeaturePlot(ST_merge, c("cDC2"), crop = T)

### cDC2 ----- 
##### Plotting cDC2 (SAHST004 mainly)-----
## Starting with SAHST004 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST004_plot <- subset(ST_merge, idents = "SAHST004")
SAHST004_plot@reductions$spatial = SAHST004_plot@reductions$umap
SAHST004_plot@reductions$spatial@key = 'spatial_'
SAHST004_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST004_plot@images$image@coordinates[,c(3,2)])
SAHST004_plot@reductions$spatial@cell.embeddings[,2] = -SAHST004_plot@reductions$spatial@cell.embeddings[,2]
colnames(SAHST004_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

### Plot
color_scale_val <- range(ST_merge$cDC2#, ST_merge$S_phase, ST_merge$Os_min
                         )
data_mat <- data.frame(SAHST004_plot$cDC2, #SAHST004_plot$S_phase, SAHST004_plot$Os_min ,
                       SAHST004_plot@reductions$spatial@cell.embeddings)

# cDC2
color_scale_val <- range(ST_merge$cDC2)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot.cDC2)) +
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

pdf('SAHST004_cDC2.pdf', height = 3, width = 3)
print(p1)
dev.off()

genelist_down_manual <- c("CD3D", "CD3E", "CD2")
#gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "cDC2_CD83hi")
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "Th")
#gene.set_v2 <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v1 <- rbind(gene.set_v1, gene.set_v2)
#rm(gene.set_v2)
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

#ST_merge <- readRDS("ST_merge_251015.rds")
ST_merge <- AddModuleScore(ST_merge, geneSet, name = names(geneSet))
#names(ST_merge@meta.data)[c(79)] <- c("cDC2_CD83hi")
names(ST_merge@meta.data)[c(80)] <- c("T_cell_sig")

##### Plotting T (SAHST004 mainly)-----
###### Starting with SAHST004 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST004_plot <- subset(ST_merge, idents = "SAHST004")
SAHST004_plot@reductions$spatial = SAHST004_plot@reductions$umap
SAHST004_plot@reductions$spatial@key = 'spatial_'
SAHST004_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST004_plot@images$image@coordinates[,c(3,2)])
SAHST004_plot@reductions$spatial@cell.embeddings[,2] = -SAHST004_plot@reductions$spatial@cell.embeddings[,2]
colnames(SAHST004_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

### Plot
color_scale_val <- range(ST_merge$T_cell_sig#, ST_merge$S_phase, ST_merge$Os_min
)
data_mat <- data.frame(SAHST004_plot$T_cell_sig, #SAHST004_plot$S_phase, SAHST004_plot$Os_min ,
                       SAHST004_plot@reductions$spatial@cell.embeddings)

color_scale_val <- range(ST_merge$T_cell_sig)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot.T_cell_sig)) +
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

pdf('SAHST004_T.pdf', height = 3, width = 3)
print(p1)
dev.off()

## Plotting expressions on SAHST004 -----
Idents(ST_merge) <- "orig.ident"
SAHST004_plot <- subset(ST_merge, idents = "SAHST004")
SAHST004_plot@reductions$spatial = SAHST004_plot@reductions$umap
SAHST004_plot@reductions$spatial@key = 'spatial_'
SAHST004_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST004_plot@images$image@coordinates[,c(3,2)])
SAHST004_plot@reductions$spatial@cell.embeddings[,2] = -SAHST004_plot@reductions$spatial@cell.embeddings[,2]
colnames(SAHST004_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

SAHST004_data = GetAssayData(SAHST004_plot, assay="SCT", layer='data') %>% as.data.frame()
genes <- c("IL7R", "FOXP3", "GZMB")

### CD8 ----- 
CD8_sigs <- readxl::read_xlsx(here("data/Re-clustering_fine_cell_types_21_05_25/T_and_NK/CD8_Sigs_NM.xlsx"))
colnames(CD8_sigs)
CD8_sigs_list <- list()
CD8_sigs_list
for (i in colnames(CD8_sigs)) {
  j <- CD8_sigs[i] %>% na.omit() %>% as.vector()
  CD8_sigs_list[i] <- j
}
#CD8_sigs_list[["Oxidative phosphorylation"]] <- gsub("\\.", "-", CD8_sigs_list)
names(CD8_sigs_list)[2] <- "Activation/Effector function"
names(CD8_sigs_list)[4] <- "TCR signaling"
names(CD8_sigs_list)[14] <- "IFN response"
names(CD8_sigs_list)[15] <- "OXPHOS"
names(CD8_sigs_list)[17] <- "FA/lipid metabolism"
#CD8 <- AddModuleScore(CD8, features = CD8_sigs_list, name = names(CD8_sigs_list))
#names(CD8@meta.data)[57:75] <- names(CD8_sigs_list)

### CD4 ----- 
CD4_sigs <- readxl::read_xlsx(here("data/Re-clustering_fine_cell_types_21_05_25/T_and_NK/CD4_sigs_NM.xlsx"))
colnames(CD4_sigs)
CD4_sigs_list <- list()
CD4_sigs_list
for (i in colnames(CD4_sigs)) {
  j <- CD4_sigs[i] %>% na.omit() %>% as.vector()
  CD4_sigs_list[i] <- j
}
names(CD4_sigs_list)[2] <- "Activation/Effector function"
names(CD4_sigs_list)[4] <- "TCR signaling"
names(CD4_sigs_list)[10] <- "IFN response"
names(CD4_sigs_list)[15] <- "FA/lipid metabolism"
#CD4 <- AddModuleScore(CD4, features = CD4_sigs_list, name = names(CD4_sigs_list))
#names(CD4@meta.data)[57:73] <- names(CD4_sigs_list)

CD4_to_add_to_CD8 <- CD4_sigs_list[names(CD4_sigs_list)[!names(CD4_sigs_list) %in% intersect(names(CD8_sigs_list), names(CD4_sigs_list))]]
CD8_union <- union(CD8_sigs_list[], CD4_to_add_to_CD8[])
names(CD8_union) <- union(names(CD8_sigs_list), names(CD4_to_add_to_CD8))
CD8_to_add_to_CD4 <- CD8_sigs_list[names(CD8_sigs_list)[!names(CD8_sigs_list) %in% intersect(names(CD4_sigs_list), names(CD8_sigs_list))]]
CD4_union <- union(CD4_sigs_list[], CD8_to_add_to_CD4[])
names(CD4_union) <- union(names(CD4_sigs_list), names(CD8_to_add_to_CD4))

ST_merge <- AddModuleScore(ST_merge, features = CD4_union, name = names(CD4_union))
names(ST_merge@meta.data)[81:101] <- names(CD4_union)

#CD8 <- AddModuleScore(CD8, features = CD8_union, name = names(CD8_union))
#names(CD8@meta.data)[57:77] <- names(CD8_union)

### TAM (or general M2) ------
seuratObj <- readRDS(here("data/Correlation_and_crosstalk_analysis_18_07_25/All_cells_for_ecoTyper_cell_type_adjusted_250911.rds"))
Idents(seuratObj) <- "Cell_type_fine_harmony"
markers <- FindMarkers(seuratObj, ident.1 = "Macro_APOE+", only.pos = T, logfc.threshold = 1#, min.diff.pct = 0.25
) # This is the preferred method 
rm(seuratObj); gc()
M2_sig <- list(M2_sig = rownames(markers)[order(markers$p_val_adj, decreasing = F)][1:50])

ST_merge <- AddModuleScore(ST_merge, features = M2_sig, name = names(M2_sig))
names(ST_merge@meta.data)[102] <- names(M2_sig)

### Plot for SAHST004 for crosstalk figure ---- 
Idents(ST_merge) <- "orig.ident"
SAHST004_plot <- subset(ST_merge, idents = "SAHST004")
SAHST004_plot@reductions$spatial = SAHST004_plot@reductions$umap
SAHST004_plot@reductions$spatial@key = 'spatial_'
SAHST004_plot@reductions$spatial@cell.embeddings = as.matrix(SAHST004_plot@images$image@coordinates[,c(3,2)])
SAHST004_plot@reductions$spatial@cell.embeddings[,2] = -SAHST004_plot@reductions$spatial@cell.embeddings[,2]
colnames(SAHST004_plot@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')

# Activation/Effector function 
color_scale_val <- range(SAHST004_plot$`Activation/Effector function`#, ST_merge$S_phase, ST_merge$Os_min
)
data_mat <- data.frame(SAHST004_plot$`Activation/Effector function`, #SAHST004_plot$S_phase, SAHST004_plot$Os_min ,
                       SAHST004_plot@reductions$spatial@cell.embeddings)

colnames(data_mat)[1]
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot..Activation.Effector.function.)) +
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

pdf('SAHST004_activation_effector.pdf', height = 3, width = 3)
print(p1)
dev.off()

# M2 
color_scale_val <- range(SAHST004_plot$`M2_sig`#, ST_merge$S_phase, ST_merge$Os_min
)
data_mat <- data.frame(SAHST004_plot$`M2_sig`, #SAHST004_plot$S_phase, SAHST004_plot$Os_min ,
                       SAHST004_plot@reductions$spatial@cell.embeddings)

colnames(data_mat)[1]
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot.M2_sig)) +
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

pdf('SAHST004_M2.pdf', height = 3, width = 3)
print(p1)
dev.off()

# Treg 
color_scale_val <- range(SAHST004_plot$`Treg signature`#, ST_merge$S_phase, ST_merge$Os_min
)
data_mat <- data.frame(SAHST004_plot$`Treg signature`, #SAHST004_plot$S_phase, SAHST004_plot$Os_min ,
                       SAHST004_plot@reductions$spatial@cell.embeddings)

colnames(data_mat)[1]
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot..Treg.signature.)) +
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

pdf('SAHST004_Treg.pdf', height = 3, width = 3)
print(p1)
dev.off()

# Cytotoxicity 
color_scale_val <- range(SAHST004_plot$Cytotoxicity#, ST_merge$S_phase, ST_merge$Os_min
)
data_mat <- data.frame(SAHST004_plot$Cytotoxicity, #SAHST004_plot$S_phase, SAHST004_plot$Os_min ,
                       SAHST004_plot@reductions$spatial@cell.embeddings)

colnames(data_mat)[1]
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST004_plot.Cytotoxicity)) +
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

pdf('SAHST004_Cytotoxicity.pdf', height = 3, width = 3)
print(p1)
dev.off()