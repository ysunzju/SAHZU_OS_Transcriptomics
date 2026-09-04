library(here) # project-root-relative paths; run scripts from repository root
# Updated deconv methods benchmarking 
rm(list = ls())
library(Seurat)
#devtools::install_github('YingMa0107/CARD')
library(CARD)
#BiocManager::install("TOAST")
#devtools::install_github('xuranw/MuSiC')
#install.packages("ggplot2")
library(MuSiC)
library(tidyverse)
library(ggplot2)

# disabled, run from repo root: setwd(here("data/New_ST_251011/Output"))
ST_merge <- readRDS("ST_merge_251015.rds")

# Reference ------ 
#seuratObj <- readRDS(here("data/Correlation_and_crosstalk_analysis_18_07_25/All_cells_for_ecoTyper_cell_type_adjusted_250911.rds"))
seuratObj <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))
seuratObj <- NormalizeData(seuratObj) %>% FindVariableFeatures() %>% ScaleData()
names(seuratObj@meta.data)

Idents(seuratObj) <- "cluster" # Change to major cell type 
DimPlot(seuratObj,
        label = T)
seuratObj$cluster <- Idents(seuratObj)

# CARD for major cell types ----- 
## SAHST001 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST001 <- subset(ST_merge, idents = "SAHST001")

# Count matrix
spatial_count <-  SAHST001@assays$Spatial@counts
spatial_count[1:4,1:4]

# Spatial location
spatial_loca <- SAHST001@images$image@coordinates
spatial_location <- spatial_loca[,2:3]
colnames(spatial_location) <- c("x","y")
spatial_location[1:4,]

# scRNA count 
sc_count[1:4,1:4]

# scRNA meta data 
sc_meta <- readRDS("sc_meta_251115_all_cells.rds")

# Object creation 
CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cluster_level_adjusted",
  ct.select = unique(sc_meta$cluster_level_adjusted),
  sample.varname = "sample",
  minCountGene = 100,
  minCountSpot = 5) 

### Deconv ----- 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
x <- CARD_obj@spatial_location$x 
y <- CARD_obj@spatial_location$y
new_x <- y 
new_y <- -x 

CARD_obj@spatial_location$x <- new_x
CARD_obj@spatial_location$y <- new_y

colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

## select the cell type that we are interested
ct.visualize = rev(levels(sc_meta$cluster_level_adjusted))

### Refined spatial map ----- 
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)

## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)

library(ggplot2)
p5 <- ggplot(location_imputation, 
             aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p5)

range(location_imputation$x) #0.85 96.05
range(location_imputation$y) #-68.0 -13.6

p6 <- CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
             "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
             "#d73027", "#a50026"),    
  pointSize = 0.5, 
  NumCols = 16)                                  
print(p6)

pp6 <- p6 + coord_fixed(ratio = 95.2/54.4) + NoLegend()
pp6

pdf("SAHST001_All_cell_disect.pdf", height = 3, width = 16)
print(pp6)
dev.off()

rm(SAHST004)

## SAHST002 ----- 
Idents(ST_merge) <- "orig.ident" 
SAHST002 <- subset(ST_merge, idents = "SAHST002") 

# Count matrix
spatial_count <-  SAHST002@assays$Spatial@counts
spatial_count[1:4,1:4]

# Spatial location
spatial_loca <- SAHST002@images$image@coordinates
spatial_location <- spatial_loca[,2:3]
colnames(spatial_location) <- c("x","y")
spatial_location[1:4,]

### scRNA count matrix 
sc_count <- seuratObj@assays$RNA@counts
sc_count[1:4,1:4]

# sc_meta 
sc_meta <- seuratObj@meta.data %>% 
  rownames_to_column("cellID") %>%
  dplyr::select(cellID,sample,cluster_level_adjusted) %>% 
  mutate(CB = cellID) %>% 
  column_to_rownames("CB")
head(sc_meta)

saveRDS(sc_meta, "sc_meta_251115_all_cells.rds")


# Object creation 
CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cluster_level_adjusted",
  ct.select = unique(sc_meta$cluster_level_adjusted),
  sample.varname = "sample",
  minCountGene = 100,
  minCountSpot = 5) 

### Deconv ----- 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

colors <- c('#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',
            '#CC6677', '#882255', '#AA4499', '#662288', '#BB5566',
            '#5D3A9B', '#997700', '#DDDDDD', '#004488',
            '#6699CC', '#66BBAA')  
names(colors) <- levels(seuratObj$cluster)

range(CARD_obj@spatial_location$x)
#  3 56
range(CARD_obj@spatial_location$y)
#  38 103

x <- CARD_obj@spatial_location$x 
y <- CARD_obj@spatial_location$y
new_x <- y
new_y <- -x

CARD_obj@spatial_location$x <- new_x
CARD_obj@spatial_location$y <- new_y

p1 <- CARD.visualize.pie(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location, 
  colors = colors) ### You can choose radius = NULL or your own radius number
print(p1)

p1 + coord_fixed(1/0.56) #+ scale_y_reverse()

## select the cell type that we are interested
ct.visualize = rev(levels(sc_meta$cluster_level_adjusted))

## visualize the spatial distribution of the cell type proportion
## using ggplot 2 
data_mat <- data.frame(CARD_obj@Proportion_CARD, CARD_obj@spatial_location)

# OB for instance 
p1 <- ggplot(data_mat, aes(x = x, y = y, colour = Osteoblasts)) +
  geom_point(shape = 16) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026")) + 
  coord_fixed(ratio = 1/0.56) + 
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf("SAHST002_OB.pdf", height = 5, width = 5)
print(p1)
dev.off()

## Using Seurat in-built functions (Related to Figure 1E) -----
SAHST002 <- AddMetaData(SAHST002, data_mat)
#OB 
SpatialFeaturePlot(SAHST002, "Osteoblasts", crop = F) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026"))

#Endo
SpatialFeaturePlot(SAHST002, "ECs", crop = F) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026"))

#TNK 8x8 
SpatialFeaturePlot(SAHST002, "TandNK", crop = F) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026"))

#Myeloid
SpatialFeaturePlot(SAHST002, "Myeloid", crop = F) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026"))

### Refined spatial map ----- 
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)

## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)

p5 <- ggplot(location_imputation, 
             aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p5)

range(location_imputation$x) #39.0 102.7 = 63.7
range(location_imputation$y) #-55.0  -4.3 =. 50.7

p6 <- CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
              "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
              "#d73027", "#a50026"),    
  NumCols = 16,
  pointSize = 0.2
)                                  
print(p6)

pp6 <- p6 + coord_fixed(63.7/50.7) + NoLegend()
pp6

pdf("SAHST002_All_cell_disect.pdf", height = 3, width = 16)
print(pp6)
dev.off()

rm(SAHST002)

## SAHST003 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST003 <- subset(ST_merge, idents = "SAHST003")

# Count matrix
spatial_count <- SAHST003@assays$Spatial@counts
spatial_count[1:4,1:4]

# Spatial location
spatial_loca <- SAHST003@images$image@coordinates
spatial_location <- spatial_loca[,2:3]
colnames(spatial_location) <- c("x","y")
spatial_location[1:4,]

# scRNA count 
sc_count[1:4,1:4]

# scRNA meta data
sc_meta <- readRDS("sc_meta_251115_all_cells.rds")

# Object creation 
CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cluster_level_adjusted",
  ct.select = unique(sc_meta$cluster_level_adjusted),
  sample.varname = "sample",
  minCountGene = 100,
  minCountSpot = 5) 

### Deconv ----- 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
x <- CARD_obj@spatial_location$x 
y <- CARD_obj@spatial_location$y
new_x <- y
new_y <- -x

CARD_obj@spatial_location$x <- new_x
CARD_obj@spatial_location$y <- new_y

## select the cell type that we are interested
ct.visualize = rev(levels(sc_meta$cluster_level_adjusted))

### Refined spatial map ----- 
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)

## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)

p5 <- ggplot(location_imputation, 
             aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p5)

range(location_imputation$x) 
range(location_imputation$y) 

p6 <- CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
             "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
             "#d73027", "#a50026"),    
  NumCols = 16,
  pointSize = 0.2
)                                  
print(p6)

pp6 <- p6 + coord_fixed(84/58) + NoLegend()
pp6

pdf("SAHST003_All_cell_disect.pdf", height = 3, width = 16)
print(pp6)
dev.off()

rm(SAHST003)

## SAHST004 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST004 <- subset(ST_merge, idents = "SAHST004")

# Count matrix
spatial_count <-  SAHST004@assays$Spatial@counts
spatial_count[1:4,1:4]

# Spatial location
spatial_loca <- SAHST004@images$image@coordinates
spatial_location <- spatial_loca[,2:3]
colnames(spatial_location) <- c("x","y")
spatial_location[1:4,]

### scRNA count matrix 
sc_count[1:4,1:4]

# scRNA meta data
sc_meta <- readRDS("sc_meta_251115_all_cells.rds")

# Object creation 
CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cluster_level_adjusted",
  ct.select = unique(sc_meta$cluster_level_adjusted),
  sample.varname = "sample",
  minCountGene = 100,
  minCountSpot = 5) 

### Deconv ----- 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
x <- CARD_obj@spatial_location$x 
y <- CARD_obj@spatial_location$y
new_x <- -x 
new_y <- -y 

CARD_obj@spatial_location$x <- new_x
CARD_obj@spatial_location$y <- new_y

colors <- c('#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77',
                     '#CC6677', '#882255', '#AA4499', '#662288', '#BB5566',
                     '#5D3A9B', '#997700', '#DDDDDD', '#004488',
                     '#6699CC', '#66BBAA') 
names(colors) <- levels(seuratObj$cluster)
#[1] "BCells"       "Chondrocytes" "ECs"          "Erythrocytes" "Fibroblasts"  "HSCs"         "Myeloid"     
#[8] "MSCs"         "MastCells"    "MuralCells"   "Osteoblasts"  "Osteoclasts"  "PlasmaCells"  "Pro_B"       
#[15] "TandNK"       "pDCs" 

p1 <- CARD.visualize.pie(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location, 
  colors = colors#, 
  #radius = 0.5
  ) ### You can choose radius = NULL or your own radius number
print(p1)
pp1 <- p1 + coord_fixed(0.6062992126) #+ scale_y_reverse()
pp1

pdf("SAHST004_deconv.pdf", width = 5, height = 5)
print(pp1)
dev.off()

## select the cell type that we are interested
ct.visualize = rev(levels(sc_meta$cluster_level_adjusted))
## visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = ct.visualize,                 
  colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
             "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
             "#d73027", "#a50026"),
  NumCols = 8,                                 
  pointSize = 0.2)                          
print(p2)
p3 <- p2 + 
  coord_fixed(0.7) 
p3

pdf("SAHST004_decon_proportions.pdf", width = 15, height = 4) 
print(p3)
dev.off()

## Manual Vis of proportions ----- 
#color_scale_val <- range(ST_merge$ECM, ST_merge$S_phase, ST_merge$Os_min)
CARD_obj@Proportion_CARD
data_mat <- data.frame(CARD_obj@Proportion_CARD, CARD_obj@spatial_location)
#write.csv(data_mat, file = "SAHST003_MP_sigs.csv")
range(data_mat$x)
#[1]  0 77
range(data_mat$y)
#[1]   0 127

#1.6493506494
#0.6062992126
# OB 
#color_scale_val <- range(ST_merge$ECM)
p1 <- ggplot(data_mat, aes(x = x, y = y, colour = Osteoblasts)) +
  geom_point(shape = 16) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026")) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 0.6062992126) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf("SAHST004_OB.pdf", height = 5, width = 5)
print(p1)
dev.off()

## EC
p1 <- ggplot(data_mat, aes(x = x, y = y, colour = ECs)) +
  geom_point(shape = 16) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026")) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 0.6062992126) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf("SAHST004_EC.pdf", height = 5, width = 5)
print(p1)
dev.off()

### Refined spatial map ----- 
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)

## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)

p5 <- ggplot(location_imputation, 
             aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p5)

p6 <- CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
             "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
             "#d73027", "#a50026"),    
  pointSize = 0.5, 
  NumCols = 16)                                  
print(p6)

pp6 <- p6 + coord_fixed(127/77) + NoLegend()
pp6

pdf("SAHST004_All_cell_disect.pdf", height = 3, width = 16)
print(pp6)
dev.off()

rm(SAHST004)

## SAHST005 ----- 
Idents(ST_merge) <- "orig.ident"
SAHST005 <- subset(ST_merge, idents = "SAHST005")

# Count matrix
spatial_count <-  SAHST005@assays$Spatial@counts
spatial_count[1:4,1:4]

# Spatial location
spatial_loca <- SAHST005@images$image@coordinates
spatial_location <- spatial_loca[,2:3]
colnames(spatial_location) <- c("x","y")
spatial_location[1:4,]

# scRNA meta data 
sc_meta <- readRDS("sc_meta_251115_all_cells.rds")

CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "cluster_level_adjusted",
  ct.select = unique(sc_meta$cluster_level_adjusted),
  sample.varname = "sample",
  minCountGene = 100,
  minCountSpot = 5) 

### Deconv ----- 
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
x <- CARD_obj@spatial_location$x 
y <- CARD_obj@spatial_location$y
new_x <- y 
new_y <- -x 

CARD_obj@spatial_location$x <- new_x
CARD_obj@spatial_location$y <- new_y

colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1 <- CARD.visualize.pie(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location, 
  colors = colors, 
  radius = 0.5)
print(p1)
p1 + coord_fixed(0.7) 

## select the cell type that we are interested
ct.visualize = rev(levels(sc_meta$cluster_level_adjusted))
## visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = ct.visualize,                 
  colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
             "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
             "#d73027", "#a50026"), 
  NumCols = 4,                               
  pointSize = 0.2)                        
print(p2)
p3 <- p2 + 
  coord_fixed(0.5)  
p3

### Refined spatial map ----- 
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)

## Visualize the newly grided spatial locations to see if the shape is correctly detected. If not, the user can provide the row names of the excluded spatial location data into the CARD.imputation function
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)

p5 <- ggplot(location_imputation, 
             aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p5)

range(location_imputation$x) #  0.9 126.9
range(location_imputation$y) #  -64.0 -17.2

p6 <- CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
             "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
             "#d73027", "#a50026"),    
  pointSize = 0.5, 
  NumCols = 16)                                  
print(p6)

pp6 <- p6 + coord_fixed(ratio = 126/46.8) + NoLegend()
pp6

pdf("SAHST005_All_cell_disect.pdf", height = 3, width = 16)
print(pp6)
dev.off()

rm(SAHST005)
