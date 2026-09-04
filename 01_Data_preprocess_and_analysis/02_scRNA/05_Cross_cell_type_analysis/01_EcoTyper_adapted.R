library(here) # project-root-relative paths; run scripts from repository root
##################### EcoTyper discovery and assigment (related to Figure 7) #####################
rm(list = ls())
library(Seurat)
library(tidyverse)
library(reshape2)
library(plyr)
library(Biobase)
library(argparse)
library(matrixTests)
library(HiClimR)
library(ggalluvial)
library(svglite)
library(plotly)
library(Seurat)
library(cluster)
library(viridis)
library(ggplot2)

# disabled, run from repo root: setwd(here("data/Correlation_and_crosstalk_analysis_18_07_25"))
# NOTE: misc.R and heatmaps.R are helper files from the EcoTyper source code
# (https://github.com/digitalcytometry/ecotyper, ecotyper/pipeline/lib/).
# Place copies of both files in the lib/ directory next to this script.
source(here("01_Data_preprocess_and_analysis/02_scRNA/05_Cross_cell_type_analysis/lib/misc.R"))
source(here("01_Data_preprocess_and_analysis/02_scRNA/05_Cross_cell_type_analysis/lib/heatmaps.R"))

# Run all cells ----- 
G <- "All_cells" # name the outputfiles as need be 
#seuratObj <- readRDS("All_cells_for_ecoTyper_cell_type_adjusted_250828.rds")
#seuratObj <- readRDS("All_cells_for_ecoTyper_cell_type_adjusted_250911.rds")
seuratObj <- readRDS(here("data/Correlation_and_crosstalk_analysis_18_07_25/All_cells_for_ecoTyper_cell_type_adjusted_250911.rds"))
seuratObj$Cell_type_for_ecotyper <- Idents(seuratObj)

Idents(seuratObj) <- "Cell_type_fine_harmony"
to_keep <- as.character(unique(Idents(seuratObj))[!unique(Idents(seuratObj)) %in% c("T_Mtgenes_hi", "pDCs", "Mesen_housekeeping", "T_Ki67+", NA)]) # To further delete T_MTgenes_hi 
seuratObj <- subset(seuratObj, idents = to_keep)
Idents(seuratObj) <- "Cell_type_for_ecotyper" 

### MonoMac ----- 
MonoMacro <- subset(seuratObj, idents = "MonoMacro") 
MonoMacro <- MonoMacro@meta.data # c('Sample', 'Sample_Origin', 'Datasets', 'Datatype', 'CellType')
MonoMacro <- MonoMacro[, c("sample", "group", "Cell_type_fine_harmony")] 
colnames(MonoMacro) <- c("Sample", "Datatype", "CellType") 
#MonoMacro$cell_type <- 'MonoMacro'
MonoMacro$cell_type <- 'MonoMacro'

### Th ----- 
Th <- subset(seuratObj, idents = "Th") 
Th <- Th@meta.data # c('Sample', 'Sample_Origin', 'Datasets', 'Datatype', 'CellType')
Th <- Th[, c("sample", "group", "Cell_type_fine_harmony")] 
colnames(Th) <- c("Sample", "Datatype", "CellType") 
Th$cell_type <- 'Th' 

# Tc ----- 
Tc <- subset(seuratObj, idents = "Tc") 
Tc <- Tc@meta.data # c('Sample', 'Sample_Origin', 'Datasets', 'Datatype', 'CellType')
Tc <- Tc[, c("sample", "group", "Cell_type_fine_harmony")] 
colnames(Tc) <- c("Sample", "Datatype", "CellType") 
Tc$cell_type <- 'Tc' 

# NK ----- 
NK <- subset(seuratObj, idents = "NK") 
NK <- NK@meta.data # c('Sample', 'Sample_Origin', 'Datasets', 'Datatype', 'CellType')
NK <- NK[, c("sample", "group", "Cell_type_fine_harmony")] 
colnames(NK) <- c("Sample", "Datatype", "CellType") 
NK$cell_type <- 'NK' 

# cDC2 ------ 
cDC2 <- subset(seuratObj, idents = c("cDC1", "cDC2", "mDCs"#, "pDCs"
)) 
cDC2 <- cDC2@meta.data # c('Sample', 'Sample_Origin', 'Datasets', 'Datatype', 'CellType')
cDC2 <- cDC2[, c("sample", "group", "Cell_type_fine_harmony")] 
colnames(cDC2) <- c("Sample", "Datatype", "CellType") 
#cDC2$CellType <- as.character(cDC2$CellType) 
#cDC2$CellType[cDC2$CellType == "cDC2_CD83lo"] <- "cDC2"
#cDC2$CellType[cDC2$CellType == "cDC2_CD83hi"] <- "cDC2"
cDC2$cell_type <- 'DC' 

# Combined metadata ----- 
meta <- rbind(#ec, 
  MonoMacro, 
  Th, 
  Tc, 
  NK, 
  cDC2#, 
  #Granulocytic, 
  #Mural, 
  #Osteoclasts, 
  #Malignant, 
  #Mesenchymal#, 
  #B
)
meta$CellType <- factor(meta$CellType, levels = names(table(meta$CellType)))
meta$CellType_new <- meta$CellType
meta1 <- meta[, c('Sample', 'Datatype', 'cell_type', 'CellType_new')]

# run -----
Data_in <- meta1 
Data_in$ID <- rownames(Data_in)
#Data_in <- Data_in %>% filter(Sample %in% SampleFi) # %>% filter(Response == "R")

Data_in <- Data_in %>%
  group_by(cell_type) %>%
  mutate(
    InitialState = {
      sub_levels <- unique(CellType_new)
      mapping <- setNames(paste0("IS", sprintf("%02d", 1:length(sub_levels))), sub_levels)
      mapping[as.character(CellType_new)]
    }
  ) %>%
  ungroup()

Data_in$State <- gsub("I", "", Data_in$InitialState)
path_root <- here("data/Correlation_and_crosstalk_analysis_18_07_25/EcoTyper_rerun_250827")

for(CT_NM in unique(Data_in$cell_type)){
  path_out <- file.path(path_root, "/02_EcoPre", CT_NM)
  if (!dir.exists(path_out)) {
    dir.create(path_out, recursive = TRUE)
  }
  
  Data_CTmajor <- Data_in[Data_in$cell_type == CT_NM, ]
  
  ## mapping_to_initial_states
  CTmap_INFO <- unique(Data_CTmajor[c("State", "InitialState", "CellType_new")]) %>% arrange(InitialState)
  # names(CTmap_INFO) <- c("State", "InitialState")
  # CTmap_INFO$State <- CTmap_INFO$InitialState
  CTmap_file <- file.path(path_out, "mapping_to_initial_states.txt")
  write.table(CTmap_INFO, CTmap_file, sep = "\t", row.names = FALSE)
  
  ## initial_state_assignment
  CellCT_INFO <- Data_CTmajor[c("ID", "InitialState", "Sample")]
  names(CellCT_INFO) <- c("ID", "State", "Sample")
  CellCT_file <- file.path(path_out, "initial_state_assignment.txt")
  write.table(CellCT_INFO, CellCT_file, sep = "\t", row.names = FALSE)
}

# Discovery scRNA new -----
path_root <- here("data/Correlation_and_crosstalk_analysis_18_07_25/EcoTyper_rerun_250827")
path_state <- file.path(path_root, "/02_EcoPre")
path_Ecout <- file.path(path_root, "/03_EcoDisc")
if (!dir.exists(path_Ecout)) {
  dir.create(path_Ecout, recursive = TRUE)
}

p_val_cutoff = 0.05
min_states = 3

all_mapping = NULL
all_classes = NULL
all_classes_filt = NULL

CT_all <- unique(Data_in$cell_type)

for(cell_type in CT_all){	
  
  mapping_path = file.path(path_state, cell_type, 'mapping_to_initial_states.txt')
  classes_path = file.path(path_state, cell_type, 'initial_state_assignment.txt')
  
  mapping = read.delim(mapping_path)
  mapping = mapping[,c("State", "InitialState", "CellType_new")] 
  mapping$CellType = cell_type
  all_mapping = rbind(all_mapping, mapping)
  
  classes = read.delim(classes_path)
  #clinical = read_clinical(classes$ID, dataset = dataset)
  #classes$Sample = clinical$Sample
  
  classes = as.data.frame(table(classes$Sample, classes$State))
  colnames(classes) = c("ID", "State", "Freq")
  splits = split(classes, classes$ID)
  classes = do.call(rbind, lapply(splits, function(spl){
    spl$Frac = spl$Freq / sum(spl$Freq)
    # spl$Max = ifelse((!is.na(max(spl$Freq))) & (max(spl$Freq) > 0) & (max(spl$Freq) == spl$Freq), 1, 0) 
    # spl <- spl %>% mutate(Max = ifelse(Freq %in% sort(Freq, decreasing = TRUE)[1:2], 1, 0))
    ####
    if(cell_type %in% c("Malignant")#"Malignant"
    ){
      spl <- spl %>% mutate(Max = ifelse(Freq %in% sort(Freq, decreasing = TRUE)[1:2], 1, 0)) 
      # Median_Frac <- mean(spl$Frac)
      # spl$Max <- ifelse(spl$Freq >= Median_Frac, 1, 0)
    }else{
      spl$Max = ifelse((!is.na(max(spl$Freq))) & (max(spl$Freq) > 0) & (max(spl$Freq) == spl$Freq), 1, 0) 
    }
    
    ####
    spl
  }))
  classes$CellType = cell_type
  all_classes = rbind(all_classes, classes)
  
}

all_mapping$InitialID = paste(all_mapping$CellType, all_mapping$InitialState, sep = "_")
all_mapping$ID = paste(all_mapping$CellType, all_mapping$State, sep = "_")

write.table(all_mapping, file.path(path_Ecout, "mapping_all_states.txt"), sep = "\t", row.names = F)

casted = reshape2::dcast(all_classes, ID ~ CellType + State, value.var = "Max")
clusters = t(casted[,-1])
clusters[is.na(clusters)] = 0

clusters = clusters[match(all_mapping$InitialID, rownames(clusters)),]
rownames(clusters) = all_mapping$ID

colnames(clusters) = casted[,1]
write.table(clusters, file.path(path_Ecout, "binary_classification_all_states.txt"), sep = "\t")

jaccard = matrix(NA, nrow(clusters), nrow(clusters))
for(i in 1:(nrow(clusters)))
{
  for(j in (i):(nrow(clusters)))
  {
    int <- sum(clusters[i,] & clusters[j,])
    idx <- int / (sum(clusters[i,]) + sum(clusters[j,]) - int)
    
    p = 1 - phyper(int, sum(clusters[i,]), ncol(clusters) - sum(clusters[i,]), sum(clusters[j,]))
    if(is.na(p) | (p >= p_val_cutoff)){			
      idx = 0
    }		
    jaccard[i, j] = jaccard[j, i] = idx
  }
}
jaccard[is.na(jaccard)] = 0

rownames(jaccard) = colnames(jaccard) = rownames(clusters)
write.table(jaccard, file.path(path_Ecout, "jaccard_matrix.txt"), sep = "\t")

hclusCut <- function(x, k, ...) list(cluster = cutree(hclust(as.dist(1-x), method = "average", ...), k=k))
choose_clusters <- function(data, name, range = 2:10)
{	
  silh <- data.frame(K = range, Silhouette = sapply(range, function(k){
    sil <<- silhouette(hclusCut(data, k)$cluster, as.dist(1-data))
    tmp <<- summary(sil)
    tmp$avg.width
  }))
  
  g2 <- ggplot(silh, aes(x = K, y = Silhouette)) + 
    geom_point() +		
    geom_line() + 
    geom_vline(xintercept = silh[which.max(silh$Silhouette), 1], lty = 2, colour = "red") + 
    theme_bw() +
    theme(panel.grid = element_blank()) + 
    theme(aspect.ratio = 1) 		
  
  pdf(file.path(path_Ecout, paste0("nclusters_", name, ".pdf")), width = 13, height = 7)
  plot(g2)
  tmp = dev.off()
  png(file.path(path_Ecout, paste0("nclusters_", name, ".png")), width = 13, height = 7, res = 300, units = "in")
  plot(g2)
  tmp = dev.off()
  silh[which.max(silh$Silhouette), 1]	
}
hc = hclust(as.dist(1-jaccard), method = "average")
n_clust = choose_clusters(jaccard, "jaccard", range = 2:(nrow(jaccard) - 1))
clust =  hclusCut(jaccard, n_clust)$cluster

sil = silhouette(clust, as.dist(1-jaccard))
avg_silhouette = summary(sil)
write.table(avg_silhouette$avg.width, file.path(path_Ecout, "silhouette_initial.txt"), sep = "\t", row.names = F)

top_ann = as.data.frame(t(sapply(rownames(jaccard), function(x) {
  s= strsplit(x, "_")[[1]]
  c(paste0(s[-length(s)],collapse = "_"), s[length(s)])
})))

colnames(top_ann) = c("CellType","State")
top_ann$InitialEcotype = as.factor(sprintf("IE%02d", clust))

write.table(top_ann, file.path(path_Ecout, "initial_ecotypes.txt"), sep = "\t")

top_ann$ID = rownames(top_ann)
top_ann = top_ann[order(top_ann$InitialEcotype),]
write.table(top_ann, file.path(path_Ecout, "ecotypes.txt"), sep = "\t", row.names = F)

jaccard = jaccard[match(top_ann$ID, rownames(jaccard)), match(top_ann$ID, rownames(jaccard))]

top_ann$"Cell type" = top_ann$CellType
top_ann$"Cell type name" <- all_mapping$CellType_new[match(top_ann$ID, all_mapping$ID)]
#top_ann$"Cell type name" <- #

diag(jaccard) = 1
pdf(file.path(path_Ecout, "initial_jaccard_matrix.pdf"), width = 12, height = 7, family = "Helvetica")
h <- heatmap_simple(jaccard, name = "ht1", top_annotation = top_ann, top_columns = c("InitialEcotype", "Cell type", "Cell type name"), 
                    legend_name = "Jaccard index", width = unit(2, "in"), height = unit(2, "in"),
                    #color_palette = c("gray", viridis(4)), 
                    color_palette = c("gray", "#e9f4f6", "#8f9fab", "#4a5c73", "#1e2235"), 
                    raster_quality = 5,
                    color_range = c(0, 0.1, 0.2, 0.3))
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
     adjust_annotation_extension = T, merge_legends = T)	

ord = top_ann$InitialEcotype
dup = (which(!duplicated(ord)) - 1)
fract = dup / nrow(top_ann)
width =  c(fract[-1], 1) - fract
decorate_heatmap_body("ht1", {
  #grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 2))
})
tmp = dev.off()

initial_tb = table(top_ann$InitialEcotype)
tb = initial_tb[initial_tb >= min_states]
if(length(tb) < 2)
{
  warning(paste0("There are less than 2 ecotypes with more than ", min_states, " cell states. Including ecotypes of size ", min_states - 1, ".\n"))
  tb = initial_tb[initial_tb >= min_states - 1]
}

top_ann = top_ann[top_ann$InitialEcotype %in% names(tb),]

nm = unique(top_ann$InitialEcotype)
mapping = sprintf("E%d", 1:length(nm))
names(mapping) = nm

top_ann$Ecotype = mapping[as.character(top_ann$InitialEcotype)]
top_ann$Ecotype = ecotype_to_factor(top_ann$Ecotype)
top_ann = top_ann[order(top_ann$Ecotype),]
write.table(top_ann, file.path(path_Ecout, "ecotypes.txt"), sep = "\t", row.names = F)

jaccard = jaccard[match(top_ann$ID, rownames(jaccard)), match(top_ann$ID, rownames(jaccard))]

sil <- silhouette(as.numeric(as.character(gsub("E", "", as.character(top_ann$Ecotype)))), as.dist(1-jaccard))
avg_silhouette <<- summary(sil)
write.table(avg_silhouette$avg.width, file.path(path_Ecout, "silhouette.txt"), sep = "\t", row.names = F)

top_ann$"Cell type" = top_ann$CellType
top_ann$"Cell type name" <- all_mapping$CellType_new[match(top_ann$ID, all_mapping$ID)]

pdf(file.path(path_Ecout, "jaccard_matrix.pdf"), width = 12, height = 7, family = "Helvetica")
h <- heatmap_simple(jaccard, name = "ht1", top_annotation = top_ann, top_columns = c("Ecotype", "Cell type", "Cell type name"), 
                    legend_name = "Jaccard index", width = unit(2, "in"), height = unit(2, "in"),
                    #color_palette = c("gray", viridis(4)), 
                    color_palette = c("gray", "#e9f4f6", "#8f9fab", "#4a5c73", "#1e2235"), 
                    raster_quality = 5,
                    color_range = c(0, 0.1, 0.2, 0.3))
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
     adjust_annotation_extension = T, merge_legends = T)	

ord = top_ann$Ecotype
dup = (which(!duplicated(ord)) - 1)
fract = dup / nrow(top_ann)
width =  c(fract[-1], 1) - fract
decorate_heatmap_body("ht1", {
  grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 2))
})
tmp = dev.off()

png(file.path(path_Ecout, "jaccard_matrix.png"), width = 12, height = 7, res = 300, units = "in", family = "Helvetica")
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
     adjust_annotation_extension = T, merge_legends = T)	
decorate_heatmap_body("ht1", {
  grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 2))
})
tmp = dev.off()

# Assign samples new -----
path_root <- here("data/Correlation_and_crosstalk_analysis_18_07_25/EcoTyper_rerun_250827")
path_state <- file.path(path_root, "02_EcoPre")
path_Ecout <- file.path(path_root, "03_EcoDisc")

ecotypes = read.delim(file.path(path_Ecout, 'ecotypes.txt'), sep = '\t') 
mapping = read.delim(file.path(path_Ecout, 'mapping_all_states.txt'), sep = '\t') 
ecotypes$Ecotype = ecotype_to_factor(ecotypes$Ecotype) 
ecotypes = ecotypes[order(ecotypes$Ecotype),] 

all_H = NULL
all_classes_filt = NULL
cell_ids = c()
top_cols = c("Ecotype", "Datatype")

CT_all <- unique(Data_in$cell_type)

for(cell_type in CT_all){	
  mapping = read.delim(file.path(path_state, cell_type, 'mapping_to_initial_states.txt'))
  classes = read.delim(file.path(path_state, cell_type, 'initial_state_assignment.txt'))
  
  #clinical = read_clinical(classes$ID, dataset = dataset)
  #classes$Sample = clinical$Sample
  
  cell_ids = unique(c(cell_ids, classes$ID))
  
  classes = as.data.frame(table(classes$Sample, classes$State))
  
  colnames(classes) = c("ID", "State", "Freq")
  splits = split(classes, classes$ID)
  
  classes = do.call(rbind, lapply(splits, function(spl)
  {
    spl$Frac = spl$Freq / sum(spl$Freq)
    spl$Max= ifelse(which.max(spl$Freq) == spl$Freq, 1, 0)
    spl
  }))
  
  classes$CellType = cell_type	
  
  H = dcast(classes, State~ID, value.var = "Frac")
  rownames(H)  = H[,1]
  H = H[,-1,drop = F]
  H_raw = H
  H =H[match(mapping$InitialState, rownames(H)),,drop = F]
  rownames(H) = mapping$State
  
  rownames(H) = paste0(cell_type, "_", rownames(H))	
  keep_rowname = c(rownames(all_H), rownames(H))
  all_H = rbind.fill(all_H, H)
  rownames(all_H) = keep_rowname
  
  classes = as.data.frame(apply(H_raw, 2, function(x) {
    idx = which.max(x)
    if(length(idx)== 0)
    {
      "Unassigned"
    }else{
      rownames(H_raw)[idx]
    }
    
  }
  ))
  classes_raw = data.frame(ID = rownames(classes), InitialState = classes[,1])
  classes_raw = classes_raw[classes_raw$InitialState %in% mapping$InitialState,]
  classes_raw$State = mapping[match(classes_raw$InitialState, mapping$InitialState), "State"]
  classes = classes_raw
  classes = classes[,c("ID", "State")]
  clusters = ecotypes[ecotypes$CellType == cell_type,] 
  classes = classes[classes$State %in% clusters$State,]
  
  colnames(classes) = c('ID', cell_type)
  
  if(is.null(all_classes_filt))
  {
    all_classes_filt = classes
  }else{
    all_classes_filt = merge(all_classes_filt, classes, by = 'ID', all = T)
  }
} 

all_H = all_H[match(ecotypes$ID, rownames(all_H)),]
write.table(all_H, file.path(path_Ecout, "combined_state_abundances.txt"), sep = "\t")

H = do.call(rbind, lapply(levels(ecotypes$Ecotype), function(clst){
  clst <<- clst
  #print(clst)
  inc <<- ecotypes[ecotypes$Ecotype == clst,]$ID
  apply(all_H[rownames(all_H) %in% inc,,drop = F], 2, mean, na.rm = T)
}))

rownames(H) = levels(ecotypes$Ecotype)
#write.table(H, file.path(path_Ecout, "ecotype_abundance.txt"), sep = "\t")
H = apply(H, 2, function(x) x / sum(x, na.rm = T))
write.table(H, file.path(path_Ecout, "ecotype_abundance.txt"), sep = "\t")

p_vals = do.call(rbind, lapply(levels(ecotypes$Ecotype), function(clst){
  clst <<- clst
  #print(clst)
  inc <<- ecotypes[ecotypes$Ecotype == clst,]$ID
  
  apply(all_H, 2, function(x){
    x <<- x
    err <<- F
    p <<- 1
    tryCatch({
      p <<- t.test(x[rownames(all_H) %in% inc], x[!(rownames(all_H) %in% inc)])$p.value
    }, error = function(x) err <<- T)
    if(err)
    {
      p = NA 
    }
    p
  })
  
})) 
rownames(p_vals) = levels(ecotypes$Ecotype)
write.table(p_vals, file.path(path_Ecout, "assignment_p_vals.txt"), sep = "\t")

assignment = as.data.frame(apply(H, 2, function(x) {
  idx = which.max(x)
  if(length(idx)==0)
  {
    "Unassigned"
  }else{
    rownames(H)[idx]
  }
}))

clinical = data.frame(ID = rownames(assignment), MaxEcotype = assignment[,1])
clinical$AssignmentP = sapply(1:ncol(H), function(i) {
  idx = which.max(H[,i])
  if(length(idx)==0)
  {
    NA
  }else{
    p_vals[idx, i] 
  }		
})
clinical$AssignmentQ = p.adjust(clinical$AssignmentP, method = "BH")

clinical$AssignedToEcotypeStates = clinical$ID %in% all_classes_filt$ID

#clinical$Ecotype = ifelse((clinical$AssignmentQ < 0.25) & clinical$AssignedToEcotypeStates, as.character(clinical$MaxEcotype), "Unassigned")
clinical$Ecotype = ifelse( clinical$AssignedToEcotypeStates, as.character(clinical$MaxEcotype), "Unassigned")
clinical$Ecotype = factor(as.character(clinical$Ecotype), levels = c(levels(ecotypes$Ecotype), "Unassigned"))

annotation <- as.data.frame(Data_in)
annotation <- unique(annotation[c("Sample", "Datatype")])
additional_columns = top_cols[top_cols %in% colnames(annotation)]

agg_ann = do.call(cbind, sapply(additional_columns, function(x){
  tb = as.data.frame(table(annotation$Sample, annotation[,x]))
  tb = tb[tb[,3] > 0,]
  if(sum(duplicated(tb[,1])) > 0){
    warning(paste0("Not plotting column '", x, "', as it has multiple distinct values within the same sample!"))
    return(NULL)
  }
  tb = tb[match(clinical$ID, tb[,1]),]
  as.character(tb[,2])
}, simplify = F))

if(!is.null(agg_ann))
{
  clinical = cbind(clinical, agg_ann)
  top_cols = unique(c(colnames(agg_ann), "Ecotype"))	
}else{
  top_cols = "Ecotype"
}

#tmp = read_clinical(clinical$ID, dataset = dataset)
#to_rem = colnames(tmp)[colnames(tmp) %in% colnames(clinical)]
#to_rem = to_rem[to_rem != "ID"]
#tmp = tmp[,!colnames(tmp) %in% to_rem]
#clinical = merge(clinical, tmp, by = "ID", all.x = T)

clinical = clinical[order(clinical$Ecotype),]
H = H[,match(clinical$ID, colnames(H))]
all_H = all_H[,match(clinical$ID, colnames(all_H))]
all_H = all_H[match(ecotypes$ID, rownames(all_H)),]

write.table(clinical, file.path(path_Ecout, "initial_ecotype_assignment.txt"), sep = "\t")

rownames(clinical) = clinical$ID
rownames(ecotypes) = ecotypes$ID

h <- heatmap_simple(all_H, top_annotation = clinical, top_columns = top_cols, 
                    left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State", "Cell.type.name"),
                    column_split = ifelse(clinical$Ecotype == "Unassigned", "Unassigned", "Assigned"),
                    width = unit(7, "in"), height = unit(4, "in"),
                    legend_name = "State abundance",
                    color_range = seq(0, quantile(as.matrix(all_H), .9, na.rm = T), 
                                      length.out = 8), 
                    color_palette = c("gray", "#e9f4f6","#c8d8df","#a9b9c6","#8f9fab",
                                      "#748293","#5a697c","#4a5c73","#1e2235"), 
                    #color_palette = c("gray", viridis(8)), 
                    raster_quality = 5)

pdf(file.path(path_Ecout, "heatmap_all_samples.pdf"), width = 12, height = 9)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
tmp = dev.off()

clinical_filt = clinical[clinical$Ecotype != "Unassigned",]
clinical_filt$Ecotype = factor(as.character(clinical_filt$Ecotype), levels = levels(ecotypes$Ecotype))
clinical_filt = clinical_filt[order(clinical_filt$Ecotype),]
write.table(clinical_filt, file.path(path_Ecout, "ecotype_assignment.txt"), sep = "\t")
small_H = as.matrix(all_H[,match(clinical_filt$ID, colnames(all_H))])
if(is.null(small_H) || nrow(small_H) == 0 || ncol(small_H) == 0)
{
  stop(paste("No samples were assigned to ecotypes!"))
}

h = heatmap_simple(small_H, top_annotation = clinical_filt, top_columns = top_cols, 
                   left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State", "Cell.type.name"),
                   width = unit(5, "in"), height = unit(3, "in"),
                   legend_name = "State abundance",
                   color_palette = c("gray", "#e9f4f6","#c8d8df","#a9b9c6","#8f9fab",
                                     "#748293","#5a697c","#4a5c73","#1e2235"), 
                   color_range = seq(0, quantile(as.matrix(all_H), .9, na.rm = T), length.out = 8), 
                   #color_palette = c("gray", viridis(8)), 
                   raster_quality = 20)

pdf(file.path(path_Ecout, "heatmap_assigned_samples_viridis.pdf"), width = 12, height = 8)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
#draw(h1, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
suppressWarnings({
  rect = rectangle_annotation_coordinates(ecotypes$Ecotype, clinical_filt$Ecotype)
})
decorate_heatmap_body("hmap", {
  grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3))
})
tmp = dev.off()

png(file.path(path_Ecout, "heatmap_assigned_samples_viridis.png"), width = 12, height = 8, units = "in", res = 300)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
  grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3)) 
})
tmp = dev.off()

h = heatmap_simple(small_H, top_annotation = clinical_filt, top_columns = top_cols, 
                   left_annotation = ecotypes, left_columns = c("Ecotype", "CellType", "State", "Cell.type.name"),
                   width = unit(5, "in"), height = unit(3, "in"), 
                   legend_name = "State abundance",
                   color_range = c(seq(0, quantile(small_H, .8, na.rm = T), length.out = 8)), color_palette = c("gray", brewer.pal(8, "YlGnBu")), raster_quality = 20)

pdf(file.path(path_Ecout, "heatmap_assigned_samples_YlGnBu.pdf"), width = 12, height = 8)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
  grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3))
})
tmp = dev.off()

png(file.path(path_Ecout, "heatmap_assigned_samples_YlGnBu.png"), width = 12, height = 8, units = "in", res = 300)
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)	
decorate_heatmap_body("hmap", {
  grid.rect(x = unit(rect$x, "native"), y = unit(rect$y, "native"), width = unit(rect$w, "native"), height = unit(rect$h, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 3))
})
tmp = dev.off()

# EcoTyper align -----
path_root <- here("data/Correlation_and_crosstalk_analysis_18_07_25/EcoTyper_rerun_250827")
path_state <- file.path(path_root, "02_EcoPre")
path_Ecout <- file.path(path_root, "03_EcoDisc")

ecotypes = read.delim(file.path(path_Ecout, 'ecotypes.txt'), sep = '\t')
ecotypes$Ecotype = ecotype_to_factor(ecotypes$Ecotype)
ecotypes = ecotypes[order(ecotypes$Ecotype),]

CT_all <- unique(Data_in$cell_type)
# CT_all <- c("Bcell")

for(cell_type in CT_all){	
  mapping = read.delim(file.path(path_state, cell_type, 'mapping_to_initial_states.txt'))
  classes = read.delim(file.path(path_state, cell_type, 'initial_state_assignment.txt'))
  
  State_list <- unique(ecotypes[ecotypes$CellType == cell_type, "State"])
  
  for(StateNM in State_list){
    
    ecotypes[(ecotypes$CellType == cell_type) & (ecotypes$State == StateNM), "CellType_new"] = mapping[mapping$State == StateNM, "CellType_new"]
    
  }
  
}

ecotypes <- ecotypes[c("CellType", "State", "InitialEcotype", "ID", "Cell.type", "Cell.type.name", "Ecotype")]
file_out <- file.path(path_root, "03_EcoDisc", "ecotypes_alignCT.txt")
write.table(ecotypes, file_out, sep = "\t", quote = FALSE, row.names = FALSE)

# Ecotype abundance boxplot ----- 
library(ggsignif)
path_root <- here("data/Correlation_and_crosstalk_analysis_18_07_25/EcoTyper_rerun_250827")

Adund_file <- file.path(path_root, "03_EcoDisc", "ecotype_abundance.txt")
Adund_data <- read.table(Adund_file, sep = "\t")
Adund_data <- as.data.frame(t(Adund_data))
Ecotype_list <- names(Adund_data)
Adund_data["SampleID"] <- row.names(Adund_data)

#SampleMeta_file <- file.path(path_root, "09_Ecotyper/01_preData", "SampleINFO.tsv")
#SampleMeta_data <- read.table(SampleMeta_file, sep = "\t", header = TRUE)
SampleMeta_data <- as.data.frame(Data_in)
SampleMeta_data <- unique(SampleMeta_data[c("Sample", "Datatype")])
names(SampleMeta_data) <- c("SampleID", "Group")

MergeData <- merge(Adund_data, SampleMeta_data, by = "SampleID", all.x = TRUE)
MergeData$Group <- factor(MergeData$Group, levels = c("A", "B", "C", "D"))

pplist <- list()
for(EcotypeNM in Ecotype_list){
  
  temp_data <- MergeData[c(EcotypeNM, "Group")]
  names(temp_data) <- c("Value", "Group")
  
  P1 <- ggplot(temp_data, mapping = aes(Group, Value, fill = Group)) +
    theme_classic(base_size = 20) +
    geom_boxplot() +
    ylab(EcotypeNM) +
    #scale_fill_manual(values = c("#66c2a5", "#8da0cb")) +
    theme(axis.text = element_text(colour = "black")) +
    geom_signif(comparisons = list(c("A", "B"), c("A", "C"), c("B", "D")), 
                map_signif_level=T,  
                label="p.signif", 
                test = "t.test", step_increase = 0.1) 
  
  print(P1)
  pplist[[EcotypeNM]] <- P1
}

cowplot::plot_grid(pplist[[1]], 
                   pplist[[2]], 
                   pplist[[3]]#, 
                   #pplist[[4]]#, 
                   #pplist[[5]]#, 
                   #pplist[[6]], 
                   #pplist[[7]], 
                   #pplist[[8]], 
                   #pplist[[9]], 
                   #pplist[[10]], 
                   #pplist[[11]], 
                   #pplist[[12]], 
                   #pplist[[13]], 
                   #pplist[[14]]#,# 
                   #pplist[[15]], 
)

# install.packages("ggtern") 
library(ggtern)
library(tidyverse)

## Ternary plot  ----- 
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

TNR_maping <- data.frame(seuratObj$sample, seuratObj$TNR_val) %>% unique()
data_for_tri <- MergeData
data_for_tri$TNR_val = TNR_maping$seuratObj.TNR_val[match(data_for_tri$SampleID, TNR_maping$seuratObj.sample)]

data_for_tri$TNR_val <- as.numeric(as.character(data_for_tri$TNR_val))
data_for_tri$Treatment_status <- ifelse(data_for_tri$Group %in% c("A", "C"), "Pre", "Post")
p1 <- ggtern(data = data_for_tri, aes(x = E1, y = E2, z = E3)) +
  geom_point(aes(size = ifelse(Treatment_status == "Post", TNR_val, 30),
                 color = Group, 
                 stroke = NA),
            alpha = 0.8) + 
  scale_size_continuous(range = c(3, 7)) + 
  scale_color_manual(values = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))) +
  theme_rgbw() #+ geom_crosshair_tern(size=1)
p1

pdf(file.path(path_Ecout, "ecotype_tern_plot_all_scRNA_samples.pdf"), width = 8, height = 6)
p1
tmp = dev.off()
