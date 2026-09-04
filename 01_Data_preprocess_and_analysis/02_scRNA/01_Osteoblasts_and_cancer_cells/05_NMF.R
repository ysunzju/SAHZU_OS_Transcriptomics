library(here) # project-root-relative paths; run scripts from repository root
##################### NMF related ##################### 
##################### Metaprogram generation (related to Figure 3A) ##################### 
## Adapted from https://github.com/tiroshlab/3ca/tree/main/ITH_hallmarks
### Source functions can be found in the 03_Doc folder or from the original Github repo 
### Courtesy of Dr Z Hu, SYSU
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/NMF_results_28_06_25"))
list_names <- dir()
list_results <- list()
for(i in list_names){
  name1 <- strsplit(i, '\\.')
  name2 <- paste0(name1[[1]][1], '_rank4_9_nrun10.RDS')	
  tmp <- readRDS(i)
  f <- (is.na((tmp$fit)))
  if(length(f[f==FALSE])==6){
    tmp4 <- tmp$fit[['4']]@fit@W
    tmp5 <- tmp$fit[['5']]@fit@W
    tmp6 <- tmp$fit[['6']]@fit@W
    tmp7 <- tmp$fit[['7']]@fit@W
    tmp8 <- tmp$fit[['8']]@fit@W
    tmp9 <- tmp$fit[['9']]@fit@W
    tmp_all <- cbind(tmp4, tmp5, tmp6, tmp7, tmp8, tmp9)
    colnames(tmp_all) <- paste(name2, c(seq(4.1,4.4,0.1), seq(5.1,5.5,0.1), seq(6.1,6.6,0.1), seq(7.1,7.7,0.1), seq(8.1,8.8,0.1), seq(9.1,9.9,0.1)), sep = '.')
    list_results[[name2]] <- tmp_all
  }
}

################
library(reshape2)
library(NMF)
library(ggplot2)
library(scales)
source(here("01_Data_preprocess_and_analysis/02_scRNA/01_Osteoblasts_and_cancer_cells/custom_magma.R"))
source(here("01_Data_preprocess_and_analysis/02_scRNA/01_Osteoblasts_and_cancer_cells/robust_nmf_programs.R"))

## Parameters 
intra_min_parameter <- 35 
intra_max_parameter <- 5 
inter_min_parameter <- 10 

Genes_nmf_w_basis <- list_results
# get top 50 genes for each NMF program 
nmf_programs <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 

# for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
nmf_filter_ccle <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter=T, inter_min = inter_min_parameter)  
nmf_programs <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs <- do.call(cbind, nmf_programs)

# calculate similarity between programs
nmf_intersect <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# ----------------------------------------------------------------------------------------------------
# Cluster selected NMF programs to generate MPs
# ----------------------------------------------------------------------------------------------------

### Parameters for clustering
Min_intersect_initial <- 5    # the minimal intersection cutoff for defining the first NMF program in a cluster
Min_intersect_cluster <- 10    # the minimal intersection cutoff for adding a new NMF to the forming cluster 
Min_group_size <- 5     # the minimal group size to consider for defining the first NMF program in a cluster 

#Min_intersect_initial <- 15    # the minimal intersection cutoff for defining the first NMF program in a cluster
#Min_intersect_cluster <- 15    # the minimal intersection cutoff for adding a new NMF to the forming cluster 
#Min_group_size <- 10     # the minimal group size to consider for defining the first NMF program in a cluster 
# Was 10, 10, 5

Sorted_intersection <- sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list <- list()   ### Every entry contains the NMFs of a chosen cluster
MP_list <- list()
k <- 1
Curr_cluster <- c()

nmf_intersect_original <- nmf_intersect

while (Sorted_intersection[1]>Min_group_size) {  
  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP                    <- nmf_programs[,names(Sorted_intersection[1])] # Genes in the forming MP are first chosen to be those in the first NMF. Genes_MP always has only 50 genes and evolves during the formation of the cluster
  nmf_programs                <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]  # remove selected NMF
  Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history                 <- Genes_MP  # has genes in all NMFs in the current cluster, for redefining Genes_MP after adding a new NMF 
  
  ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.    
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {  
    
    Curr_cluster <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    
    Genes_MP_temp <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene
    
    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste( strsplit(i , "[.]")[[1]][1 : which(strsplit(i , "[.]")[[1]] == "RDS")]   , collapse = "."  )
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   
      
      Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    
    NMF_history <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP <- Genes_MP_temp[1:50]
    
    nmf_programs <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]  # remove selected NMF
    
    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]] <- Genes_MP
  k <- k+1
  
  nmf_intersect <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # Remove current chosen cluster
  
  Sorted_intersection <- sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   # Sort intersection of remaining NMFs not included in any of the previous clusters
  
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

####  Sort Jaccard similarity plot according to new clusters: 723
#### Manual rearrangement of MP and cluster sequence 
MP_list <- MP_list[c(3, 9, 1, 2, 8, 5, 6, 4, 7)]
Cluster_list <- Cluster_list[c(3, 9, 1, 2, 8, 5, 6, 4, 7)]

inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  
  inds_sorted <- c(inds_sorted ,match(Cluster_list[[j]] ,colnames(nmf_intersect_original)))
  
}

inds_new <- c(inds_sorted, which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered

nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new]) 

p <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile(color = NA) + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  coord_fixed() + 
  theme(axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  #guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  theme(legend.position = "none")
#ggsave("Jaccard.pdf", plot = p, width = 12, height = 10)
p

pdf("Cancer_cell_NMF.pdf", height = 5, width = 5)
print(p)
dev.off()

anno <- data.frame(levels = levels(p$data$Var1), group = substr(levels(p$data$Var1), 1, 4))
anno$group <- factor(anno$group, levels = c("H_BF", "H_AF", "L_BF", "L_AF"))
anno$levels <- factor(anno$levels, levels = levels(nmf_intersect_meltI_NEW$Var1))
anno$strat <- "dummy"
p1 <- ggplot(anno, aes(x = strat, y = levels, fill = group)) + 
  geom_tile(aes(x = strat, y = levels, fill = group), color = NA) + 
  scale_fill_manual(values = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))) + 
  coord_equal() + 
  theme_bw() + 
  theme(panel.border = element_rect(size = 1, color = "black", fill = NA), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none")
#p1 <- ggplotify::as.ggplot(p1)
p1

pdf("Cancer_cell_NMF_anno.pdf", height = 5, width = 5)
print(p1)
dev.off()

### Save updated MP list and Cluster list 
saveRDS(MP_list, here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
#MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
#saveRDS(Cluster_list, here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_Cluster_list.rds"))

##################### MP assignment (related to Figure 3F) ##################### 
library(ggplot2)
library(scalop) 
library(gridExtra)
library(ggpubr)
library(Seurat)

### Define the following:
SeuratObj <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Malig_OB_for_NMF_sample_w_too_few_cells_removed.rds"))
DefaultAssay(SeuratObj) <- 'RNA'
SeuratObj.list <- SplitObject(SeuratObj, split.by = "sample")
list1 <- list()
for(i in names(SeuratObj.list)){
  counts <- as.matrix(SeuratObj.list[[i]][["RNA"]]$counts)
  list1[[i]] <- counts
}

My_study = list1
MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))

L = My_study
cell_type = "Cancer"
MinGenes = 25 # MinGenes  = the min number of MP genes that should exist in the study (default is 25)
MinScore = 0.05 # MinScore  = the minimal score for assigning a cell to a MP (default is 1). A cell is assigned to the MP with the maximal score, given that it exceeded MinScore. If the maximal score is below MinScore, the cell is unassigned. 
#MinScore = 1 
MinCells = 0   #0.05 # MinCells  = the minimal % of cells that should be assigned to a MP in order to account for the MP (default is 0.05)

#MP_list = MP_list 

#MP_list <- MP_list[[cell_type]]

L <- lapply(L, function(x) log2((x/10) + 1) ) # Log normalize before scoring 

MP_scores_per_sample <- lapply(L, function(x)  scalop::sigScores(as.matrix(x), sigs = MP_list, conserved.genes = MinGenes/50 ))
#saveRDS(MP_scores_per_sample, "MP_scores_per_sample.RDS")
#MP_scores_per_sample <- readRDS("MP_scores_per_sample.RDS")

remove_cells <- function(x,MinScore){         # remove cells whose max score was below MinScore
  max_score <- apply(x, 1, function(y) max(y))
  cells_rm  <- which(max_score < MinScore)
  if (length(cells_rm) > 0 ){
    x <- x[-cells_rm , ]
  }
  return(x)
}
MP_scores_per_sample <- lapply(MP_scores_per_sample, function(x) remove_cells(x,MinScore))

Assign_MP_per_cell   <- lapply(MP_scores_per_sample, function(x) apply(x, 1, function(y) colnames(x)[which(y==max(y))] ) ) 

filter_cells <- function(x,MinCells){         # remove MP that were assassin to less than MinCells in each sample
  MP_frequency <- as.numeric(ave(x, x, FUN = length))
  MP_rm        <- which(MP_frequency/length(x) < MinCells)  # MPs to be removed
  if (length(MP_rm)>0){
    x <- x[-MP_rm]
  }
  return(x)
}
Assign_MP_per_cell_filtered <- lapply(Assign_MP_per_cell, function(x) filter_cells(x,MinCells)) 

Assign_MP_per_cell_filtered1 <- data.frame()
for(i in names(Assign_MP_per_cell_filtered)){
  tmp <- data.frame(Assign_MP_per_cell_filtered[[i]])
  rownames(tmp) <- names(Assign_MP_per_cell_filtered[[i]])
  Assign_MP_per_cell_filtered1 <- rbind(Assign_MP_per_cell_filtered1, tmp)
}
colnames(Assign_MP_per_cell_filtered1) <- 'MP_group'

diffsample <- setdiff(rownames(SeuratObj@meta.data), rownames(Assign_MP_per_cell_filtered1))
Assign_MP_per_cell_filtered2 <- data.frame(rep('Unresolved', length(diffsample)))
rownames(Assign_MP_per_cell_filtered2) <- diffsample
colnames(Assign_MP_per_cell_filtered2) <- 'MP_group'
Assign_MP_per_cell_filtered3 <- rbind(Assign_MP_per_cell_filtered1, Assign_MP_per_cell_filtered2)
table(Assign_MP_per_cell_filtered3$MP_group)

Cancer.OS.combined <- SeuratObj
Cancer.OS.combined <- AddMetaData(Cancer.OS.combined, Assign_MP_per_cell_filtered3)
Cancer.OS.combined1 <- Cancer.OS.combined

Idents(Cancer.OS.combined1) <- 'MP_group'
levels(Cancer.OS.combined1) <- c(paste0('MP_', c(1:9)), 'Unresolved')

Cancer.OS.combined1$MP_group <- factor(Cancer.OS.combined1$MP_group, levels = c(paste0('MP_', c(1:9)), 'Unresolved'))

library(scop)
p1 <- SCP::CellDimPlot(
  srt = Cancer.OS.combined1, group.by = "MP_group", 
  #stat.by = "group_anno",
  split.by = "MP_group",
  #cols.highlight = new_colors,
  reduction = "umap", theme_use = "theme_blank", label = F, label_repel = F, 
  palcolor = c("#a97c50", "#2e3192", "#010101", "#972D15","#FBB4AE","#CCEBC5","#DECBE4","#FED9A6",
               "#7f3f98", "#B3CDE3"# "#339999", 
  ), 
  raster = T, 
  #ncol = 5,
  pt.size = 3, 
  combine = T, 
  #split.by = "treatment_status", 
  label_insitu = F) #+ theme(legend.position = "none")

ggsave("MP_group_split.pdf", plot = pps, width = 30, height = 20, unit = 'in')
# Please note that the diffusely distributed unassigned cells (i.e., expressions not meeting the minumum score for assignment to any of the MPs) are not displayed in the figure panel in the manuscript

rm(My_study)
rm(L)
rm(list1)
gc()

##################### Cellular state plot (related to Figure 3E) ##################### 
List_for_plot <- MP_list[c(1:4)]
names(List_for_plot) <- c("G2_M", "S", "Os_min", "ECM")
Idents(SeuratObj) <- "sample"

L_idents <- c("L_BF_P2", "L_BF_P5", "L_BF_P6", "L_AF_P2", "L_AF_P5", "L_AF_P6")
H_idents <- c("H_BF_P1", "H_BF_P3", "H_BF_P5", "H_AF_P1", "H_AF_P3", "H_AF_P5")
scRNA_for_plot <- subset(SeuratObj, idents = c(L_idents, H_idents))
out <- SCpubr::do_CellularStatesPlot(sample = scRNA_for_plot,
                                     input_gene_list = List_for_plot,
                                     x1 = "G2_M",
                                     y1 = "S",
                                     x2 = "Os_min",
                                     y2 = "ECM",
                                     pt.size = 0.1, 
                                     plot_cell_borders = TRUE,
                                     enforce_symmetry = T,
                                     group.by = "group_anno", 
                                     #plot_features = TRUE,
                                     #plot_enrichment_scores = T, 
                                     plot_features = F#,
                                     #features = c("PC_1", "nFeature_RNA")
)
out

#install.packages("ggdensity")
#install.packages("ggblanket")
library(ggdensity)
library(ggblanket)
library(ggsci)

data_for_plot <- out$main$data
data_for_plot <- out$data
custom_magma <- colorRampPalette(c("#e9f4f6", "#8f9fab", "#1e2235"))(333)

p1 <- ggplot(data_for_plot,aes(x = set_x, y = set_y, fill = group.by)) +
  geom_hdr() + 
  scale_fill_manual(values = rep("#1e2235", 4)) + 
  #scale_fill_manual(values = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))) + 
  #geom_hdr_lines() + 
  #scale_fill_discrete("black") + 
  #geom_point(shape = 16, alpha = 0.1) +
  #scale_fill_aaas() +
  #scale_fill_continuous(custom_magma) + 
  #scale_fill_continuous(custom_magma) + 
  scale_color_gradientn(colors = custom_magma) + 
  coord_fixed() + 
  facet_wrap(vars(group.by)) +
  coord_equal() + 
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line())

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "_Cancer_cell_state_v1.pdf"), 
    width = 8.5, # The width of the plot in inches
    height = 7) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()

##################### CytoTRACE score boxplot and slingshot analysis (related to Figure 3G,H) ##################### 
## CytoTRACE confirmation ----- 
library(CytoTRACE2)
Cancer_OB_Cytotrace_metadata <- readRDS(here("data/NMF/Cancer_OB_Cytotrace_metadata_15_04_25.rds"))
Cancer.OS.combined1 <- AddMetaData(Cancer.OS.combined1, Cancer_OB_Cytotrace_metadata)
names(Cancer.OS.combined1@meta.data)

table(Cancer.OS.combined1$CytoTRACE2_Potency, Cancer.OS.combined1$Cell_type_fine_harmony)

# generate prediction and phenotype association plots with plotData function
expression_data <- GetAssayData(Cancer.OS.combined1)
annotation <- as.data.frame(Cancer.OS.combined1$MP_group)

plots <- plotData(cytotrace2_result = Cancer.OS.combined1, 
                  annotation = annotation,
                  expression_data = expression_data, 
                  is_seurat = T)

# plot
p1  <- plots$CytoTRACE2_Boxplot_byPheno +
  scale_fill_gradient2(
    low = "#50859f", high = "#d66692", mid = "white",
    midpoint = 0.5, limit = c(0, 1), #name = "Correlation",
    labels = scales::number_format(accuracy = 0.1)
  )
p1@layers$geom_jitter <- NULL
p1

pdf("Cancer_cell_CytoTRACE_Cell_potency_cat_260219.pdf", height = 3, width = 4)
print(p1)
dev.off()

#rm(Cancer.OS.combined1)
rm(expression_data)
gc()

## Slingshot ----- 
library(scop)
library(slingshot)
library(Seurat)
library(SingleCellExperiment)
library(RColorBrewer)

sc <- as.SingleCellExperiment(Cancer.OS.combined1)
sc <- slingshot(sc, clusterLabels = "MP_group", 
                start.clus = c("MP_3")#, reducedDim = "umap"
)

lin1 <- getLineages(sc, 
                    clusterLabels = "MP_group", 
                    reducedDim = "UMAP")
plot(reducedDims(sc)$UMAP,col = c("#a97c50", "#2e3192", "#010101", "#972D15","#FBB4AE","#CCEBC5","#DECBE4","#FED9A6",
                                  "#7f3f98", "#B3CDE3")[Cancer.OS.combined1$MP_group],pch=16,asp=1)
lines(SlingshotDataSet(lin1), lwd=2,col = 'black',type = 'lineages')

## Inferred lineage-specific trajectories plotting ---- 
library(ggplot2)
library(viridis)
library(ggrastr)

### L1 -----
sling_pseudotime <- sc$slingPseudotime_1
sling_meta <- data.frame(
  row.names = colnames(sc),
  slingshot_pseudotime = sling_pseudotime
)
identical(rownames(sling_meta), rownames(Cancer.OS.combined1@meta.data))
Cancer.OS.combined1 <- AddMetaData(Cancer.OS.combined1, metadata = sling_meta)
head(Cancer.OS.combined1@meta.data)

# curv info
sds <- SlingshotDataSet(sc)
lineage1_coords <- as.data.frame(sds@curves$Lineage1$s)
colnames(lineage1_coords) <- c("umap_1", "umap_2")

cell_coords <- as.data.frame(Embeddings(Cancer.OS.combined1, reduction = "umap"))
cell_coords$slingshot_pseudotime <- Cancer.OS.combined1$slingshot_pseudotime

p_slingshot_umap <- ggplot() +
  geom_point_rast(
    data = filter(cell_coords, is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2),
    size = 1,
    color = "grey80" 
  ) + 
  geom_point_rast(
    data = filter(cell_coords, !is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2, color = slingshot_pseudotime),
    size = 1
  ) +
  scale_color_gradientn(colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                   "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                   "#d73027", "#a50026"), 
                        na.value = "grey80") + 
  geom_path(
    data = lineage1_coords,
    aes(x = umap_1, y = umap_2),
    color = "black", 
    linewidth = 1.2    
  ) +
  labs(
    title = "Slingshot Trajectory on UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1, 
    legend.position = "top"
  )+
  guides(
    color = guide_colorbar(
      barwidth = unit(3, "cm"),  
      barheight = unit(0.5, "cm"), 
      title = "Pseudotime",
      title.theme = element_text(size = 11, face = "bold"), 
      title.vjust = 0.8, 
      label.theme = element_text(size = 8) 
    )
  ) 

p_slingshot_umap
ggsave("slingshot_L1.pdf", p_slingshot_umap, width = 5, height = 5.5)

# L2 -----
sling_pseudotime <- sc$slingPseudotime_2
sling_meta <- data.frame(
  row.names = colnames(sc),
  slingshot_pseudotime = sling_pseudotime
)
identical(rownames(sling_meta), rownames(Cancer.OS.combined1@meta.data))
Cancer.OS.combined1 <- AddMetaData(Cancer.OS.combined1, metadata = sling_meta)
head(Cancer.OS.combined1@meta.data)

# curv info
sds <- SlingshotDataSet(sc)
lineage1_coords <- as.data.frame(sds@curves$Lineage2$s)
colnames(lineage1_coords) <- c("umap_1", "umap_2")

cell_coords <- as.data.frame(Embeddings(Cancer.OS.combined1, reduction = "umap"))
cell_coords$slingshot_pseudotime <- Cancer.OS.combined1$slingshot_pseudotime

p_slingshot_umap <- ggplot() +
  geom_point_rast(
    data = filter(cell_coords, is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2),
    size = 1,
    color = "grey80" 
  ) + 
  geom_point_rast(
    data = filter(cell_coords, !is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2, color = slingshot_pseudotime),
    size = 1
  ) +
  scale_color_gradientn(colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                   "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                   "#d73027", "#a50026"), 
                        na.value = "grey80") + 
  geom_path(
    data = lineage1_coords,
    aes(x = umap_1, y = umap_2),
    color = "black",  
    linewidth = 1.2      
  ) +
  labs(
    title = "Slingshot Trajectory on UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"   
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1, 
    legend.position = "top"
  )+
  guides(
    color = guide_colorbar(
      barwidth = unit(3, "cm"),  
      barheight = unit(0.5, "cm"), 
      title = "Pseudotime",
      title.theme = element_text(size = 11, face = "bold"),   
      title.vjust = 0.8, 
      label.theme = element_text(size = 8)   
    )
  ) #+
#SeuratExtend::theme_umap_arrows()

p_slingshot_umap
ggsave("slingshot_L2.pdf", p_slingshot_umap, width = 5, height = 5.5)

# L3 -----
sling_pseudotime <- sc$slingPseudotime_3
sling_meta <- data.frame(
  row.names = colnames(sc),
  slingshot_pseudotime = sling_pseudotime
)
identical(rownames(sling_meta), rownames(Cancer.OS.combined1@meta.data))
Cancer.OS.combined1 <- AddMetaData(Cancer.OS.combined1, metadata = sling_meta)
head(Cancer.OS.combined1@meta.data)

# curv info
sds <- SlingshotDataSet(sc)
lineage1_coords <- as.data.frame(sds@curves$Lineage3$s)
colnames(lineage1_coords) <- c("umap_1", "umap_2")

cell_coords <- as.data.frame(Embeddings(Cancer.OS.combined1, reduction = "umap"))
cell_coords$slingshot_pseudotime <- Cancer.OS.combined1$slingshot_pseudotime

p_slingshot_umap <- ggplot() +
  geom_point_rast(
    data = filter(cell_coords, is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2),
    size = 1,
    color = "grey80" 
  ) +
  geom_point_rast(
    data = filter(cell_coords, !is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2, color = slingshot_pseudotime),
    size = 1
  ) +
  scale_color_gradientn(colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                   "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                   "#d73027", "#a50026"), 
                        na.value = "grey80") + 
  geom_path(
    data = lineage1_coords, 
    aes(x = umap_1, y = umap_2),
    color = "black",  
    linewidth = 1.2      
  ) +
  labs(
    title = "Slingshot Trajectory on UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"   
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1,   
    legend.position = "top"
  )+
  guides(
    color = guide_colorbar(
      barwidth = unit(3, "cm"),  
      barheight = unit(0.5, "cm"), 
      title = "Pseudotime",
      title.theme = element_text(size = 11, face = "bold"),   
      title.vjust = 0.8, 
      label.theme = element_text(size = 8)   
    )
  ) #+
#SeuratExtend::theme_umap_arrows()

p_slingshot_umap
ggsave("slingshot_L3.pdf", p_slingshot_umap, width = 5, height = 5.5)

# L4 -----
sling_pseudotime <- sc$slingPseudotime_4
sling_meta <- data.frame(
  row.names = colnames(sc),
  slingshot_pseudotime = sling_pseudotime
)
identical(rownames(sling_meta), rownames(Cancer.OS.combined1@meta.data))
Cancer.OS.combined1 <- AddMetaData(Cancer.OS.combined1, metadata = sling_meta)
head(Cancer.OS.combined1@meta.data)

# curv info
sds <- SlingshotDataSet(sc)
lineage1_coords <- as.data.frame(sds@curves$Lineage4$s)
colnames(lineage1_coords) <- c("umap_1", "umap_2")

cell_coords <- as.data.frame(Embeddings(Cancer.OS.combined1, reduction = "umap"))
cell_coords$slingshot_pseudotime <- Cancer.OS.combined1$slingshot_pseudotime

p_slingshot_umap <- ggplot() +
  geom_point_rast(
    data = filter(cell_coords, is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2),
    size = 1,
    color = "grey80" 
  ) +
  geom_point_rast(
    data = filter(cell_coords, !is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2, color = slingshot_pseudotime),
    size = 1
  ) +
  scale_color_gradientn(colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                   "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                   "#d73027", "#a50026"), 
                        na.value = "grey80") + 
  geom_path(
    data = lineage1_coords,
    aes(x = umap_1, y = umap_2),
    color = "black",  
    linewidth = 1.2      
  ) +
  labs(
    title = "Slingshot Trajectory on UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"   
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1,   
    legend.position = "top"
  )+
  guides(
    color = guide_colorbar(
      barwidth = unit(3, "cm"),  
      barheight = unit(0.5, "cm"), 
      title = "Pseudotime",
      title.theme = element_text(size = 11, face = "bold"),   
      title.vjust = 0.8, 
      label.theme = element_text(size = 8)   
    )
  ) #+
#SeuratExtend::theme_umap_arrows()

p_slingshot_umap
ggsave("slingshot_L4.pdf", p_slingshot_umap, width = 5, height = 5.5)

# L5 -----
sling_pseudotime <- sc$slingPseudotime_5
sling_meta <- data.frame(
  row.names = colnames(sc),
  slingshot_pseudotime = sling_pseudotime
)
identical(rownames(sling_meta), rownames(Cancer.OS.combined1@meta.data))
Cancer.OS.combined1 <- AddMetaData(Cancer.OS.combined1, metadata = sling_meta)
head(Cancer.OS.combined1@meta.data)

# curv info
sds <- SlingshotDataSet(sc)
lineage1_coords <- as.data.frame(sds@curves$Lineage5$s)
colnames(lineage1_coords) <- c("umap_1", "umap_2")

cell_coords <- as.data.frame(Embeddings(Cancer.OS.combined1, reduction = "umap"))
cell_coords$slingshot_pseudotime <- Cancer.OS.combined1$slingshot_pseudotime

p_slingshot_umap <- ggplot() +
  geom_point_rast(
    data = filter(cell_coords, is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2),
    size = 1,
    color = "grey80" 
  ) + 
  geom_point_rast(
    data = filter(cell_coords, !is.na(slingshot_pseudotime)),
    aes(x = umap_1, y = umap_2, color = slingshot_pseudotime),
    size = 1
  ) +
  scale_color_gradientn(colors = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                   "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                   "#d73027", "#a50026"), 
                        na.value = "grey80") + 
  geom_path(
    data = lineage1_coords,
    aes(x = umap_1, y = umap_2),
    color = "black",  
    linewidth = 1.2      
  ) +
  labs(
    title = "Slingshot Trajectory on UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Pseudotime"   
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio = 1,   
    legend.position = "top"
  )+
  guides(
    color = guide_colorbar(
      barwidth = unit(3, "cm"),  
      barheight = unit(0.5, "cm"), 
      title = "Pseudotime",
      title.theme = element_text(size = 11, face = "bold"),   
      title.vjust = 0.8, 
      label.theme = element_text(size = 8)   
    )
  ) #+
#SeuratExtend::theme_umap_arrows()

p_slingshot_umap
ggsave("slingshot_L5.pdf", p_slingshot_umap, width = 5, height = 5.5)

# Line segs -----
slingMST(sc)
mst <- slingMST(sc)
edges <- igraph::as_data_frame(mst)                                      
umap <- reducedDim(sc, "UMAP")
cluster <- colData(sc)$MP_group
umap_df <- data.frame(
  UMAP_1 = umap[,1],
  UMAP_2 = umap[,2],
  cluster = cluster
)

library(dplyr)
cluster_centers <- umap_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    UMAP_1 = mean(UMAP_1),
    UMAP_2 = mean(UMAP_2)
  )

edges_plot <- edges %>%
  dplyr::left_join(cluster_centers, by = c("from" = "cluster")) %>%
  dplyr::rename(x = UMAP_1, y = UMAP_2) %>%
  dplyr::left_join(cluster_centers, by = c("to" = "cluster")) %>%
  dplyr::rename(xend = UMAP_1, yend = UMAP_2)

pp2 <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = cluster), size = 0.5, 
             color = "grey80") + 
  geom_segment(
    data = edges_plot,
    aes(x = x, y = y, xend = xend, yend = yend#, 
        #linewidth = weight
    ),
    color = "black"#,
    #linewidth = 1
  )  +
  geom_point(
    data = cluster_centers,
    aes(UMAP_1, UMAP_2, colour = cluster),
    size = 5,
    shape = 16, 
    #linewidth = 1,
    color = c("#a97c50", "#2e3192", "#010101", "#972D15","#FBB4AE","#CCEBC5","#DECBE4","#FED9A6",
              "#7f3f98", "#B3CDE3")
  ) + 
  theme_classic() + 
  theme(aspect.ratio = 1)

pp2
ggsave("UMAP_graph_structure_260220.pdf", pp2, width = 5, height = 6)
