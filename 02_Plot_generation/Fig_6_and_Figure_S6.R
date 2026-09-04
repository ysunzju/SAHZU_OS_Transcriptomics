library(here) # project-root-relative paths; run scripts from repository root
##################### Fig 6 and extended data fig 6 #####################
# DC TF regulon (Related to Fig 6k) ------ 
rm(list = ls())
# disabled, run from repo root: setwd(here("data/PySCENIC_10_03_25/run1_myeloid_10_03_25/Vis_in_R"))
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

# Load scRNA object
DC <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/DC_251025/Myeloid_w_DC_idents_active_251211.rds"))
#Idents(DC)
New_Idents_DC <- subset(DC, cells = WhichCells(DC, idents = c("cDC2_CD83hi", "cDC2_CD83lo")))
#dim(New_Idents_DC)
#dim(cDC2)
#Idents(cDC2) <- Idents(DC)
New_Idents_DC$Cell_type_fine_harmony <- Idents(New_Idents_DC)

#sce <- readRDS(here("data/CytoTRACE_Myeloid_Cells_07_03_25/CytoTRACE2_score_for_Mye_TIL_08_03_25.rds"))
sce <- New_Idents_DC
rm(DC)
gc()

Idents(sce) <- "group"
sce <- RenameIdents(sce, "A" = "H_Pre", "B" = "H_Post", "C" = "L_Pre", "D" = "L_Post")
sce$group_anno <- Idents(sce)
#Idents(sce) <- "Cell_type_med_v1_07_03_25"
#sce <- subset(sce, idents = c("CD8T", "NK", "T", "B"))

# Load loom 
loom <- open_loom(here("data/PySCENIC_10_03_25/run1_myeloid_10_03_25/New_Vis_run_251214/cDC2_hi_vs_lo_251214/aucell.loom"))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4]

regulons <- regulonsToGeneLists(regulons_incidMat) 
class(regulons)

regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
head(regulonAUC)[1:3,1:3]

regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)
embeddings

#Idents(scrna) <- "Cell_type_fine_harmony"
#sce <- subset(scrna, idents = c("Macro_FOLR2"))
sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
dim(sub_regulonAUC) 
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
dim(sub_regulonAUC)

# Regulon Specificity Score, RSS
selectedResolution <- "celltype"
calcRSS <- function(AUC, cellAnnotation, cellTypes=NULL)
{
  if(any(is.na(cellAnnotation))) stop("NAs in annotation")
  if(any(class(AUC)=="aucellResults")) AUC <- getAUC(AUC)
  normAUC <- AUC/rowSums(AUC)
  if(is.null(cellTypes)) cellTypes <- unique(cellAnnotation)
  # 
  ctapply <- lapply
  if(require('BiocParallel')) ctapply <- bplapply
  
  rss <- ctapply(cellTypes, function(thisType)
    sapply(rownames(normAUC), function(thisRegulon)
    {
      pRegulon <- normAUC[thisRegulon,]
      pCellType <- as.numeric(cellAnnotation==thisType)
      pCellType <- pCellType/sum(pCellType)
      .calcRSS.oneRegulon(pRegulon, pCellType)
    })
  )
  rss <- do.call(cbind, rss)
  colnames(rss) <- cellTypes
  return(rss)
}

.calcRSS.oneRegulon <- function(pRegulon, pCellType)
{
  jsd <- calcJSD(pRegulon, pCellType)
  1 - sqrt(jsd)
}

calcJSD <- function(pRegulon, pCellType)
{
  (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
}

.H <- function(pVect){
  pVect <- pVect[pVect>0] # /sum(pVect) ??
  - sum(pVect * log2(pVect))
}

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),                
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),selectedResolution]) 
rss=na.omit(rss) 

plotRSS_heatmap <- plotRSS_heatmap <- function(rss, thr=NULL, row_names_gp=gpar(fontsize=5), order_rows=TRUE, cluster_rows=FALSE, name="RSS", verbose=TRUE, ...)
{
  if(is.null(thr)) thr <- signif(quantile(rss, p=.97),2)
  
  library(ComplexHeatmap)
  rssSubset <- rss[rowSums(rss > thr)>0,]
  rssSubset <- rssSubset[,colSums(rssSubset > thr)>0]
  
  if(verbose) message("Showing regulons and cell types with any RSS > ", thr, " (dim: ", nrow(rssSubset), "x", ncol(rssSubset),")")
  
  if(order_rows)
  {
    maxVal <- apply(rssSubset, 1, which.max)
    rss_ordered <- rssSubset[0,]
    for(i in 1:ncol(rssSubset))
    {
      tmp <- rssSubset[which(maxVal==i),,drop=F]
      tmp <- tmp[order(tmp[,i], decreasing=FALSE),,drop=F]
      rss_ordered <- rbind(rss_ordered, tmp)
    }
    rssSubset <- rss_ordered
    cluster_rows=FALSE
  }
  
  Heatmap(rssSubset, name=name, row_names_gp=row_names_gp, cluster_rows=cluster_rows, ...)
} 

plotRSS <- function(rss, labelsToDiscard=NULL, zThreshold=1,
                    cluster_columns=FALSE, order_rows=TRUE, thr=0.01, varName="cellType",
                    col.low="grey90", col.mid="darkolivegreen3", col.high="darkgreen",
                    revCol=FALSE, verbose=TRUE)
{
  varSize="RSS"
  varCol="Z"
  if(revCol) {
    varSize="Z"
    varCol="RSS"
  }
  
  rssNorm <- scale(rss) # scale the full matrix...
  rssNorm <- rssNorm[,which(!colnames(rssNorm) %in% labelsToDiscard)] # remove after calculating...
  rssNorm[rssNorm < 0] <- 0
  
  ## to get row order (easier...)
  rssSubset <- rssNorm
  if(!is.null(zThreshold)) rssSubset[rssSubset < zThreshold] <- 0
  tmp <- plotRSS_heatmap(rssSubset, thr=thr, cluster_columns=cluster_columns, order_rows=order_rows, verbose=verbose)
  rowOrder <- rev(tmp@row_names_param$labels)
  rm(tmp)
  
  
  ## Dotplot
  rss.df <- reshape2::melt(rss)
  head(rss.df)
  colnames(rss.df) <- c("Topic", varName, "RSS")
  rssNorm.df <- reshape2::melt(rssNorm)
  colnames(rssNorm.df) <- c("Topic", varName, "Z")
  rss.df <- base::merge(rss.df, rssNorm.df)
  
  rss.df <- rss.df[which(!rss.df[,varName] %in% labelsToDiscard),] # remove after calculating...
  if(nrow(rss.df)<2) stop("Insufficient rows left to plot RSS.")
  
  rss.df <- rss.df[which(rss.df$Topic %in% rowOrder),]
  rss.df[,"Topic"] <- factor(rss.df[,"Topic"], levels=rowOrder)
  p <- dotHeatmap(rss.df, 
                  var.x=varName, var.y="Topic", 
                  var.size=varSize, min.size=.5, max.size=5,
                  var.col=varCol, col.low=col.low, col.mid=col.mid, col.high=col.high)
  
  invisible(list(plot=p, df=rss.df, rowOrder=rowOrder))
}

dotHeatmap <- function (enrichmentDf,
                        var.x="Topic", var.y="ID", 
                        var.col="FC", col.low="dodgerblue", col.mid="floralwhite", col.high="brown1", 
                        var.size="p.adjust", min.size=1, max.size=8,
                        ...)
{
  require(data.table)
  require(ggplot2)
  
  colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
  p <- ggplot(data=enrichmentDf, mapping=aes_string(x=var.x, y=var.y)) + 
    geom_point(mapping=aes_string(size=var.size, color=var.col)) +
    scale_radius(range=c(min.size, max.size)) +
    scale_colour_gradientn(colors=colorPal(10)) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y=element_blank(), 
          axis.text.x=element_text(angle=90, hjust=1),
          ...)
  return(p)
}

##### rssPlot #####
rss <- rss[, c(2,1)]
rssPlot <- plotRSS(rss,                   
                   labelsToDiscard = NULL,                    
                   zThreshold = 1,                    
                   cluster_columns = FALSE,                    
                   order_rows = T,                    
                   thr = 0.01,                    
                   varName = "cellType",                   
                   col.low = '#e9f4f6',                     
                   col.mid = '#8f9fab',                     
                   col.high= '#1e2235',                   
                   revCol = T,                   
                   verbose = TRUE)

p1 <- rssPlot$plot +
  coord_fixed()
p1

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_rssPlot_regulon_cDC2_hi_v_lo.pdf"), 
    width = 3, # The width of the plot in inches
    height = 7) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()

# Cell_chat plots (Related to Extended data fig 6g) ------
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Correlation_and_crosstalk_analysis_18_07_25/CellChat_251218"))

library(CellChat)
library(patchwork)

c_levels <- c("NK_KIT+", "NK_FCGR3A+", "NK_KLRC1+", "Tn_CCR7+", "Tmem_IL7R+", "Tmem_GPR183+FOSB+", 
              "Trm_ZNF683+", "Treg", "CD8_Teff_IFNG+", "CD8_Teff_GZMK+", "Tex_CXCL13+", "T_Ki67+", "T_Mtgenes_hi", "gdT", "CD8_MAIT",
              "Pro_B", "BCells", "PlasmaCells","Mono_FCN1+", "Macro_CCL4+", "Macro_APOE+", "Macro_ISGs+", "Macro_MTs+", 
              "Macro_Ki67+", "Promyelocytes", "Neu", "cDC1_CLEC9A+", "cDC2_CD83hi", "cDC2_CD83lo", "Mature_DC", "MastCells", "pDCs")
no_show_b_cells <- c("NK_KIT+", "NK_FCGR3A+", "NK_KLRC1+", "Tn_CCR7+", "Tmem_IL7R+", "Tmem_GPR183+FOSB+", 
                     "Trm_ZNF683+", "Treg", "CD8_Teff_IFNG+", "CD8_Teff_GZMK+", "Tex_CXCL13+", "T_Ki67+", "T_Mtgenes_hi", "gdT", "CD8_MAIT",
                     "Mono_FCN1+", "Macro_CCL4+", "Macro_APOE+", "Macro_ISGs+", "Macro_MTs+", 
                     "Macro_Ki67+", "Promyelocytes", "Neu", "cDC1_CLEC9A+", "cDC2_CD83hi", "cDC2_CD83lo", "Mature_DC", "MastCells")

# Differential heatmap of interaction strengths between NR Post and R Post ----- 
## Load CellChat object of each dataset and merge them together ---- 
cellchat.NR.post <- readRDS(here("data/Correlation_and_crosstalk_analysis_18_07_25/CellChat_251218/CellChat_all_immune_minor_celltype_NR_Post_251221.rds"))
cellchat.R.post <- readRDS(here("data/Correlation_and_crosstalk_analysis_18_07_25/CellChat_251218/CellChat_all_immune_minor_celltype_H_Post_251221.rds")
)
object.list <- list(NR.post = cellchat.NR.post, R.post = cellchat.R.post)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
cellchat@idents$joint <- factor(cellchat@idents$joint, levels = c_levels)
cellchat@idents$R.post <- factor(cellchat@idents$R.post, levels = c_levels)
cellchat@idents$NR.post <- factor(cellchat@idents$NR.post, levels = c_levels)
levels(cellchat@idents$joint)

gg1 <- netVisual_heatmap(cellchat, 
                         sources.use = c_levels,
                         targets.use = c_levels, 
                         cluster.rows = T, 
                         cluster.cols = T, 
                         color.heatmap = c("#50859f","#d66692"))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",
                         sources.use = c_levels,
                         targets.use = c_levels, 
                         cluster.rows = T, 
                         cluster.cols = T, 
                         color.heatmap = c("#50859f","#d66692")) 
ggg1 <- gg1 + gg2 
ggg1

#Manual saves 4 x 9 
#All_interact_n_and_strengths_NR_v_R_Post_0208.pdf

#time <- gsub(" ", "_", Sys.time()) 
#time <- gsub("-", "_", time) 
#time <- gsub(":", "_", time) 
#pdf(file = paste0(time, "All_interact_n_and_strengths_NR_v_R_Post_0208.pdf"), 
#    width = 4, # The width of the plot in inches
#    height = 9) # The height of the plot in inches
#print(ggg1, newpage = FALSE)
