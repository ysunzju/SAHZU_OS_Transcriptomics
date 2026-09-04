library(here) # project-root-relative paths; run scripts from repository root
##################### Bulk RNA-seq data preprocessing and analysis ##################### 
##################### Section I - Batch effect removal ##################### 
rm(list = ls())
# disabled, run from repo root: setwd(here("data/New_bulk_ana_250907/New_Survival_250907"))

# Load bulk-seq data (fpkm matrices) ----- 
## 144 patients (Batch 1) ----- 
fpkm_mat_complete <- read.csv(here("data/25_04_25_New_OS_Bulk/OS_bulk_RNA_seq/P22051002_result_144_complete/3.Quant/gene_fpkm.txt"), 
                              sep = "\t", header = T)
fpkm_mat_144 <- fpkm_mat_complete
#data = data[data$gene_biotype == "protein_coding",]
fpkm_mat_144 = fpkm_mat_144[!duplicated(fpkm_mat_144$gene_name),]
rownames(fpkm_mat_144) = fpkm_mat_144$gene_name
colnames(fpkm_mat_144)
fpkm_mat_144 <- fpkm_mat_144[, c(146, 2:145)]
colnames(fpkm_mat_144)[1] <- "GeneSymbol"

## 38 patients (Batch 2) ----- 
fpkm_mat_38 <- readxl::read_xlsx(here("data/EcoTyper_run_16_06_25/39_Pt_OS_fpkm.xlsx"))
colnames(fpkm_mat_38)
fpkm_mat_38 <- fpkm_mat_38[, c(1, 4:42)] 

# Clinical meta data loading ----- 
#group <- data.frame(readxl::read_xlsx(here("data/EcoTyper_run_16_06_25/Meta_182_v1.csv.xlsx")))
group <- read.csv(here("data/Re-clustering_fine_cell_types_21_05_25/Metadata_260802/updated_bulk_metadta_260802.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")
colnames(group)[3] <- "Sample"

meta38 <- group[which(group$Group == "Seq_Batch_2"), ]
for (i in colnames(fpkm_mat_38)) {
  if (i %in% meta38$Name) {
    colnames(fpkm_mat_38)[which(colnames(fpkm_mat_38) == i)] <-  meta38$Sample[which(meta38$Name == i)]
  }
}

# Merge ----- 
fpkm_merge <- merge(fpkm_mat_144, fpkm_mat_38, by = "GeneSymbol")

# PCA for combined dataset, not batch-effect corrected ----- 
fpkm_merge = fpkm_merge[!duplicated(fpkm_merge$GeneSymbol),]
mat_for_ana <- fpkm_merge
rownames(mat_for_ana)  <- mat_for_ana$GeneSymbol
mat_for_ana <- mat_for_ana[, -1]
mat_for_ana <- t(mat_for_ana)
#data.pca <- prcomp(mat_for_ana)

# PCA {FactoMineR}
library(ggplot2)
library(ggrepel)
library(FactoMineR)

gene.pca <- PCA(mat_for_ana, ncp = 2, scale.unit = TRUE, graph = FALSE)

# Extract sample embeddings 
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
pca_sample$Sample=row.names(pca_sample)
# Extract PC1 and PC2 
pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2], 2)
group_anno <- data.frame(matrix(ncol = 2, nrow = nrow(group)))
names(group_anno) <- c("Sample", "Group")
group_anno$Sample <- group$Sample
group_anno$Group <- group$Group
group_anno$Treatment <- group$Treatment

pca_sample <- merge(pca_sample, group_anno, by="Sample")
head(pca_sample)

# Plot 
library(ggplot2)
p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 2) + # "Group" refers to sequencing batches
  scale_color_manual(values = c('orange', 'purple', "grey")) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  

# View - This would demonstrate significant batch effect 
p + stat_ellipse(aes(color = Group), level = 0.95, show.legend = FALSE)

## Remove batch effect -----
#BiocManager::install("sva")
#install.packages("sva")
library(sva)
mat_for_ana <- fpkm_merge 
rownames(mat_for_ana)  <- mat_for_ana$GeneSymbol 
mat_for_ana <- mat_for_ana[, -1] 
#identical(group$ID, colnames(mat_for_ana)) # FALSE 
group_anno$ID <- factor(group_anno$Sample, levels = colnames(mat_for_ana)) 
#order(group$ID)
#class(as.character(group[order(group$ID), "ID"]))
#class(colnames(mat_for_ana))
#length(as.character(group[order(group$ID), "ID"]))
#length(colnames(mat_for_ana))
identical(as.character(group_anno[order(group_anno$ID), "ID"]), colnames(mat_for_ana)) # TRUE
#group[order(group$ID), "ID"] %in% intersect(as.character(group[order(group$ID), "ID"]), colnames(mat_for_ana)) 
#colnames(mat_for_ana)[!colnames(mat_for_ana) %in% intersect(as.character(group[order(group$ID), "ID"]), colnames(mat_for_ana))]
group_anno <- group_anno[order(group_anno$ID), ]
rownames(group_anno) <- group_anno$ID
identical(rownames(group_anno), colnames(mat_for_ana)) # TRUE

#install.packages("tinyarray")
library(tinyarray)
draw_pca(exp = mat_for_ana, group_list = factor(group_anno$Group)) # Same as before 

library(sva)
expr_combat <- ComBat(dat = mat_for_ana, batch = group_anno$Group)
## Found5batches
## Adjusting for0covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
## Finding parametric adjustments
## Adjusting the Data

expr_combat <- as.data.frame(expr_combat)

# Mitigation of batch effect 
draw_pca(exp = expr_combat, group_list = factor(group_anno$Group))

##################### Section II - Survival (Related to Figures 3I and 7O)##################### 
library(TCGAbiolinks) 
library(dplyr) 
library(SummarizedExperiment) 
library(GSVA) 
library(limma) 
library(survival) 
library(survminer) 
library(Seurat) 

expr_combat <- read.table(here("data/EcoTyper_run_16_06_25/OS_182_post_batch_rem_260731.txt"))
group <- read.csv(here("data/Re-clustering_fine_cell_types_21_05_25/Metadata_260802/updated_bulk_metadta_260802.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")

## Post with available OS days only  ----- 
## For analysis based on pre-chemotherapy data, subsetting data to "Pre"
select_group <- group[group$Treatment %in% "Post", ] 
select_group$OS <- as.numeric(select_group$OS) 
select_group <- select_group[!is.na(select_group$OS), ] 
select_group$Death <- as.numeric(select_group$Death)    

# Keep matrix conforming to metadata
processed_mat_df <- expr_combat
dat_temp <- processed_mat_df[, colnames(processed_mat_df) %in% select_group$Sample]

# GeneSet construction (Change as need be, here generation) 
### S Phase ----- 
MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
genelist_down_manual <- MP_list[["MP_9"]] # S Phase 
#genelist_down_manual <- MP_list[["MP_1"]] # Ossification/mineralisation 
#genelist_down_manual <- MP_list[["MP_2"]] # ECM 

# Formatting geneset for ssGSEA
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes") # Pseudonames
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssGSEA 
dat <- as.matrix(dat_temp)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)

# Construct survival matrix 
c1 <- as.numeric(which(colnames(select_group) == "Sample")) 
c2 <- as.numeric(which(colnames(select_group) == "Name")) 
c3 <- as.numeric(which(colnames(select_group) == "Death")) 
c4 <- as.numeric(which(colnames(select_group) == "OS")) 
Surv_matrix = select_group[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character" 
class(Surv_matrix[, 2]) <- "character" 
colnames(Surv_matrix)[3] <- "vital_status" 
colnames(Surv_matrix)[4] <- "months_to_death" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "1")] <- "2" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "0")] <- "1" 
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

class(Surv_matrix$vital_status) <- "numeric" 
class(Surv_matrix$months_to_death) <- "numeric"  
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5
class(Surv_matrix$months_to_death) <- "numeric"

# Merging ssGSEA results with survival data matrix 
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores (NOT run)
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.5)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

# For best cutoff 
res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

# Plot generation
plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "S",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#7f3f98", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("IN_house_S_Phase", '_Surv_plot.pdf'),   
    width = 5, 
    height = 5) 
print(plot, newpage = FALSE)
dev.off()

# Saving the merged bulk expression matrix ffile 
#saveRDS(fpkm_merge, file = "Bulk_183_exp_mat_post_batch_eff_rem_251119.rds")

### OS/Min ----- 
MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
#genelist_down_manual <- MP_list[["MP_9"]] # S Phase 
genelist_down_manual <- MP_list[["MP_1"]] # Ossification/mineralisation 
#genelist_down_manual <- MP_list[["MP_2"]] # ECM 

# Formatting geneset for ssGSEA
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes") # Pseudonames
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssGSEA 
dat <- as.matrix(dat_temp)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)

# Construct survival matrix 
c1 <- as.numeric(which(colnames(select_group) == "Sample")) 
c2 <- as.numeric(which(colnames(select_group) == "Name")) 
c3 <- as.numeric(which(colnames(select_group) == "Death")) 
c4 <- as.numeric(which(colnames(select_group) == "OS")) 
Surv_matrix = select_group[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character" 
class(Surv_matrix[, 2]) <- "character" 
colnames(Surv_matrix)[3] <- "vital_status" 
colnames(Surv_matrix)[4] <- "months_to_death" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "1")] <- "2" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "0")] <- "1" 
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

class(Surv_matrix$vital_status) <- "numeric" 
class(Surv_matrix$months_to_death) <- "numeric"  
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5
class(Surv_matrix$months_to_death) <- "numeric"

# Merging ssGSEA results with survival data matrix 
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores (NOT run)
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.5)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

# For best cutoff 
res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

# Plot generation
plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "OSMIN",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#a97c51", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("IN_house_OS_MIN", '_Surv_plot.pdf'),   
    width = 5, 
    height = 5) 
print(plot, newpage = FALSE)
dev.off()

### ECM -----
MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
#genelist_down_manual <- MP_list[["MP_9"]] # S Phase 
#genelist_down_manual <- MP_list[["MP_1"]] # Ossification/mineralisation 
genelist_down_manual <- MP_list[["MP_2"]] # ECM 

# Formatting geneset for ssGSEA
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes") # Pseudonames
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssGSEA 
dat <- as.matrix(dat_temp)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)

# Construct survival matrix 
c1 <- as.numeric(which(colnames(select_group) == "Sample")) 
c2 <- as.numeric(which(colnames(select_group) == "Name")) 
c3 <- as.numeric(which(colnames(select_group) == "Death")) 
c4 <- as.numeric(which(colnames(select_group) == "OS")) 
Surv_matrix = select_group[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character" 
class(Surv_matrix[, 2]) <- "character" 
colnames(Surv_matrix)[3] <- "vital_status" 
colnames(Surv_matrix)[4] <- "months_to_death" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "1")] <- "2" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "0")] <- "1" 
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

class(Surv_matrix$vital_status) <- "numeric" 
class(Surv_matrix$months_to_death) <- "numeric"  
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5
class(Surv_matrix$months_to_death) <- "numeric"

# Merging ssGSEA results with survival data matrix 
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores (NOT run)
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.5)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

# For best cutoff 
res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

# Plot generation
plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "ECM",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#2e368f", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("IN_house_ECM", '_Surv_plot.pdf'),   
    width = 5, 
    height = 5) 
print(plot, newpage = FALSE)
dev.off()

### cDC2 ----- 
#MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
#genelist_down_manual <- MP_list[["MP_9"]] # S Phase 
#genelist_down_manual <- MP_list[["MP_1"]] # Ossification/mineralisation 
genelist_down_manual <- c("FCER1A","CD1C","CLEC10A","HLA-DQA1","HLA-DQB1","JAML","NR4A3","GPR183",
                          "CD83","LST1","HLA-DMB","HLA-DQA2","LYZ","BCL2A1","IL1B","C1orf162","FCGR2B",
                          "HLA-DMA","INSIG1","PLEK","PHACTR1","PLAUR","NAMPT","AIF1","HLA-DRB6","NR4A2",
                          "LCP1","DUSP2","HLA-DPA1","FCER1G","LSP1","C15orf48","MS4A6A","MNDA","HLA-DPB1",
                          "NAPSB","REL","TYROBP","CDKN1A","RGS1","CXCL8","HLA-DRB1","PTPRE","LGALS2",
                          "CD86","FGL2","CXCR4","CYTIP","G0S2","RGS2")

# Formatting geneset for ssGSEA
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes") # Pseudonames
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssGSEA 
dat <- as.matrix(dat_temp)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)

# Construct survival matrix 
c1 <- as.numeric(which(colnames(select_group) == "Sample")) 
c2 <- as.numeric(which(colnames(select_group) == "Name")) 
c3 <- as.numeric(which(colnames(select_group) == "Death")) 
c4 <- as.numeric(which(colnames(select_group) == "OS")) 
Surv_matrix = select_group[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character" 
class(Surv_matrix[, 2]) <- "character" 
colnames(Surv_matrix)[3] <- "vital_status" 
colnames(Surv_matrix)[4] <- "months_to_death" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "1")] <- "2" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "0")] <- "1" 
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

class(Surv_matrix$vital_status) <- "numeric" 
class(Surv_matrix$months_to_death) <- "numeric"  
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5
class(Surv_matrix$months_to_death) <- "numeric"

# Merging ssGSEA results with survival data matrix 
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores (NOT run)
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.5)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

# For best cutoff 
res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

# Plot generation
plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "cDC2",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#5c696e", "darkgrey")) # "#00214F" for Th FOSB+ 

plot #
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("IN_house_cDC2", '_post_surv_plot.pdf'),   
    width = 5, 
    height = 5) 
print(plot, newpage = FALSE)
dev.off()
## Pre with available OS days only  ----- 
## For analysis based on pre-chemotherapy data, subsetting data to "Pre"
select_group <- group[group$Treatment %in% "Pre", ] 
select_group$OS <- as.numeric(select_group$OS) 
select_group <- select_group[!is.na(select_group$OS), ] 
select_group$Death <- as.numeric(select_group$Death)    

# Keep matrix conforming to metadata
processed_mat_df <- expr_combat
dat_temp <- processed_mat_df[, colnames(processed_mat_df) %in% select_group$Sample]

# GeneSet construction (Change as need be, here generation) 
### S Phase ----- 
MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
genelist_down_manual <- MP_list[["MP_9"]] # S Phase 
#genelist_down_manual <- MP_list[["MP_1"]] # Ossification/mineralisation 
#genelist_down_manual <- MP_list[["MP_2"]] # ECM 

# Formatting geneset for ssGSEA
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes") # Pseudonames
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssGSEA 
dat <- as.matrix(dat_temp)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)

# Construct survival matrix 
c1 <- as.numeric(which(colnames(select_group) == "Sample")) 
c2 <- as.numeric(which(colnames(select_group) == "Name")) 
c3 <- as.numeric(which(colnames(select_group) == "Death")) 
c4 <- as.numeric(which(colnames(select_group) == "OS")) 
Surv_matrix = select_group[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character" 
class(Surv_matrix[, 2]) <- "character" 
colnames(Surv_matrix)[3] <- "vital_status" 
colnames(Surv_matrix)[4] <- "months_to_death" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "1")] <- "2" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "0")] <- "1" 
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

class(Surv_matrix$vital_status) <- "numeric" 
class(Surv_matrix$months_to_death) <- "numeric"  
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5
class(Surv_matrix$months_to_death) <- "numeric"

# Merging ssGSEA results with survival data matrix 
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores (NOT run)
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.5)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

# For best cutoff 
res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

# Plot generation
plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "S",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   #axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#7f3f98", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("IN_house_S_Phase_Pre", '_Surv_plot.pdf'),   
    width = 5, 
    height = 5) 
print(plot, newpage = FALSE)
dev.off()

# Saving the merged bulk expression matrix ffile 
#saveRDS(fpkm_merge, file = "Bulk_183_exp_mat_post_batch_eff_rem_251119.rds")

### OS/Min ----- 
MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
#genelist_down_manual <- MP_list[["MP_9"]] # S Phase 
genelist_down_manual <- MP_list[["MP_1"]] # Ossification/mineralisation 
#genelist_down_manual <- MP_list[["MP_2"]] # ECM 

# Formatting geneset for ssGSEA
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes") # Pseudonames
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssGSEA 
dat <- as.matrix(dat_temp)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)

# Construct survival matrix 
c1 <- as.numeric(which(colnames(select_group) == "Sample")) 
c2 <- as.numeric(which(colnames(select_group) == "Name")) 
c3 <- as.numeric(which(colnames(select_group) == "Death")) 
c4 <- as.numeric(which(colnames(select_group) == "OS")) 
Surv_matrix = select_group[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character" 
class(Surv_matrix[, 2]) <- "character" 
colnames(Surv_matrix)[3] <- "vital_status" 
colnames(Surv_matrix)[4] <- "months_to_death" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "1")] <- "2" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "0")] <- "1" 
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

class(Surv_matrix$vital_status) <- "numeric" 
class(Surv_matrix$months_to_death) <- "numeric"  
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5
class(Surv_matrix$months_to_death) <- "numeric"

# Merging ssGSEA results with survival data matrix 
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores (NOT run)
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.5)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

# For best cutoff 
res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

# Plot generation
plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "OSMIN",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#a97c51", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("IN_house_OS_MIN_Pre", '_Surv_plot.pdf'),   
    width = 5, 
    height = 5) 
print(plot, newpage = FALSE)
dev.off()

### ECM -----
MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
#genelist_down_manual <- MP_list[["MP_9"]] # S Phase 
#genelist_down_manual <- MP_list[["MP_1"]] # Ossification/mineralisation 
genelist_down_manual <- MP_list[["MP_2"]] # ECM 

# Formatting geneset for ssGSEA
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes") # Pseudonames
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssGSEA 
dat <- as.matrix(dat_temp)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)

# Construct survival matrix 
c1 <- as.numeric(which(colnames(select_group) == "Sample")) 
c2 <- as.numeric(which(colnames(select_group) == "Name")) 
c3 <- as.numeric(which(colnames(select_group) == "Death")) 
c4 <- as.numeric(which(colnames(select_group) == "OS")) 
Surv_matrix = select_group[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character" 
class(Surv_matrix[, 2]) <- "character" 
colnames(Surv_matrix)[3] <- "vital_status" 
colnames(Surv_matrix)[4] <- "months_to_death" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "1")] <- "2" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "0")] <- "1" 
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

class(Surv_matrix$vital_status) <- "numeric" 
class(Surv_matrix$months_to_death) <- "numeric"  
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5
class(Surv_matrix$months_to_death) <- "numeric"

# Merging ssGSEA results with survival data matrix 
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores (NOT run)
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.5)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

# For best cutoff 
res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

# Plot generation
plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "ECM",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#2e368f", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("IN_house_ECM_Pre", '_Surv_plot.pdf'),   
    width = 5, 
    height = 5) 
print(plot, newpage = FALSE)
dev.off()

### cDC2 ----- 
#MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
#genelist_down_manual <- MP_list[["MP_9"]] # S Phase 
#genelist_down_manual <- MP_list[["MP_1"]] # Ossification/mineralisation 
genelist_down_manual <- c("FCER1A","CD1C","CLEC10A","HLA-DQA1","HLA-DQB1","JAML","NR4A3","GPR183",
                          "CD83","LST1","HLA-DMB","HLA-DQA2","LYZ","BCL2A1","IL1B","C1orf162","FCGR2B",
                          "HLA-DMA","INSIG1","PLEK","PHACTR1","PLAUR","NAMPT","AIF1","HLA-DRB6","NR4A2",
                          "LCP1","DUSP2","HLA-DPA1","FCER1G","LSP1","C15orf48","MS4A6A","MNDA","HLA-DPB1",
                          "NAPSB","REL","TYROBP","CDKN1A","RGS1","CXCL8","HLA-DRB1","PTPRE","LGALS2",
                          "CD86","FGL2","CXCR4","CYTIP","G0S2","RGS2")

# Formatting geneset for ssGSEA
gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes") # Pseudonames
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssGSEA 
dat <- as.matrix(dat_temp)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)

# Construct survival matrix 
c1 <- as.numeric(which(colnames(select_group) == "Sample")) 
c2 <- as.numeric(which(colnames(select_group) == "Name")) 
c3 <- as.numeric(which(colnames(select_group) == "Death")) 
c4 <- as.numeric(which(colnames(select_group) == "OS")) 
Surv_matrix = select_group[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character" 
class(Surv_matrix[, 2]) <- "character" 
colnames(Surv_matrix)[3] <- "vital_status" 
colnames(Surv_matrix)[4] <- "months_to_death" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "1")] <- "2" 
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "0")] <- "1" 
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

class(Surv_matrix$vital_status) <- "numeric" 
class(Surv_matrix$months_to_death) <- "numeric"  
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5
class(Surv_matrix$months_to_death) <- "numeric"

# Merging ssGSEA results with survival data matrix 
ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores (NOT run)
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.5)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

# For best cutoff 
res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

# Plot generation
plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "cDC2",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#5c696e", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("IN_house_cDC2", '_pre_Surv_plot.pdf'),   
    width = 5, 
    height = 5) 
print(plot, newpage = FALSE)
dev.off()

## Surv w/ TARGET-OS ----- 
### S phase ----- 
i <- "TARGET-OS"
load(paste0(here("data/TCGA_TARGET_OS_Surv_23_04_25/03_09_Other_cancers_bulk/"), i, "_mRNA.Rdata"))
# Subset mRNA data to protein-coding 
rownames(data) = rowData(data)$gene_name
data = data[!duplicated(rownames(data)),]
data = data[rowData(data)$gene_type=="protein_coding",]

# Data selection
unique(colData(data)$tissue_type)
data = data[, colData(data)$tissue_type=="Tumor"] # Keep only tumor
if (length(which(is.na(colData(data)$vital_status))) != 0) {
  print("NAs removed")
  data <- data[, -which(is.na(colData(data)$vital_status) == TRUE)] # Remove NA 
} else {
  print("No NAs to be removed")
}
data = data[, !colData(data)$vital_status == "Not Reported"] # Keep only with reported vital_status

# Obtain tpm data
tpm <- assay(data, "tpm_unstrand")
logTPM <- log2(tpm+1) #logNorm

MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
# GeneSet construction (Change as need be)
genelist_down_manual <- MP_list[["MP_9"]]
#genelist_down_manual <- MP_list[["MP_1"]]
#genelist_down_manual <- MP_list[["MP_2"]]

gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v2 <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v1 <- rbind(gene.set_v1, gene.set_v2)
#rm(gene.set_v2)
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssgsea
dat <- as.matrix(logTPM)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)
# Construct survival matrix 
Surv_matrix <- as.matrix(colData(data))
Surv_matrix <- as.data.frame(Surv_matrix)
c1 <- as.numeric(which(colnames(Surv_matrix) == "barcode"))
c2 <- as.numeric(which(colnames(Surv_matrix) == "patient"))
c3 <- as.numeric(which(colnames(Surv_matrix) == "vital_status"))
c4 <- as.numeric(which(colnames(Surv_matrix) == "days_to_death"))
Surv_matrix = Surv_matrix[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character"
class(Surv_matrix[, 2]) <- "character"

Surv_matrix$vital_status[which(Surv_matrix$vital_status == "Alive")] <- "1"
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "Dead")] <- "2"
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

colData(data)$vital_status

class(Surv_matrix$vital_status) <- "numeric"
class(Surv_matrix$days_to_death) <- "numeric" 
Surv_matrix$months_to_death <- Surv_matrix$days_to_death / 30
dd1 <- max(Surv_matrix$days_to_death, na.rm = T) + 100
Surv_matrix$days_to_death[is.na(Surv_matrix$days_to_death)] <- dd1
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5

class(Surv_matrix$days_to_death) <- "numeric"
class(Surv_matrix$months_to_death) <- "numeric"

ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores 
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.75)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(days_to_death, vital_status) ~ Iron_score_status, data = Surv_matrix)
#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") 
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "S",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   #axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#7f3f98", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("S_Phase_TARGET_OS", '_Surv_plot_16_12.pdf'),   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(plot, newpage = FALSE)
dev.off()

### Osteo/Min  ----- 
i <- "TARGET-OS"
load(paste0(here("data/TCGA_TARGET_OS_Surv_23_04_25/03_09_Other_cancers_bulk/"), i, "_mRNA.Rdata"))
# Subset mRNA data to protein-coding 
rownames(data) = rowData(data)$gene_name
data = data[!duplicated(rownames(data)),]
data = data[rowData(data)$gene_type=="protein_coding",]

# Data selection
unique(colData(data)$tissue_type)
data = data[, colData(data)$tissue_type=="Tumor"] # Keep only tumor
if (length(which(is.na(colData(data)$vital_status))) != 0) {
  print("NAs removed")
  data <- data[, -which(is.na(colData(data)$vital_status) == TRUE)] # Remove NA 
} else {
  print("No NAs to be removed")
}
data = data[, !colData(data)$vital_status == "Not Reported"] # Keep only with reported vital_status

# Obtain tpm data
tpm <- assay(data, "tpm_unstrand")
logTPM <- log2(tpm+1) #logNorm

MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
# GeneSet construction (Change as need be)
#genelist_down_manual <- MP_list[["MP_9"]]
genelist_down_manual <- MP_list[["MP_1"]]
#genelist_down_manual <- MP_list[["MP_2"]]

gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v2 <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v1 <- rbind(gene.set_v1, gene.set_v2)
#rm(gene.set_v2)
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssgsea
dat <- as.matrix(logTPM)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)
# Construct survival matrix 
Surv_matrix <- as.matrix(colData(data))
Surv_matrix <- as.data.frame(Surv_matrix)
c1 <- as.numeric(which(colnames(Surv_matrix) == "barcode"))
c2 <- as.numeric(which(colnames(Surv_matrix) == "patient"))
c3 <- as.numeric(which(colnames(Surv_matrix) == "vital_status"))
c4 <- as.numeric(which(colnames(Surv_matrix) == "days_to_death"))
Surv_matrix = Surv_matrix[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character"
class(Surv_matrix[, 2]) <- "character"

Surv_matrix$vital_status[which(Surv_matrix$vital_status == "Alive")] <- "1"
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "Dead")] <- "2"
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

colData(data)$vital_status

class(Surv_matrix$vital_status) <- "numeric"
class(Surv_matrix$days_to_death) <- "numeric" 
Surv_matrix$months_to_death <- Surv_matrix$days_to_death / 30
dd1 <- max(Surv_matrix$days_to_death, na.rm = T) + 100
Surv_matrix$days_to_death[is.na(Surv_matrix$days_to_death)] <- dd1
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5

class(Surv_matrix$days_to_death) <- "numeric"
class(Surv_matrix$months_to_death) <- "numeric"

ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores 
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.75)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(days_to_death, vital_status) ~ Iron_score_status, data = Surv_matrix)
#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") # 找这3个变量的最佳切点
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "S",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#a97d52", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("Osteo_TARGET_OS", '_Surv_plot_16_12.pdf'),   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(plot, newpage = FALSE)
dev.off()

### ECM  ----- 
i <- "TARGET-OS"
load(paste0(here("data/TCGA_TARGET_OS_Surv_23_04_25/03_09_Other_cancers_bulk/"), i, "_mRNA.Rdata"))
# Subset mRNA data to protein-coding 
rownames(data) = rowData(data)$gene_name
data = data[!duplicated(rownames(data)),]
data = data[rowData(data)$gene_type=="protein_coding",]

# Data selection
unique(colData(data)$tissue_type)
data = data[, colData(data)$tissue_type=="Tumor"] # Keep only tumor
if (length(which(is.na(colData(data)$vital_status))) != 0) {
  print("NAs removed")
  data <- data[, -which(is.na(colData(data)$vital_status) == TRUE)] # Remove NA 
} else {
  print("No NAs to be removed")
}
data = data[, !colData(data)$vital_status == "Not Reported"] # Keep only with reported vital_status

# Obtain tpm data
tpm <- assay(data, "tpm_unstrand")
logTPM <- log2(tpm+1) #logNorm

MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
# GeneSet construction (Change as need be)
#genelist_down_manual <- MP_list[["MP_9"]]
#genelist_down_manual <- MP_list[["MP_1"]]
genelist_down_manual <- MP_list[["MP_2"]]

gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v2 <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v1 <- rbind(gene.set_v1, gene.set_v2)
#rm(gene.set_v2)
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssgsea
dat <- as.matrix(logTPM)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)
# Construct survival matrix 
Surv_matrix <- as.matrix(colData(data))
Surv_matrix <- as.data.frame(Surv_matrix)
c1 <- as.numeric(which(colnames(Surv_matrix) == "barcode"))
c2 <- as.numeric(which(colnames(Surv_matrix) == "patient"))
c3 <- as.numeric(which(colnames(Surv_matrix) == "vital_status"))
c4 <- as.numeric(which(colnames(Surv_matrix) == "days_to_death"))
Surv_matrix = Surv_matrix[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character"
class(Surv_matrix[, 2]) <- "character"

Surv_matrix$vital_status[which(Surv_matrix$vital_status == "Alive")] <- "1"
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "Dead")] <- "2"
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

colData(data)$vital_status

class(Surv_matrix$vital_status) <- "numeric"
class(Surv_matrix$days_to_death) <- "numeric" 
Surv_matrix$months_to_death <- Surv_matrix$days_to_death / 30
dd1 <- max(Surv_matrix$days_to_death, na.rm = T) + 100
Surv_matrix$days_to_death[is.na(Surv_matrix$days_to_death)] <- dd1
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5

class(Surv_matrix$days_to_death) <- "numeric"
class(Surv_matrix$months_to_death) <- "numeric"

ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores 
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.75)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(days_to_death, vital_status) ~ Iron_score_status, data = Surv_matrix)
#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") # 找这3个变量的最佳切点
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "S",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#2e398e", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("ECM_TARGET_OS", '_Surv_plot_16_12.pdf'),   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(plot, newpage = FALSE)
dev.off()

### cDC2  ----- 
i <- "TARGET-OS"
load(paste0(here("data/TCGA_TARGET_OS_Surv_23_04_25/03_09_Other_cancers_bulk/"), i, "_mRNA.Rdata"))
# Subset mRNA data to protein-coding 
rownames(data) = rowData(data)$gene_name
data = data[!duplicated(rownames(data)),]
data = data[rowData(data)$gene_type=="protein_coding",]

# Data selection
unique(colData(data)$tissue_type)
data = data[, colData(data)$tissue_type=="Tumor"] # Keep only tumor
if (length(which(is.na(colData(data)$vital_status))) != 0) {
  print("NAs removed")
  data <- data[, -which(is.na(colData(data)$vital_status) == TRUE)] # Remove NA 
} else {
  print("No NAs to be removed")
}
data = data[, !colData(data)$vital_status == "Not Reported"] # Keep only with reported vital_status

# Obtain tpm data
tpm <- assay(data, "tpm_unstrand")
logTPM <- log2(tpm+1) #logNorm

#MP_list <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))
# GeneSet construction (Change as need be)
#genelist_down_manual <- MP_list[["MP_9"]]
#genelist_down_manual <- MP_list[["MP_1"]]
genelist_down_manual <- c("FCER1A","CD1C","CLEC10A","HLA-DQA1","HLA-DQB1","JAML","NR4A3","GPR183",
                          "CD83","LST1","HLA-DMB","HLA-DQA2","LYZ","BCL2A1","IL1B","C1orf162","FCGR2B",
                          "HLA-DMA","INSIG1","PLEK","PHACTR1","PLAUR","NAMPT","AIF1","HLA-DRB6","NR4A2",
                          "LCP1","DUSP2","HLA-DPA1","FCER1G","LSP1","C15orf48","MS4A6A","MNDA","HLA-DPB1",
                          "NAPSB","REL","TYROBP","CDKN1A","RGS1","CXCL8","HLA-DRB1","PTPRE","LGALS2",
                          "CD86","FGL2","CXCR4","CYTIP","G0S2","RGS2")

gene.set <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v2 <- data.frame("Metagene" = genelist_down_manual, `Cell type` = "DOWN_genes")
#gene.set_v1 <- rbind(gene.set_v1, gene.set_v2)
#rm(gene.set_v2)
colnames(gene.set) <- c("Metagene", "Cell type")
geneSet = split(1:dim(gene.set)[1], gene.set$`Cell type`) %>%
  lapply(function(x) {
    gene.set[unlist(x), 1]
  })
head(geneSet, 4)

# ssgsea
dat <- as.matrix(logTPM)
#ssgsea <- GSVA::gsva(dat, geneSet, method='ssgsea', kcdf="Gaussian", abs.ranking=T) #geneSet is preset
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = geneSet,
                       normalize = TRUE)
ssgsea <- gsva(gsvaPar, verbose = FALSE)
# Construct survival matrix 
Surv_matrix <- as.matrix(colData(data))
Surv_matrix <- as.data.frame(Surv_matrix)
c1 <- as.numeric(which(colnames(Surv_matrix) == "barcode"))
c2 <- as.numeric(which(colnames(Surv_matrix) == "patient"))
c3 <- as.numeric(which(colnames(Surv_matrix) == "vital_status"))
c4 <- as.numeric(which(colnames(Surv_matrix) == "days_to_death"))
Surv_matrix = Surv_matrix[, c(c1, c2, c3, c4)] # Select to columns needed 

class(Surv_matrix[, 1]) <- "character"
class(Surv_matrix[, 2]) <- "character"

Surv_matrix$vital_status[which(Surv_matrix$vital_status == "Alive")] <- "1"
Surv_matrix$vital_status[which(Surv_matrix$vital_status == "Dead")] <- "2"
#Surv_matrix[-which(Surv_matrix$vital_status == "Not Reported"), ]

colData(data)$vital_status

class(Surv_matrix$vital_status) <- "numeric"
class(Surv_matrix$days_to_death) <- "numeric" 
Surv_matrix$months_to_death <- Surv_matrix$days_to_death / 30
dd1 <- max(Surv_matrix$days_to_death, na.rm = T) + 100
Surv_matrix$days_to_death[is.na(Surv_matrix$days_to_death)] <- dd1
dd2 <- max(Surv_matrix$months_to_death, na.rm = T) + 1
Surv_matrix$months_to_death[is.na(Surv_matrix$months_to_death)] <- dd2
dd3 <- dd2 - .5

class(Surv_matrix$days_to_death) <- "numeric"
class(Surv_matrix$months_to_death) <- "numeric"

ssgsea <- t(ssgsea)
ssgsea <- as.data.frame(ssgsea)
ssgsea$barcode <- rownames(ssgsea)
Surv_matrix <- cbind(Surv_matrix, ssgsea)

# by arbitary grouping based on scores 
#for (j in c("DOWN_genes")) {
#  Count50 <- quantile(Surv_matrix[[j]], 0.75)
#  Surv_matrix[[paste0(j, "_status", sep = "")]] <- NA
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] > Count50)] <- "High"
#  Surv_matrix[[paste0(j, "_status", sep = "")]][which(Surv_matrix[[j]] < Count50)] <- "Low"
#}

#fit <- survfit(Surv(days_to_death, vital_status) ~ Iron_score_status, data = Surv_matrix)
#fit <- survfit(Surv(months_to_death, vital_status) ~ DOWN_genes_status, data = Surv_matrix)

res.cut <- surv_cutpoint(Surv_matrix, time = "months_to_death", event = "vital_status",
                         variables = c("DOWN_genes") # 找这3个变量的最佳切点
)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(months_to_death, vital_status) ~DOWN_genes, data = res.cat)

plot <- ggsurvplot(fit, 
                   #           surv.median.line = "hv",        
                   pval = T,                         
                   #conf.int = T,                    
                   risk.table = T,                    
                   risk.table.col = "strata",         
                   xlab = "Time (months)",             
                   ylab = "Survival Probability",    
                   xlim = c(0, dd3), 
                   legend.title = "cCD2",          
                   legend.labs = c("High", "Low"), 
                   ggtheme = theme(legend.position = "right", legend.box = "vertical", 
                                   panel.background = element_blank(),
                                   panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid"),
                                   legend.margin=margin(t= 0, unit='cm'), 
                                   panel.grid = element_blank(),
                                   title = element_blank(), 
                                   axis.ticks.x = element_blank(),
                                   legend.spacing = unit(0,"in"), 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),#x轴
                                   #        axis.text.y  = element_blank(), 
                                   #        axis.ticks = element_blank(),
                                   legend.text = element_text(size =12,color="black"), 
                                   legend.title = element_text(size =12,color="black"), 
                                   axis.title.y = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold"), 
                                   axis.title.x = element_text(vjust=1,  
                                                               size=0,
                                                               face = "bold")
                   ),
                   break.x.by = 12,               
                   palette = c("#5c696e", "darkgrey")) # "#00214F" for Th FOSB+ 

plot # 
#ggsave(plot, filename = paste(i, '_top_genes_surv.pdf'), width = 10, height = 7)

##### Save PDF
print("Printing PDF...")
pdf(file = paste0("cDC2_TARGET_OS", '_TARGET_OS_Surv_plot.pdf'),   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(plot, newpage = FALSE)
dev.off()



