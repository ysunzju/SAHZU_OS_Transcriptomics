library(here) # project-root-relative paths; run scripts from repository root
##### Contains records of manual identification of malignant cells #####
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer"))
library(infercnv)
library(Seurat)
library(SeuratObject)
#rm(list=ls())
#seurat_obj <- readRDS("seurat_obj.RDS")
#cancer.epi <- readRDS("All_OB_cells_reintegrated.rds")
DefaultAssay(cancer.epi) <- "RNA"
sample_names <- names(table(cancer.epi$sample))
cluster_group_all <- data.frame()

# Manual selection of malignant cells (Related to Figure S2A-D) ----- 
## Sample H_BF_P1 ----- 
sample <- sample_names[which(sample_names == "H_BF_P1")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}
cluster1
cluster2 <- "Malignant" # As discussed, all cells as malignant for now 
cluster3
cluster4
cluster5

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_AF_P1 ----- 
sample <- sample_names[which(sample_names == "H_AF_P1")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}
cluster1 <- "Malignant" # As discussed, all cells as malignant for now 
cluster2 <- "Malignant" # As discussed, all cells as malignant for now 
cluster3 <- "Malignant" # As discussed, all cells as malignant for now 
cluster4
cluster5

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_BF_P2 ----- 
sample <- sample_names[which(sample_names == "H_BF_P2")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}
cluster1 <- "Malignant" # As discussed 
cluster2 
cluster3 
cluster4
cluster5 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_AF_P2 ----- 
sample <- sample_names[which(sample_names == "H_AF_P2")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}
cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_BF_P3 ----- 
sample <- sample_names[which(sample_names == "H_BF_P3")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}
cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_AF_P3 ----- 
sample <- sample_names[which(sample_names == "H_AF_P3")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed 
# Checked and aligns with pre

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_BF_P4 ----- 
sample <- sample_names[which(sample_names == "H_BF_P4")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed 
# Checked and aligns with pre

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_AF_P4 ----- 
sample <- sample_names[which(sample_names == "H_AF_P4")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed, too few malignant cells, treat as non-malignant 
# Checked and aligns with pre

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_BF_P5 ----- 
sample <- sample_names[which(sample_names == "H_BF_P5")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed, too few malignant cells, treat as non-malignant 
# Checked and aligns with pre

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_AF_P5 ----- 
sample <- sample_names[which(sample_names == "H_AF_P5")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed, too few malignant cells, treat as non-malignant 
# Checked and aligns with pre

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_BF_P6 ----- 
sample <- sample_names[which(sample_names == "H_BF_P6")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed
# As is as discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_AF_P6 ----- 
sample <- sample_names[which(sample_names == "H_AF_P6")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_BF_P7 ----- 
sample <- sample_names[which(sample_names == "H_BF_P7")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_AF_P7 ----- 
sample <- sample_names[which(sample_names == "H_AF_P7")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_BF_P8 ----- 
sample <- sample_names[which(sample_names == "H_BF_P8")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample H_AF_P8 ----- 
sample <- sample_names[which(sample_names == "H_AF_P8")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_BF_P1 ----- 
sample <- sample_names[which(sample_names == "L_BF_P1")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_AF_P1 ----- 
sample <- sample_names[which(sample_names == "L_AF_P1")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_BF_P2 ----- 
sample <- sample_names[which(sample_names == "L_BF_P2")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_AF_P2 ----- 
sample <- sample_names[which(sample_names == "L_AF_P2")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_BF_P3 ----- 
sample <- sample_names[which(sample_names == "L_BF_P3")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_AF_P3 ----- 
sample <- sample_names[which(sample_names == "L_AF_P3")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_BF_P4 ----- 
sample <- sample_names[which(sample_names == "L_BF_P4")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_AF_P4 ----- 
sample <- sample_names[which(sample_names == "L_AF_P4")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_BF_P5 ----- 
sample <- sample_names[which(sample_names == "L_BF_P5")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 #<- "Non-Malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_AF_P5 ----- 
sample <- sample_names[which(sample_names == "L_AF_P5")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_BF_P6 ----- 
sample <- sample_names[which(sample_names == "L_BF_P6")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_AF_P6 ----- 
sample <- sample_names[which(sample_names == "L_AF_P6")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_BF_P7 ----- 
sample <- sample_names[which(sample_names == "L_BF_P7")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 #<- "Non-malignant"
cluster2 #<- "Non-malignant"
cluster3 #<- "Non-malignant" # As discussed 
cluster4 #<- "Non-malignant"
cluster5 <- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

## Sample L_AF_P7 ----- 
sample <- sample_names[which(sample_names == "L_AF_P7")] 
dir <- paste0(here("data/InferCNV_31_01_25/"), sample, "_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal")
dd <- sample
cnv_table1 <- read.table(paste(dir, "General_HCL_1_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table2 <- read.table(paste(dir, "General_HCL_2_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table3 <- read.table(paste(dir, "General_HCL_3_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table4 <- read.table(paste(dir, "General_HCL_4_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table5 <- read.table(paste(dir, "General_HCL_5_members.txt", sep="/"), header=T, sep = " ", check.names = FALSE)
cnv_table1_m <- as.matrix(cnv_table1)
cnv_table2_m <- as.matrix(cnv_table2)
cnv_table3_m <- as.matrix(cnv_table3)
cnv_table4_m <- as.matrix(cnv_table4)
cnv_table5_m <- as.matrix(cnv_table5)

cnv_score_table1 <- cnv_table1_m
cnv_score_mat1 <- cnv_table1_m

# Scoring
# CNV??5????ϵͳ
cnv_score_table1[cnv_score_mat1 == min(cnv_score_mat1)] <- 2 #complete loss. 2pts
cnv_score_table1[cnv_score_mat1 != min(cnv_score_mat1) & cnv_score_mat1 != max(cnv_score_mat1)] <- 0 #Neutral. 0pts
cnv_score_table1[cnv_score_mat1 == max(cnv_score_mat1)] <- 2 #addition of two copies. 2pts

cnv_score_table2 <- cnv_table2_m
cnv_score_mat2 <- cnv_table2_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table2[cnv_score_mat2 == min(cnv_score_mat2)] <- 2 #complete loss. 2pts
cnv_score_table2[cnv_score_mat2 != min(cnv_score_mat2) & cnv_score_mat2 != max(cnv_score_mat2)] <- 0 #Neutral. 0pts
cnv_score_table2[cnv_score_mat2 == max(cnv_score_mat2)] <- 2 #addition of two copies. 2pts

cnv_score_table3 <- cnv_table3_m
cnv_score_mat3 <- cnv_table3_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table3[cnv_score_mat3 == min(cnv_score_mat3)] <- 2 #complete loss. 2pts
cnv_score_table3[cnv_score_mat3 != min(cnv_score_mat3) & cnv_score_mat3 != max(cnv_score_mat3)] <- 0 #Neutral. 0pts
cnv_score_table3[cnv_score_mat3 == max(cnv_score_mat3)] <- 2 #addition of two copies. 2pts

cnv_score_table4 <- cnv_table4_m
cnv_score_mat4 <- cnv_table4_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table4[cnv_score_mat4 == min(cnv_score_mat4)] <- 2 #complete loss. 2pts
cnv_score_table4[cnv_score_mat4 != min(cnv_score_mat4) & cnv_score_mat4 != max(cnv_score_mat4)] <- 0 #Neutral. 0pts
cnv_score_table4[cnv_score_mat4 == max(cnv_score_mat4)] <- 2 #addition of two copies. 2pts

cnv_score_table5 <- cnv_table5_m
cnv_score_mat5 <- cnv_table5_m
# Scoring
# CNV??5????ϵͳ
cnv_score_table5[cnv_score_mat5 == min(cnv_score_mat5)] <- 2 #complete loss. 2pts
cnv_score_table5[cnv_score_mat5 != min(cnv_score_mat5) & cnv_score_mat5 != max(cnv_score_mat5)] <- 0 #Neutral. 0pts
cnv_score_table5[cnv_score_mat5 == max(cnv_score_mat5)] <- 2 #addition of two copies. 2pts

mean1 <- c()
for (i in 1:ncol(cnv_score_table1)) {
  mean1 <- rbind(mean1, mean(cnv_score_table1[,i]))
}
mean2 <- c()
for (i in 1:ncol(cnv_score_table2)) {
  mean2 <- rbind(mean2, mean(cnv_score_table2[,i]))
}
mean3 <- c()
for (i in 1:ncol(cnv_score_table3)) {
  mean3 <- rbind(mean3, mean(cnv_score_table3[,i]))
}
mean4 <- c()
for (i in 1:ncol(cnv_score_table4)) {
  mean4 <- rbind(mean4, mean(cnv_score_table4[,i]))
}
mean5 <- c()
for (i in 1:ncol(cnv_score_table5)) {
  mean5 <- rbind(mean5, mean(cnv_score_table5[,i]))
}

if(all(c(max(mean1), max(mean2), max(mean3), max(mean4), max(mean5))>1)){
  cutoff <- 0.5
}else{
  cutoff <- 0.2
}

count <- 0
count_list1 <- c()
index1 <- which(mean1 > cutoff)
if(length(index1)>1){
  for(i in 1:(length(index1)-1)){
    if(index1[i+1]-index1[i] == 1){
      if(i == length(index1)-1){
        count <- count+2
        count_list1 <- c(count_list1, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list1 <- c(count_list1, count)
      count <- 0
    }
  }
}

count <- 0
count_list2 <- c()
index2 <- which(mean2 > cutoff)
if(length(index2)>1){
  for(i in 1:(length(index2)-1)){
    if(index2[i+1]-index2[i] == 1){
      if(i == length(index2)-1){
        count <- count+2
        count_list2 <- c(count_list2, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list2 <- c(count_list2, count)
      count <- 0
    }
  }
}

count <- 0
count_list3 <- c()
index3 <- which(mean3 > cutoff)
if(length(index3)>1){
  for(i in 1:(length(index3)-1)){
    if(index3[i+1]-index3[i] == 1){
      if(i == length(index3)-1){
        count <- count+2
        count_list3 <- c(count_list3, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list3 <- c(count_list3, count)
      count <- 0
    }
  }
}

count <- 0
count_list4 <- c()
index4 <- which(mean4 > cutoff)
if(length(index4)>1){
  for(i in 1:(length(index4)-1)){
    if(index4[i+1]-index4[i] == 1){
      if(i == length(index4)-1){
        count <- count+2
        count_list4 <- c(count_list4, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list4 <- c(count_list4, count)
      count <- 0
    }
  }
}

count <- 0
count_list5 <- c()
index5 <- which(mean5 > cutoff)
if(length(index5)>1){
  for(i in 1:(length(index5)-1)){
    if(index5[i+1]-index5[i] == 1){
      if(i == length(index5)-1){
        count <- count+2
        count_list5 <- c(count_list5, count)
        break
        
      }else{
        count <- count+1
      }
    }else{
      count <- count+1
      count_list5 <- c(count_list5, count)
      count <- 0
    }
  }
}

gene_cutoff <- 10
if(any(count_list1 >= gene_cutoff)){
  cluster1 <- "Malignant"
}else{
  cluster1 <- "Non-malignant"
}

if(any(count_list2 >= gene_cutoff)){
  cluster2 <- "Malignant"
}else{
  cluster2 <- "Non-malignant"
}

if(any(count_list3 >= gene_cutoff)){
  cluster3 <- "Malignant"
}else{
  cluster3 <- "Non-malignant"
}

if(any(count_list4 >= gene_cutoff)){
  cluster4 <- "Malignant"
}else{
  cluster4 <- "Non-malignant"
}

if(any(count_list5 >= gene_cutoff)){
  cluster5 <- "Malignant"
}else{
  cluster5 <- "Non-malignant"
}

cluster1 <- "Non-malignant"
cluster2 <- "Non-malignant"
cluster3 <- "Non-malignant" # As discussed 
cluster4 <- "Non-malignant"
cluster5 #<- "Non-malignant" # As discussed
# As discussed 

cluster_group <- read.table(paste(dir, "infercnv.observation_groupings.txt", sep = '/'), header=T, check.names=FALSE)
#table(cluster_group$`Dendrogram Group`, cluster_group$`Dendrogram Color`)
colnames(cluster_group) <- c("cell_type", "b", "c", "d")
cluster_group[cluster_group$cell_type==1,1] <- cluster1
cluster_group[cluster_group$cell_type==2,1] <- cluster2
cluster_group[cluster_group$cell_type==3,1] <- cluster3
cluster_group[cluster_group$cell_type==4,1] <- cluster4
cluster_group[cluster_group$cell_type==5,1] <- cluster5
cluster_group1 <- cbind(cell=rownames(cluster_group), cluster_group)
cluster_group2 <- cluster_group1[,1:2]
#table(cluster_group2$cell_type)

write.table(cluster_group2, file=paste0(dd, "_Malignant_cell_type.txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cluster_group3 <- cbind(cluster_group2, dd)

cluster_group_all <- rbind(cluster_group_all, cluster_group3)

# Save master file ----- 
write.table(cluster_group_all, file="all_Malignant_cell_type.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
