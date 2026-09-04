library(here) # project-root-relative paths; run scripts from repository root
##################### InferCNV ##################### 
##### Identification of tumour cells #####
# SETUP -----------------------------------------------------------------------
# Libraries -------------------------------------------------------------------
rm(list = ls())
# disabled, run from repo root: setwd(here("data/infercnv"))
library(BiocManager)
library(plyr)
library(dplyr) 
library(Matrix)
library(Seurat)
library(ggplot2)
library(Rcpp)
library(RcppZiggurat)
library(Rfast)
library("SingleCellExperiment")
library(SummarizedExperiment)
library(BiocGenerics)
library(edgeR)
library(rjags)
library(infercnv)
library(igraph)
#library(reticulate)
#library(leiden) 
library(parallel)
detectCores()
numcores <- 30

# Load in Seurat objects -------------------------------------------------------
#setwd(here("data/infercnv"))
cancer.OS <- readRDS(here("data/infercnv/Malig_OB_for_NMF_sample_w_too_few_cells_removed.rds"))
DefaultAssay(cancer.OS) <- "RNA"
#Idents(combo.reference) <- "cluster"
#cancer.OS <- subset(combo.reference, ident = "Osteoblasts")
#Ref <- subset(combo.reference, ident = "Osteoclasts")
n_cells <- 2000
Ref1 <- readRDS(here("data/infercnv/Myeloid_reident_updated_doublets_removed_27_05_25.rds"))
Ref2 <- readRDS(here("data/infercnv/TNK_doublets_removed_reintegrated_v1_25_05_25.rds"))
Ref1 <- subset(Ref1, idents = c("Macro_FOLR2+", "Macro_TREM2+", "MoMac_ISGs", "Macro_SLC48A1+", 
                                "Macro_KLF2+", "cDC2_steady-state", "MoMac_late_activated", "Macro_ISGs", 
                                "Macro_SELENOP+", "Macro_SLC40A1hi", "cDC2_activated", "Macro_APOC1+", 
                                "moDC", "Macro_metallothionein-enriched", "MoMac_early_activated", 
                                "Mono_CD16+", "mDC", "Macro_glycolytic/inflammatory", "cDC1", "Macro_LYVE1", 
                                "Mono_classical", "Neu"))
Ref2 <- subset(Ref2, idents = c("HelperT_C03_stressed/late_activated", "Treg", "HelperT_C02_late_activated", "HelperT_C01_early_activated", 
                                "MAIT", "NaiveT"))
Ref <- merge(Ref1, y = Ref2)
Idents(Ref) <- "Sham_ident"
Ref <- subset(Ref, downsample = n_cells)
Ref$CellType <- "ImmuneCells"
Ref[["RNA"]] <- as(object = Ref[["RNA"]], Class = "Assay")
HLA_genes <- rownames(Ref)[(grep("^HLA-", rownames(Ref), invert=F))] # Obtain names of HLA genes 
Ref <- subset(Ref, features = setdiff(rownames(Ref), HLA_genes))
rm(combo.reference); rm(Ref1, Ref2); gc()
# saveRDS(Ref, "ref.rds")

# inferCNV by sample ----- 
for (i in unique(cancer.OS$sample)[c(unique(cancer.OS$sample)) %in% c("H_BF_P1", "H_BF_P2", "H_BF_P3", "H_BF_P4", "H_BF_P5", 
                                                                        "H_BF_P6", "H_BF_P7", "H_BF_P8", "H_AF_P1", "H_AF_P2", 
                                                                        "H_AF_P3", "H_AF_P4", "H_AF_P5", "H_AF_P6", "H_AF_P7", 
                                                                        "H_AF_P8", "L_BF_P1", "L_BF_P2", "L_BF_P3", "L_BF_P4", 
                                                                        "L_BF_P5", "L_BF_P6", "L_BF_P7", 
                                                                        "L_AF_P1", "L_AF_P2", 
                                                                        "L_AF_P3", "L_AF_P4", "L_AF_P5", "L_AF_P6", "L_AF_P7"
)]
) { 
  # generate input matrix and metadata files -------
  # load seurat object:
  scObject <- subset(cancer.OS, subset = sample == i)
  scObject$CellType <- "Tumor Candidate"
  scObject <- merge(x=scObject, y=Ref)
  scObject$anno <- as.character(scObject$CellType)
  write.table(scObject@meta.data[,"anno",drop=FALSE], paste0(i,"_annotations_file.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t" )
  ref_group_names <- setdiff(names(table(as.character(scObject$anno))),"Tumor Candidate")
  
  # Create infercnv_obj
  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix= GetAssayData(scObject, layer="counts"),
    annotations_file= paste0(i,"_annotations_file.txt"),
    delim="\t",
    #gene_order_file= here("data/InferCNV_31_01_25/infercnv_gene_order_nodup.txt"),
    gene_order_file= here("data/infercnv/hg38_gencode_v27.txt"), 
    #gene_order_file= here("data/InferCNV_31_01_25/infercnv_gene_order.txt"),
    #gene_order_file= here("data/InferCNV_31_01_25/gencode_v21_gen_pos.complete.txt"), 
    ref_group_names= ref_group_names,
    max_cells_per_group = NULL,
    min_max_counts_per_cell = c(100, +Inf),
    chr_exclude = c("chrX", "chrY", "chrM")
  )
  
  new_gene_order = data.frame()
  for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                     "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                     "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
    new_gene_order = rbind(new_gene_order, 
                           infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]]
                                                         == chr_name) , , drop=FALSE])
  }
  
  names(new_gene_order) <- c("chr", "start", "stop") 
  infercnv_obj@gene_order = new_gene_order 
  infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]
  
  infercnv_obj = infercnv::run( 
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10X Genomics
    out_dir= paste0(here("data/infercnv/"), i,"_inferCNV_rerun_28_05_25_immune_cells_w_T_HLA_removal"), 
    cluster_by_groups=FALSE, # If observations are defined according to groups (ie. patients), each group of cells will be clustered separately.
    denoise=TRUE,
    HMM=TRUE, # when set to True, runs HMM to predict CNV level (default: FALSE)
    analysis_mode="subclusters", # options(samples|subclusters|cells), Grouping level for image filtering or HMM predictions. default: samples (fastest, but subclusters is ideal)
    tumor_subcluster_partition_method="random_trees",
    scale_data = TRUE,
    leiden_resolution = 0.01,  
    k_obs_groups = 5,     
    output_format="pdf",  
    num_threads = numcores, 
    write_expr_matrix = T
  ) 
  
  gc()
}

##### Rerun for phylogeny (Now w/ malig cells of sufficient numbers only) #####
# Run as a seperate script 
# inferCNV script 28/06/25
rm(list = ls())
library(infercnv)
library(Seurat)

pat <- commandArgs(TRUE)[1]
numcores = 30

# Step1 : PrepareFiles
scObject1 <- readRDS( paste0("/data/s01020/InferCNV_phylo_rerun_28_06_25/tOS.",pat,".scObject.rds"))
Idents(scObject1) <- "cluster"
scObject1 <- subset(scObject1, idents = c("Osteoblasts"))
scObject1$CellType <- "Osteoblasts"
Ref <- readRDS(here("data/infercnv/ref.rds"))
Ref$CellType <- 'ImmuneCells'
scObject <- merge(x=scObject1, y=Ref)
scObject$anno <- as.character(scObject$CellType)
scObject$anno[scObject$anno %in% c("Osteoblasts")] <- "Tumor Candidate"
write.table(scObject@meta.data[,"anno",drop=FALSE], paste0(pat,"_annotations_file.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t" )
#table(scObject$CellType)
#table(scObject$anno)
ref_group_names <- setdiff(names(table(as.character(scObject$anno))),"Tumor Candidate")
ref_group_names

# Step2 : CreateInfercnvObject
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix= GetAssayData(scObject, layer="counts"),
  annotations_file= paste0(pat,"_annotations_file.txt"),
  delim="\t",
  #gene_order_file= here("data/InferCNV_31_01_25/infercnv_gene_order_nodup.txt"),
  gene_order_file= here("data/infercnv/hg38_gencode_v27.txt"), 
  #gene_order_file= here("data/InferCNV_31_01_25/infercnv_gene_order.txt"),
  #gene_order_file= here("data/InferCNV_31_01_25/gencode_v21_gen_pos.complete.txt"), 
  ref_group_names= ref_group_names,
  max_cells_per_group = NULL,
  min_max_counts_per_cell = c(100, +Inf),
  chr_exclude = c("chrX", "chrY", "chrM")
)

new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
  new_gene_order = rbind(new_gene_order, 
                         infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]]
                                                       == chr_name) , , drop=FALSE])
}

names(new_gene_order) <- c("chr", "start", "stop") 
infercnv_obj@gene_order = new_gene_order 
infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]

# Step3 : infercnv
infercnv_obj = infercnv::run( 
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10X Genomics
    out_dir= paste0(here("data/infercnv/"), pat), 
    cluster_by_groups=FALSE, # If observations are defined according to groups (ie. patients), each group of cells will be clustered separately.
    denoise=TRUE,
    HMM=TRUE, # when set to True, runs HMM to predict CNV level (default: FALSE)
    analysis_mode="subclusters", # options(samples|subclusters|cells), Grouping level for image filtering or HMM predictions. default: samples (fastest, but subclusters is ideal)
    tumor_subcluster_partition_method="random_trees",
    scale_data = TRUE,
    leiden_resolution = 0.01,  
    k_obs_groups = 5,     
    output_format="pdf",  
    num_threads = numcores, 
    write_expr_matrix = T
  ) 



