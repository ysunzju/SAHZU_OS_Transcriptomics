library(here) # project-root-relative paths; run scripts from repository root
##################### Figure 2 and extended data fig 2 & 3 panels ##################### 
##################### Annotated chr and gene CNV heatmaps (Fig 2b and 2c) ##################### 
# gene level CNV heatmap ----- 
rm(list = ls()) 
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps"))
library(Seurat)
library(stringr)
scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Malig_OB_for_NMF_sample_w_too_few_cells_removed.rds"))

scrna$cell <- colnames(scrna)
master_meta_matrix <- data.frame(scrna@meta.data[c(47, 4)])
master_meta_matrix$patient <- master_meta_matrix$sample
master_meta_matrix$patient <- gsub("_BF_", "", master_meta_matrix$patient)
master_meta_matrix$patient <- gsub("_AF_", "", master_meta_matrix$patient)

#OS_genes <- readxl::read_xlsx(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/CNV_CNV_heatmaps/Curated_OS_genes.xlsx"))
#Oncogenes <- OS_genes$Oncogenes
#Suppressors <- OS_genes$`Tumour suppressors`[!is.na(OS_genes$`Tumour suppressors`)]

Oncogenes <- c(
  "AKT1", "ALK", "ALOX12B", "AURKB", "BRD4", "CARD11", "CALR",
  "CCND1", "CCND2", "CCND3", "CCNE1", "CD36", "CDK4", "CKDN1A",
  "COPS3", "CSMD3", "DNMT1", "FGF3", "FGF4", "FGF19", "FGFR1",
  "FLCN", "FOXM1", "GLI1", "HSP90AB1", "IL7R", "IGF1", "IGF1R",
  "INSR", "JUN", "KDR", "KEAP1", "KIT", "MAP2K4", "MCL1", "MDM2",
  "MET", "MYC", "NCOR1", "NOTCH3", "PDGFRA", "PIM1", "PRKDC",
  "PTPRD", "RAC1", "RAD21", "RICTOR", "RUNX2", "STAT6", "TERT",
  "TFDP1", "TMEM127", "VEGFA"
)

Suppressors <- c(
  "ATRX", "BRCA2", "CDKN2A", "CDKN2B", "CIC", "ERBB4", "EXT1",
  "FAS", "FAT1", "NF1", "NSD1", "PTEN", "RB1", "SPRED1", "TP53"
)

All_genes <- union(Oncogenes, Suppressors)

Patients <- unique(master_meta_matrix$patient)

## Loop for all patients (significant genes) ----- 
for (i in Patients) {
  ### Calc CNV event type ----- 
  gene_dat <- read.delim(paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), i, "/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.pred_cnv_genes.dat"))
  gene_dat_oi <- gene_dat[which(gene_dat$gene %in% All_genes), ]
  if (length(grep("ImmuneCells.ImmuneCells", gene_dat_oi$cell_group_name)) > 0) {
    gene_dat_oi <- gene_dat_oi[-grep("ImmuneCells.ImmuneCells", gene_dat_oi$cell_group_name), ]
  }
  cell_allocation <- read.delim(paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), i, "/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings"))
  cell_allocation <- cell_allocation[-grep("ImmuneCells.ImmuneCells", cell_allocation$cell_group_name), ] 
  
  patient_meta <- master_meta_matrix[which(master_meta_matrix$patient == i), ]
  patient_meta$sample <- droplevels(patient_meta$sample)
  patient_matrix <- merge(patient_meta, cell_allocation, by = "cell")
  patient_matrix_exten <- merge(patient_matrix, gene_dat_oi, by = "cell_group_name")
  patient_cell_count <- data.frame(table(patient_matrix$cell_group_name, patient_matrix$sample))
  
  patient_matrix_for_plot <- data.frame(matrix(ncol = length(All_genes), nrow = length(unique(patient_meta$sample))))
  colnames(patient_matrix_for_plot) <- All_genes
  rownames(patient_matrix_for_plot) <- unique(patient_meta$sample)
  
  # Nil CNV event = 3 
  for (j in All_genes[!All_genes %in% gene_dat_oi$gene]) {
    patient_matrix_for_plot[, j] <- 3
  }
  
  genes_for_score <- All_genes[All_genes %in% gene_dat_oi$gene]
  
  ## Pre 
  HL <- str_sub(string = i,start = 1,end = 1)
  Pt_num <- str_sub(string = i,start = 2,end = 3)
  Pre_label_name <- paste0(HL,"_BF_", Pt_num)
  patient_cell_count_pre <- patient_cell_count[which(patient_cell_count$Var2 == Pre_label_name), ]
  for (g in genes_for_score) {
    cell_line <- gene_dat_oi$cell_group_name[gene_dat_oi$gene == g] 
    pos_percent <- sum(patient_cell_count_pre$Freq[which(patient_cell_count_pre$Var1 %in% cell_line)]) / 
      sum(patient_cell_count_pre$Freq)
    pos_num <- sum(patient_cell_count_pre$Freq[which(patient_cell_count_pre$Var1 %in% cell_line)])
    max <- max(gene_dat_oi$state[gene_dat_oi$gene == g])
    max_cell_line <- gene_dat_oi$cell_group_name[which(gene_dat_oi$gene == g & gene_dat_oi$state == max)]
    max_percent <- sum(patient_cell_count_pre$Freq[which(patient_cell_count_pre$Var1 %in% max_cell_line)]) / 
      sum(patient_cell_count_pre$Freq[which(patient_cell_count_pre$Var1 %in% cell_line)])
    min <- min(gene_dat_oi$state[gene_dat_oi$gene == g])
    min_cell_line <- gene_dat_oi$cell_group_name[which(gene_dat_oi$gene == g & gene_dat_oi$state == min)]
    min_percent <- sum(patient_cell_count_pre$Freq[which(patient_cell_count_pre$Var1 %in% min_cell_line)]) / 
      sum(patient_cell_count_pre$Freq[which(patient_cell_count_pre$Var1 %in% cell_line)])
    states <- unique(gene_dat_oi$state[gene_dat_oi$gene == g])
    dom_count <- max(patient_cell_count_pre$Freq[patient_cell_count_pre$Var1 %in% cell_line])
    dom_cell_line <- as.character(patient_cell_count_pre$Var1[which(patient_cell_count_pre$Freq == dom_count)])
    dom_state <- gene_dat_oi$state[which(gene_dat_oi$cell_group_name == dom_cell_line & gene_dat_oi$gene == g)]
    
    if (all(pos_num > 10 & pos_percent > 0.05)) { 
      if (g %in% Oncogenes) { 
        if (all(max == states & length(states) == 1)) {
          patient_matrix_for_plot[Pre_label_name, g] <- max
        } else if (all(max_percent > 0.1 & length(states) != 1)) {
          patient_matrix_for_plot[Pre_label_name, g] <- as.numeric(paste0(max, ".", dom_state))
        } else if (max_percent <= 0.1) { 
          patient_matrix_for_plot[Pre_label_name, g] <- gene_dat_oi$state[which(gene_dat_oi$cell_group_name == dom_cell_line & 
                                                                                  gene_dat_oi$gene == g)]
        }
      } else if (g %in% Suppressors) { 
        if (all(min == states & length(states) == 1)) {
          patient_matrix_for_plot[Pre_label_name, g] <- min
        } else if (all(min_percent > 0.1 & length(states) != 1)) { 
          patient_matrix_for_plot[Pre_label_name, g] <- as.numeric(paste0(min, ".", dom_state))
        } else if (min_percent <= 0.1) { 
          patient_matrix_for_plot[Pre_label_name, g] <- dom_state
        }
      }
    } else { 
      patient_matrix_for_plot[, g] <- 3
    }
  }
  
  ## Post
  HL <- str_sub(string = i,start = 1,end = 1)
  Pt_num <- str_sub(string = i,start = 2,end = 3)
  Post_label_name <- paste0(HL,"_AF_", Pt_num)
  if (Post_label_name %in% rownames(patient_matrix_for_plot)) { # Only proceed if Post-treatment sample exists
    patient_cell_count_post <- patient_cell_count[which(patient_cell_count$Var2 == Post_label_name), ]
    for (g in genes_for_score) {
      cell_line <- gene_dat_oi$cell_group_name[gene_dat_oi$gene == g] 
      pos_percent <- sum(patient_cell_count_post$Freq[which(patient_cell_count_post$Var1 %in% cell_line)]) / 
        sum(patient_cell_count_post$Freq)
      pos_num <- sum(patient_cell_count_post$Freq[which(patient_cell_count_post$Var1 %in% cell_line)])
      max <- max(gene_dat_oi$state[gene_dat_oi$gene == g])
      max_cell_line <- gene_dat_oi$cell_group_name[which(gene_dat_oi$gene == g & gene_dat_oi$state == max)]
      max_percent <- sum(patient_cell_count_post$Freq[which(patient_cell_count_post$Var1 %in% max_cell_line)]) / 
        sum(patient_cell_count_post$Freq[which(patient_cell_count_post$Var1 %in% cell_line)])
      min <- min(gene_dat_oi$state[gene_dat_oi$gene == g])
      min_cell_line <- gene_dat_oi$cell_group_name[which(gene_dat_oi$gene == g & gene_dat_oi$state == min)]
      min_percent <- sum(patient_cell_count_post$Freq[which(patient_cell_count_post$Var1 %in% min_cell_line)]) / 
        sum(patient_cell_count_post$Freq[which(patient_cell_count_post$Var1 %in% cell_line)])
      states <- unique(gene_dat_oi$state[gene_dat_oi$gene == g])
      dom_count <- max(patient_cell_count_post$Freq[patient_cell_count_post$Var1 %in% cell_line])
      dom_cell_line <- as.character(patient_cell_count_post$Var1[which(patient_cell_count_post$Freq == dom_count)])
      dom_state <- gene_dat_oi$state[which(gene_dat_oi$cell_group_name == dom_cell_line & gene_dat_oi$gene == g)]
      
      if (all(pos_num > 10 & pos_percent > 0.01)) { 
        if (g %in% Oncogenes) { 
          if (all(max == states & length(states) == 1)) { 
            patient_matrix_for_plot[Post_label_name, g] <- max
          } else if (all(max_percent > 0.1 & length(states) != 1)) { 
            patient_matrix_for_plot[Post_label_name, g] <- as.numeric(paste0(max, ".", dom_state))
          } else if (max_percent <= 0.1) { 
            patient_matrix_for_plot[Post_label_name, g] <- gene_dat_oi$state[which(gene_dat_oi$cell_group_name == dom_cell_line & 
                                                                                     gene_dat_oi$gene == g)]
          }
        } else if (g %in% Suppressors) {
          if (all(min == states & length(states) == 1)) { 
            patient_matrix_for_plot[Post_label_name, g] <- min
          } else if (all(min_percent > 0.1 & length(states) != 1)) { 
            patient_matrix_for_plot[Post_label_name, g] <- as.numeric(paste0(min, ".", dom_state))
          } else if (min_percent <= 0.1) { 
            patient_matrix_for_plot[Post_label_name, g] <- dom_state
          }
        }
      } else { 
        patient_matrix_for_plot[, g] <- 3
      }
    }
  }
  
  if (i == "HP1") {
    All_patient_matrix_for_plot <- patient_matrix_for_plot
  } else if (i != "HP1") {
    All_patient_matrix_for_plot <- rbind(All_patient_matrix_for_plot, patient_matrix_for_plot)
  }
  
  ### Calc CNV percent ----- 
  patient_matrix_for_plot_percent <- data.frame(matrix(ncol = length(All_genes), nrow = length(unique(patient_meta$sample))))
  colnames(patient_matrix_for_plot_percent) <- All_genes
  rownames(patient_matrix_for_plot_percent) <- unique(patient_meta$sample)
  
  # Find most extreme event values 
  genes_for_score <- All_genes[All_genes %in% gene_dat_oi$gene] # Find genes with CNV 
  
  # If nil gene event then = 0 
  for (j in All_genes[!All_genes %in% gene_dat_oi$gene]) {
    patient_matrix_for_plot_percent[, j] <- 0
  }
  
  ## Pre 
  HL <- str_sub(string = i,start = 1,end = 1)
  Pt_num <- str_sub(string = i,start = 2,end = 3)
  Pre_label_name <- paste0(HL,"_BF_", Pt_num)
  patient_cell_count_pre <- patient_cell_count[which(patient_cell_count$Var2 == Pre_label_name), ]
  for (g in genes_for_score) {
    cell_line <- gene_dat_oi$cell_group_name[gene_dat_oi$gene == g] 
    pos_percent <- sum(patient_cell_count_pre$Freq[which(patient_cell_count_pre$Var1 %in% cell_line)]) / 
      sum(patient_cell_count_pre$Freq)
    patient_matrix_for_plot_percent[Pre_label_name, g] <- pos_percent
  }
  
  ## Post 
  HL <- str_sub(string = i,start = 1,end = 1)
  Pt_num <- str_sub(string = i,start = 2,end = 3)
  Post_label_name <- paste0(HL,"_AF_", Pt_num)
  if (Post_label_name %in% rownames(patient_matrix_for_plot_percent)) { # Only proceed if Post-treatment sample exists
    patient_cell_count_post <- patient_cell_count[which(patient_cell_count$Var2 == Post_label_name), ]
    for (g in genes_for_score) {
      cell_line <- gene_dat_oi$cell_group_name[gene_dat_oi$gene == g] 
      pos_percent <- sum(patient_cell_count_post$Freq[which(patient_cell_count_post$Var1 %in% cell_line)]) / 
        sum(patient_cell_count_post$Freq)
      patient_matrix_for_plot_percent[Post_label_name, g] <- pos_percent
    }
  }
  
  if (i == "HP1") {
    All_patient_matrix_for_plot_percent <- patient_matrix_for_plot_percent
  } else if (i != "HP1") {
    All_patient_matrix_for_plot_percent <- rbind(All_patient_matrix_for_plot_percent, patient_matrix_for_plot_percent)
  }
  
}

## Plot attempt ---- 
library(reshape2)

### Assemble CNV event matrix ----- 
All_patient_matrix_for_plot <- t(All_patient_matrix_for_plot)
ht <- reshape2::melt(All_patient_matrix_for_plot, value.name="State",na.rm = F)
colnames(ht)[1:2] <- c("sample","gene")
head(ht, 12)

# Change to factor
ht$gene <- factor(ht$gene,
                  levels = unique(ht$gene),
                  ordered = T)
ht$sample <- factor(ht$sample,
                    levels = rev(unique(ht$sample)),
                    ordered = T)
library(ggplot2)
ht$State_R <- ht$State
ht$State_R <- as.character(ht$State_R)
ht$State_R <- gsub("1", "I", ht$State_R)
ht$State_R <- gsub("2", "II", ht$State_R)
ht$State_R <- gsub("3", "III", ht$State_R)
ht$State_R <- gsub("4", "IV", ht$State_R)
ht$State_R <- gsub("5", "V", ht$State_R)
ht$State_R <- gsub("6", "VI", ht$State_R)
ht$State_R <- gsub("\\.", "-", ht$State_R)
ht$State_R <- gsub("I-I", "I", ht$State_R)
ht$State_R <- gsub("II-II", "II", ht$State_R)
ht$State_R <- gsub("III-III", "III", ht$State_R)
ht$State_R <- gsub("IV-IV", "IV", ht$State_R)
ht$State_R <- gsub("V-V", "V", ht$State_R)
ht$State_R <- gsub("VI-VI", "VI", ht$State_R)

### Assemble CNV percent matrix ----- 
All_patient_matrix_for_plot_percent <- t(All_patient_matrix_for_plot_percent)
ht1 <- reshape2::melt(All_patient_matrix_for_plot_percent, value.name="Percent",na.rm = F)
colnames(ht1)[1:2] <- c("sample","gene")
head(ht1, 12)

# Change to factor
ht1$gene <- factor(ht1$gene,
                   levels = unique(ht1$gene),
                   ordered = T)
ht1$sample <- factor(ht1$sample,
                     levels = rev(unique(ht1$sample)),
                     ordered = T)

### Combine and plot ----- 
ht$meta <- paste0(ht$sample, "_", ht$gene)
ht1$meta <- paste0(ht1$sample, "_", ht1$gene)
ht_comb <- merge(ht, ht1, by = "meta")
colnames(ht_comb)
ht_comb <- ht_comb[, -(6:7)]
colnames(ht_comb)[2:3] <- c("sample", "gene")
ht_comb$State <- floor(ht_comb$State)
ht_comb$patient <- gsub("_AF_", "", ht_comb$gene) %>% gsub("_BF_", "", .)

mat_v1 <- as.data.frame(table(ht_comb$sample[ht_comb$State_R != "III"]))
as.character(rev(mat_v1$Var1[order(mat_v1$Freq)]))
mat_sub <- mat_v1[mat_v1$Var1 %in% Oncogenes, ]
Onc <- as.character(rev(mat_sub$Var1[order(mat_sub$Freq)])) 
mat_sub <- mat_v1[mat_v1$Var1 %in% Suppressors, ]
Sup <- as.character(rev(mat_sub$Var1[order(mat_sub$Freq)])) 
Gene_ord <- rev(append(Onc, Sup))
ht_comb$sample <- factor(ht_comb$sample, levels = Gene_ord)

p <- ggplot(ht_comb, aes(gene, sample)) +
  geom_tile(aes(fill = State, height = Percent#scale(ht_comb$Percent, center = F)
  ), colour = "white", stat = "identity") + 
  #facet_grid(facets = ~feature.groups,  switch = "x", scales = "free_x", space = "free_x")
  #coord_equal() + 
  #geom_point(aes(fill = State, size = scale(ht_comb$Percent, center = F)), colour = "white", shape = 16) 
  geom_text(aes(label = ifelse(!State == 3, State_R, NA)), size = 3) + 
  facet_grid(~patient, scales = "free_x", space = "free_x") + 
  scale_fill_gradient2(low = "#50859f",
                       mid = "white",
                       high = "#d66692",
                       midpoint = 3) + 
  scale_size_continuous(guide = FALSE) +
  labs(x=NULL,y=NULL) + 
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  theme(axis.text.x.top = element_text(angle = 45,
                                       face = "italic",
                                       size = 8,
                                       hjust = 0,
                                       vjust = 0.01),
        axis.text.y = element_text(size = 8), 
        panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid"), 
        text = element_text(size = 8), 
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.position = "right", 
        panel.spacing.x = unit(0.2, "cm"),
        panel.spacing.y = unit(0.2, "cm"), 
        panel.background = element_blank()) + 
  NoLegend()
p

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_Updated_Tumour_gene_info_plot_260205.pdf"), 
    width = 10, 
    height = 15) 
print(p, newpage = FALSE)
dev.off()

## Per chr heatmap -----
scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Malig_OB_for_NMF_sample_w_too_few_cells_removed.rds"))

scrna$cell <- colnames(scrna)
master_meta_matrix <- data.frame(scrna@meta.data[c(47, 4)])
master_meta_matrix$patient <- master_meta_matrix$sample
master_meta_matrix$patient <- gsub("_BF_", "", master_meta_matrix$patient)
master_meta_matrix$patient <- gsub("_AF_", "", master_meta_matrix$patient)

hg38 = read.delim(here("data/Phylotree_annotation_script_28_05_25/cytoBand_hg38.txt"), header = F)
All_genes <- rep(1:22, 2) %>% as.character()
All_genes[1:22] <- paste0(All_genes[1:22], "p")
All_genes[23:44] <- paste0(All_genes[23:44], "q")

for (i in Patients) {
  ### Calc CNV event type ----- 
  gene_dat <- read.delim(paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), i, "/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.pred_cnv_regions.dat"))
  #hg19 = read.table("~/Downloads/cytoBand_hg19.txt")##change the object name to cytoband depending on which genome the annotation was done
  #gene_dat = read.table(here("data/InferCNV_rerun_for_phylo_22_06_25/HP1/HMM_CNV_predictions.HMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat"),header = T) ##from inferCNV resulting files
  if (length(grep("ImmuneCells.ImmuneCells", gene_dat$cell_group_name)) > 0) {
    gene_dat <- gene_dat[-grep("ImmuneCells.ImmuneCells", gene_dat$cell_group_name), ]
  }
  
  #cytoband = hg19
  cytoband = hg38
  cytoband <- data.frame(V1=gsub("chr", "", cytoband[,1]), V2=cytoband[,2], V3=cytoband[,3], V4=substring(cytoband$V4, 1, 1), stringsAsFactors=F)
  start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
  end <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
  cytoband <- data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
  cytoband <- cytoband [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
  cytoband$V4[grep("q", cytoband$V4)] <- "q"
  cytoband$V4[grep("p", cytoband$V4)] <- "p"
  rownames(cytoband) <- NULL
  names(cytoband) = c("chr", "start", "end", "arm")
  cytoband$chr_arm = paste0(cytoband$chr, cytoband$arm)
  cytoband$chromosome = paste0("chr", cytoband$chr)
  cytoband$chr = cytoband$chromosome
  cytoband = cytoband[,1:5]
  
  library(dplyr)
  gene_dat$chr = as.character(gene_dat$chr)
  cytoband$chr = as.character(cytoband$chr)
  gene_dat$mid = (as.numeric(gene_dat$start) + as.numeric(gene_dat$end))/2
  new = inner_join(gene_dat,cytoband, by=c("chr"))%>%
    mutate(arm = ifelse(mid >= as.numeric(start.y) & mid <= as.numeric(end.y), chr_arm, NA))%>%
    group_by(chr)%>%
    filter(!is.na(arm)|n()==1)
  final_cnv = new[,c(1:6, 10)]
  names(final_cnv) = c("cell_group_name", "cnv_name", "state", "chr", "start", "end", "arm")
  final_cnv$event = ifelse(final_cnv$state >=4, "gain", "loss")
  final_cnv$large_event = paste(final_cnv$arm, final_cnv$event, sep = " ")
  final_cnv = final_cnv[!duplicated(final_cnv[,c('cell_group_name', 'large_event')]),]
  cell_allocation <- read.delim(paste0(here("data/InferCNV_rerun_for_phylo_22_06_25/"), i, "/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings"))
  
  if (length(grep("ImmuneCells.ImmuneCells", cell_allocation$cell_group_name)) > 0) {
    cell_allocation <- cell_allocation[-grep("ImmuneCells.ImmuneCells", cell_allocation$cell_group_name), ] 
  }
  
  patient_meta <- master_meta_matrix[which(master_meta_matrix$patient == i), ]
  patient_meta$sample <- droplevels(patient_meta$sample)
  patient_matrix <- merge(patient_meta, cell_allocation, by = "cell")
  patient_matrix_exten <- merge(patient_matrix, final_cnv, by = "cell_group_name")
  patient_cell_count <- data.frame(table(patient_matrix$cell_group_name, patient_matrix$sample))
  sum(patient_cell_count$Freq)
  patient_matrix_for_plot <- data.frame(matrix(ncol = length(All_genes), nrow = length(unique(patient_meta$patient))))
  colnames(patient_matrix_for_plot) <- All_genes
  rownames(patient_matrix_for_plot) <- unique(patient_meta$patient)
  
  events_for_mat <- data.frame(matrix(ncol = 2, nrow = length(unique(patient_matrix_exten$large_event))))
  colnames(events_for_mat) <- c("Event", "YN")
  events_for_mat$Event <- unique(patient_matrix_exten$large_event)
  for (k in unique(patient_matrix_exten$large_event)) {
    clone_pres <- unique(patient_matrix_exten$cell_group_name[which(patient_matrix_exten$large_event == k)])
    cell_count_pres <- patient_cell_count$Freq[which(patient_cell_count$Var1 %in% clone_pres)]
    prop_pres <- sum(cell_count_pres)/sum(patient_cell_count$Freq)
    if (prop_pres > 0.2) {
      events_for_mat$YN[which(events_for_mat$Event == k)] <- T
    } else {
      events_for_mat$YN[which(events_for_mat$Event == k)] <- F
    }
  }
  
  sig_events <- events_for_mat$Event[events_for_mat$YN == T]
  patient_matrix_exten <- patient_matrix_exten[patient_matrix_exten$large_event %in% sig_events, ]
  
  for (j in unique(patient_matrix_exten$large_event)) {
    arm <- unique(patient_matrix_exten$arm[which(patient_matrix_exten$large_event == j)])
    event <- unique(patient_matrix_exten$event[which(patient_matrix_exten$arm == arm)])
    if (length(unique(patient_matrix_exten$event[which(patient_matrix_exten$arm == arm)])) > 1) {
      event <- "dual"
    }
    patient_matrix_for_plot[, which(colnames(patient_matrix_for_plot) == arm)] <- event
    patient_matrix_for_plot[is.na(patient_matrix_for_plot)] <- "Nil"
  }
  
  if (i == "HP1") {
    All_patient_matrix_for_plot <- patient_matrix_for_plot
  } else if (i != "HP1") {
    All_patient_matrix_for_plot <- rbind(All_patient_matrix_for_plot, patient_matrix_for_plot)
  }
  
}

# Plot ----- 
library(reshape2)
All_patient_matrix_for_plot <- t(All_patient_matrix_for_plot)
ht <- reshape2::melt(All_patient_matrix_for_plot, value.name="State",na.rm = F)
colnames(ht)[1:2] <- c("sample","gene")
head(ht, 12)

ht$arm <- substr(ht$sample, 2, 3)
for (i in 1:length(ht$arm)) {
  if (nchar(ht$arm[i]) > 1) {
    ht$arm[i] <- substr(ht$arm[i], 2, 2)
  }
}

ht$chr <- substr(ht$sample, 1, 2)
for (i in 1:length(ht$chr)) {
  if (grepl("p", substr(ht$sample, 1, 2))[i] == TRUE) {
    ht$chr[i] <- substr(ht$chr[i], 1, 1)
  }
  if (grepl("q", substr(ht$sample, 1, 2))[i] == TRUE) {
    ht$chr[i] <- substr(ht$chr[i], 1, 1)
  }
}

library(ggplot2)
library(cols4all)

mycol1 <- c4a('ag_sunset',10)
ht$chr <- factor(ht$chr, levels = rev(c(1:22)))
p <- ggplot(ht, aes(arm, chr)) +
  geom_tile(aes(fill = State, height = 0.3, width = 0.8
  ), colour = "white", stat = "identity") + 
  scale_fill_discrete(type = c("black", "#F09BA0", "#9BBBE1", "white")) +
  facet_grid(~gene, scales = "fixed", space = "fixed") + 
  scale_size_continuous(guide = FALSE) +
  labs(x=NULL,y=NULL) + 
  coord_equal(ratio = 0.8) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  theme(axis.text.x.top = element_text(
  ), 
  axis.ticks.x = element_blank(), 
  strip.background = element_blank(), 
  axis.ticks.y = element_blank(), 
  axis.text.y = element_text(size = 8), 
  panel.border = element_rect(fill = NA, color="white", linewidth = 1, linetype="solid"), 
  text = element_text(size = 8), 
  legend.title = element_text(size = 5),
  legend.text = element_text(size = 5),
  legend.position = "right", 
  panel.spacing.x = unit(0.08, "cm"),
  panel.spacing.y = unit(0.08, "cm"), 
  panel.background = element_blank()) 
p

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_per_Chr_CNV_heatmap.pdf"), 
    width = 6, 
    height = 3) 
print(p, newpage = FALSE)
dev.off()
                                                
##################### MP GO enrichment analysis (Fig 3b) ##################### 
MP_OS <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/NMF/Updated_MP_list.rds"))

### Construct the file for refererence 
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(R.utils)
library(tidyverse)
library(patchwork)

org.Hs.eg.db <- org.Hs.eg.db
gene <- MP_OS[[1]]
ego_ALL <- enrichGO(gene = gene, 
                    OrgDb=org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = "ALL", 
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.25,
                    readable = TRUE)
g1 <- as.data.frame(ego_ALL)
g1$logpvalue <- -log10(g1$pvalue)
colnames(g1)
g1 <- g1[order(g1$logpvalue, decreasing = T), ]

g2 <- g1
g2$Description <- capitalize(g2$Description)
g2$Cluster <- "MP1"

for (i in 2:9) {
  gene <- MP_OS[[i]]
  GO_temp <- enrichGO(gene = gene,
                      OrgDb=org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      minGSSize = 1,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.25,
                      readable = TRUE)
  GO_temp_df <- as.data.frame(GO_temp)
  GO_temp_df$logpvalue <- -log10(GO_temp_df$pvalue)
  GO_temp_df$Cluster <- paste("MP", i, sep = "") 
  g2 <- rbind(g2, GO_temp_df)
}

g2$Description <- capitalize(g2$Description)
which(g2$Description == c("Nuclear division")) # Manual selection of the pathways for display 

g3 <- g2[c(1, 5, 9, 
           356, 367, 369, 
           736, 739, 744, 
           840, 846, 851, 
           1043, 1100, 1104, 
           1129, 1130, 1149, 
           1343, 1351, 1455, 
           1953, 2006, 2044, 
           2616, 2623, 2629
), ]

go_data <- data.frame(GOItem = g3$Description,  
                      Pvalue = g3$p.adjust, 
                      CellType = g3$Cluster)
go_data$logP <- -log10(go_data$Pvalue)
go_data$GOItem <- gsub("-", " ", go_data$GOItem)

go_data$GOItem <- factor(go_data$GOItem, levels = rev(c(
  "Nuclear division", "Mitotic sister chromatid segregation", "Spindle",
  "DNA replication", "Recombinational repair", "DNA duplex unwinding", 
  "Ossification", "External encapsulating structure organization", "Bone mineralization", 
  "Extracellular matrix organization", "Collagen fibril organization", "Connective tissue development", 
  "Regulation of RNA splicing", "Chromatin protein adaptor activity",  "MiRNA binding", 
  "Oxidative phosphorylation", "Aerobic respiration", "Reactive oxygen species metabolic process", 
  "Integrated stress response signaling", "Cellular response to chemical stress", "Cellular response to environmental stimulus", 
  "MHC protein complex assembly", "Regulation of immune effector process", "Regulation of leukocyte cell cell adhesion", 
  "Response to type I interferon", "Regulation of innate immune response", "Pattern recognition receptor signaling pathway"
  )))

Color = colorRampPalette(c("#C8B5B3", "#594C57"))(9)
Color27 <- rep(Color, each = 3)
p2 <- ggplot(go_data, aes(y = GOItem)) +  geom_bar(aes(x = logP, fill = GOItem),            
                                                   stat = "identity",            
                                                   width = 0.5,         
                                                   color = "transparent",            
                                                   alpha = 0.7) +  
  geom_text(aes(x = 0.3, label = GOItem), hjust = 0, size = 4.5) +  
  labs(title = "GO enrichment item", x = "-log10(Pvalue)", y = NULL) +  
  scale_fill_manual(values = Color27) + theme_bw(base_size = 18) +  
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),    
        axis.text.x = element_text(color = "black"), axis.text.y = element_blank(), axis.ticks.y = element_blank(),    
        panel.grid = element_blank(), legend.position = "none")
p2

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_", "GO_enrichment_plot.pdf"), 
    width = 4, 
    height = 6) 
print(p2, newpage = FALSE)
dev.off()

##################### MP AddModulescore (Fig 3c-e) ##################### 
cancer.OB <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/Malig_OB_for_NMF_sample_w_too_few_cells_removed.rds"))
names(cancer.OB@meta.data)
cancer.OB <- AddModuleScore(cancer.OB, MP_OS, name =names(MP_OS))
#length(cancer.OB@meta.data) # 11105
names(cancer.OB@meta.data)[47:55] <- names(MP_OS)
#VlnPlot(cancer.OB, features = names(MP_OS), group.by = "group_anno", pt.size = 0)

# Violon plots of signatures across groups 
library(SCP)
cancer.OB[["RNA"]] <- as(object = cancer.OB[["RNA"]], Class = "Assay")
my_comparisons <- list(c("H_Pre", "H_Post"), c("H_Pre", "L_Pre"), c("L_Pre", "L_Post"), c("H_Post", "L_Post"))
pplist <- list()
for (i in names(MP_OS)) {
  pplist[[i]] <- FeatureStatPlot(srt = cancer.OB, group.by = "group_anno",  
                                 stat.by = i, 
                                 add_box = TRUE,  
                                 palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")), 
                                 bg_palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")),
                                 comparisons = my_comparisons,  
                                 add_trend = T, 
                                 sig_label = "p.format", 
                                 bg.by = 'group_anno') + theme(legend.position = "none") + 
    theme(legend.position = "none", 
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.title.y = element_blank())
}

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

pdf("Cancer_cell_NMF_scores.pdf", height = 8, width = 6)
print(pps)
dev.off()

# Featrue plot of MPs 
p <- FeatureDimPlot(cancer.OB, c("MP_3"#, "MP_9"
), ncol = 1, 
keep_scale = "all", 
add_density = F, #palcolor = c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124"),
theme_use = "theme_blank", 
palette = "material-grey", 
#raster = T, 
#pt.size = 5,
#raster.dpi = c(2048, 2048), 
reduction = "umap", 
color_blend_mode = "multiply"
)

pdf("Cycling_FeaturePlot.pdf", height = 3, width = 3)
print(p)
dev.off()

p <- FeatureDimPlot(cancer.OB, c("MP_9"#, "MP_9"
), ncol = 1, 
keep_scale = "all", 
add_density = F, #palcolor = c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124"),
theme_use = "theme_blank", 
palette = "material-grey", 
#raster = T, 
#pt.size = 5,
#raster.dpi = c(2048, 2048), 
reduction = "umap", 
color_blend_mode = "multiply"
)

pdf("S_Phase_FeaturePlot.pdf", height = 3, width = 3)
print(p)
dev.off()

p <- FeatureDimPlot(cancer.OB, c("MP_1"#, "MP_9"
), ncol = 1, 
keep_scale = "all", 
add_density = F, #palcolor = c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124"),
theme_use = "theme_blank", 
palette = "material-grey", 
#raster = T, 
#pt.size = 5,
#raster.dpi = c(2048, 2048), 
color_blend_mode = "multiply"
)

pdf("Osteo_FeaturePlot.pdf", height = 3, width = 3)
print(p)
dev.off()

p <- FeatureDimPlot(cancer.OB, c("MP_2"#, "MP_9"
), ncol = 1, 
keep_scale = "all", 
add_density = F, #palcolor = c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124"),
theme_use = "theme_blank", 
palette = "material-grey", 
#raster = T, 
#pt.size = 5,
#raster.dpi = c(2048, 2048), 
color_blend_mode = "multiply"
) 

pdf("Fibro_FeaturePlot.pdf", height = 3, width = 3)
print(p)
dev.off()

# MP correlation 
## Extract Cycling, S phase, OS/MIN, ECM
cellper <- data.frame(cancer.OB@meta.data)[c("MP_3", "MP_9", "MP_1", "MP_2")]

df <- cellper
res <- rcorr(as.matrix(df))
r_mat <- res$r
p_mat <- res$P

r_df <- reshape2::melt(r_mat, na.rm = FALSE)
p_df <- reshape2::melt(p_mat, na.rm = FALSE)

p_df$signif <- cut(p_df$value,
                   breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", ""))

plot_df <- merge(r_df, p_df, by = c("Var1", "Var2"))
plot_df$signif[plot_df$Var1 == plot_df$Var2] <- ""
plot_df$label <- paste0(sprintf("%.2f", plot_df$value.x)#, "\n", plot_df$signif
)
plot_df$label <- ifelse(plot_df$value.x > 0.3 | plot_df$value.x < -0.3, plot_df$label, "")

plot_df <- plot_df[as.numeric(plot_df$Var2) <= as.numeric(plot_df$Var1), ]

diagonal_labels <- subset(plot_df, Var1 == Var2)
diagonal_labels$label <- as.character(diagonal_labels$Var1)
diagonal_labels$y_pos <- as.numeric(diagonal_labels$Var1) - 1 

p <- ggplot() +
  geom_tile(data = plot_df, aes(x = Var2, y = Var1, fill = value.x), color = NA) +
  geom_text(data = plot_df, aes(x = Var2, y = Var1, label = label),
            colour = ifelse(plot_df$value.y < 0.05, "red", "black"), 
            #family = "Times New Roman", 
            size = 4) +
  #geom_text(data = diagonal_labels, aes(x = Var2, y = y_pos, label = label, hjust = 0.1, 
  #                                      vjust = 0.5
  #),
  #family = "Times New Roman",
  #size = 4) + 
  scale_fill_gradient2(
    low = "#50859f", high = "#d66692", mid = "white",
    midpoint = 0, limit = c(-1, 1), name = "Correlation",
    labels = scales::number_format(accuracy = 0.1)
  ) +
  coord_equal(ratio = 1) +
  # coord_fixed() +  
  #labs(title = "", x = "TNK", y = "Myeloid") +
  #theme_minimal(base_size = 14) +
  #expand_limits(y = max(as.numeric(plot_df_E1$Var1)) + 1)+
  theme(
    axis.text.y = element_text(),
    axis.ticks.y = element_blank(),
    # axis.text.x = element_blank(),
    panel.grid = element_blank(),
    #text = element_text(size = 4),
    plot.title = element_text(hjust = 0.5),
    #legend.position = "none"
  )
p

pdf(file = "OB_cells_correlation_260208.pdf", 
    width = 4, 
    height = 4) 
print(p, newpage = FALSE)
dev.off()
