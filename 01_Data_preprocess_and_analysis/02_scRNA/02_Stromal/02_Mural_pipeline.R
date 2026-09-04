library(here) # project-root-relative paths; run scripts from repository root
##################### Mural cell pipeline ##################### 
rm(list = ls())
# disabled, run from repo root: setwd(here("data/Re-clustering_fine_cell_types_21_05_25/MuralCells"))
library(Seurat)
library(tidyverse)

# Using existing Identification previously ----- 
scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/MuralCells/MuralCs_doublets_removed_reintegrated_v2_11_07_25.rds"))
scrna <- NormalizeData(scrna) %>% FindVariableFeatures() %>% ScaleData()

## RNA_snn_res.0.5 ---- 
DimPlot(scrna, group.by = "RNA_snn_res.0.5")
mar <- FindAllMarkers(scrna, group.by = "RNA_snn_res.0.5", only.pos = T) 
Idents(scrna) <- "RNA_snn_res.0.5"
scrna <- RenameIdents(scrna, "0" = "SMC_MYH11+")
scrna <- RenameIdents(scrna, "1" = "apPC") # MALAT, CD74
scrna <- RenameIdents(scrna, "2" = "matPC_ATF3+")
scrna <- RenameIdents(scrna, "3" = "myoPC_IGFBP2+") # SYNPO2, IGFBP2
scrna <- RenameIdents(scrna, "4" = "matPC_COL6A3+") # FOSB
scrna <- RenameIdents(scrna, "5" = "matPC_BASP+") #
scrna <- RenameIdents(scrna, "6" = "myoPC_CXCL12+") # 
scrna <- RenameIdents(scrna, "7" = "SMC_BTG2+") # BTG2, MYH11hi 
scrna <- RenameIdents(scrna, "8" = "infPC") 
#scrna <- RenameIdents(scrna, "9" = "MC_ISGs") 
scrna <- RenameIdents(scrna, "10" = "cyclMC") 
scrna <- RenameIdents(scrna, "11" = "myoPC_IGFBP2+") 
scrna <- RenameIdents(scrna, "12" = "cyclMC") 

scrna$temp <- Idents(scrna)
levels(scrna)
Idents(scrna) <- factor(Idents(scrna), levels = c("cyclMC", "matPC_BASP+", "matPC_COL6A3+", 
                                                  "matPC_ATF3+", "infPC", "apPC", "myoPC_IGFBP2+", 
                                                  "myoPC_CXCL12+", "SMC_MYH11+", "SMC_BTG2+"
))
scrna$Cell_type_fine_harmony <- Idents(scrna)
scrna <- RenameIdents(scrna, "matPC_BASP+" = "matPC", 
                      "matPC_COL6A3+" = "matPC", "matPC_ATF3+" = "matPC", 
                      "myoPC_IGFBP2+" = "MyoPC", "myoPC_CXCL12+" = "MyoPC", 
                      "SMC_MYH11+" = "SMC", "SMC_BTG2+" = "SMC")
levels(scrna)
#[1] "matPC"  "MyoPC"  "SMC"    "cyclMC" "infPC"  "apPC"  
Idents(scrna) <- factor(Idents(scrna), levels = c("cyclMC", "matPC", "infPC", "apPC", "MyoPC", "SMC"))
scrna$Cell_type_med_harmony <- Idents(scrna)
#saveRDS(scrna, "Mural_fine_and_med_updated_260105.rds") 

## Mural cell UMAP (Fig 4L) ----- 
levels(scrna)
colors <- c("#9d3b62", "#405993", "#2d3462","#3ca0cf", "#64a776", "#e0cfda",
            "#696a6c", "#927c9a", "#6c408e", "#d25774")  
mycol <- colors
names(mycol) <- levels(scrna$Cell_type_fine_harmony)

p1 <- DimPlot(scrna, reduction = "umap", group.by = c("Cell_type_fine_harmony"), 
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

pdf(paste0("Mural_cells", "_UMAP_v251231.pdf"), width = 7, height = 7, family = "Helvetica")
print(p1)
dev.off() 

## vSMC and Pericyte and other signatures plot (Figure 4F,G) ----- 
library(Nebulosa)
p <- plot_density(scrna, c("CD248", "COL4A2", "THY1", "BASP1", 
                           "COL6A3", "ATF3", "MKI67", "IFIT1", "CD74", 
                           "SYNPO2", "IGFBP2", "CXCL12", 
                           "ACTA2", "MYH11", "BTG2"
), reduction = "umap", combine = F, raster = T)

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

f <- plot_grid(p[[1]], 
               p[[2]],
               p[[3]], 
               p[[4]],
               p[[5]],
               p[[6]],
               p[[7]],
               p[[8]],
               p[[9]],
               p[[10]],
               p[[11]], 
               p[[12]], 
               p[[13]], 
               p[[14]], 
               p[[15]], 
               ncol = 5 #, nrow = 2
)
f

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "FeaturePlot_MCs_cluster_marker_sig.pdf"), 
    width = 15, 
    height = 9) 
print(f, newpage = FALSE)
dev.off()

# Plot Peri, VSMC sigs, and angiogenic genes 
library(Nebulosa)
p <- plot_density(scrna, c("Pericytes", "VSMC", "PDF", "PDGFA", "PDGFB"
), reduction = "umap", combine = F, raster = T
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

f <- plot_grid(p[[1]], 
               p[[2]],
               p[[3]], 
               p[[4]], 
               p[[5]], 
               ncol = 5 #, nrow = 2
)
f

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "FeaturePlot_MCs_pericytes_VSMC_sig.pdf"), 
    width = 15, # The width of the plot in inches
    height = 3) # The height of the plot in inches
print(f, newpage = FALSE)
dev.off()

## Mural and endothelial cell correlation (Figure 4M) ----- 
library(Seurat)
library(reshape2)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(tidyverse)

#scRNA <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))
Endo <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Endo/EC_w_Cap_updated_251231.rds"))
Mural <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/MuralCells/Mural_fine_and_med_updated_260105.rds"))
Mural$Cell_type_med_harmony <- droplevels(Mural$Cell_type_med_harmony)
Endo$Cell_type_fine_harmony <- droplevels(Endo$Cell_type_fine_harmony)
Endo$Cell_type_fine_harmony <- factor(Endo$Cell_type_fine_harmony, 
                                      levels = c("artEC", "Tip-like", "Transition_FLT1hi", "Transition_FLT1lo", 
                                                 "Stalk-like", "venEC", "LEC",  "EndoMT-I", "EndoMT-II", "Cycling_ECs"))

Cellratio_1 <- prop.table(table(Endo$Cell_type_fine_harmony, Endo$sample), margin = 2) %>% data.frame()
cellper_i <- reshape2::dcast(Cellratio_1, Var2~Var1, value.var="Freq")
rownames(cellper_i)<- cellper_i[,1]
cellper_i <- cellper_i[,-1]

Cellratio_2 <- prop.table(table(Mural$Cell_type_med_harmony, Mural$sample), margin = 2)  %>% data.frame()
cellper_j <- reshape2::dcast(Cellratio_2, Var2~Var1, value.var="Freq")
rownames(cellper_j)<- cellper_j[,1]
cellper_j <- cellper_j[,-1] 

cellper <- cbind(cellper_i, cellper_j)
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

plot_df <- plot_df[plot_df$Var1 %in% unique(Endo$Cell_type_fine_harmony), ]
plot_df <- plot_df[plot_df$Var2 %in% unique(Mural$Cell_type_med_harmony), ]

plot_df <- plot_df[plot_df$Var1 %in% c("artEC", "Tip-like", "Transition_FLT1hi", "Transition_FLT1lo", 
                                       "Stalk-like", "venEC"), ]
plot_df <- plot_df[plot_df$Var2 %in% c("matPC", "infPC", "apPC", "MyoPC", "SMC"), ]

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
  coord_equal(ratio = 0.5) +
  theme_minimal(base_size = 14) +
  expand_limits(y = max(as.numeric(plot_df$Var1)) + 1)+
  theme(
    axis.text.y = element_text(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank(),
    #text = element_text(size = 4),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "EC_and_Mural_subsetsOI_", "_co-occur_plot.pdf"), 
    width = 10, 
    height = 6) 
print(p, newpage = FALSE)
dev.off() 

## matPC cell percentages per group (Figure 4P) ----- 
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

matPC <- readRDS("Mural_fine_and_med_updated_260105.rds")
Idents(matPC) <- "group"
matPC <- RenameIdents(matPC, "A" = "H_Pre", 
                      "B" = "H_Post", 
                      "C" = "L_Pre", 
                      "D" = "L_Post")
matPC$group_anno <- Idents(matPC)
matPC$Cell_type_fine_harmony <- droplevels(matPC$Cell_type_fine_harmony)
Cellratio <- prop.table(table(matPC$Cell_type_fine_harmony, matPC$sample), margin = 2) 
Cellratio <- data.frame(Cellratio)

cellper <- reshape2::dcast(Cellratio, Var2~Var1, value.var="Freq")
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]

#write.csv()
#to_save <- t(cellper)
write.csv(to_save, file='matPC_as_mural_cells_prop_260806.csv', quote = F)

meta <- matPC@meta.data
colnames(meta)
meta <- meta[,c(32,4)] #group and sample columns 
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
#cellper$group <- NA
meta <- unique(meta)
cellper <- merge(cellper, meta, by = "sample")
colnames(cellper)
colnames(cellper)[5] <- "group" # Group anno to group 
cellper <- as.data.frame(cellper)

pplist =list()
seuratObj_groups = unique(levels(matPC$Cell_type_fine_harmony))

mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
for(group_ in seuratObj_groups){
  cellper_ = cellper %>% dplyr::select(one_of(c('sample', 'group', group_)))
  colnames(cellper_)=c('sample','group','percent')
  cellper_$percent =as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper = quantile(percent,0.75),
                                                      lower = quantile(percent,0.25),
                                                      mean = mean(percent),
                                                      median = median(percent),
                                                      lower_lim = range(percent)[1], 
                                                      upper_lim = range(percent)[2])#上下分位数
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
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
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
                 ncol = 3 #, nrow = 2
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_Mural_cell_fine_type_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 4.75, 
    height = 3) 
print(pps, newpage = FALSE)
dev.off()

## matPC vasogenic factor expressions (Fig 4Q) -----
matPC <- readRDS("Mural_fine_and_med_updated_260105.rds")
matPC <- subset(matPC, idents = "matPC")
matPC <- NormalizeData(matPC) %>% FindVariableFeatures() %>% ScaleData()

# BY CELL TYPE
matPC$Cell_type_fine_harmony <- droplevels(matPC$Cell_type_fine_harmony)
matPC[["RNA"]] <- as(object = matPC[["RNA"]], Class = "Assay")
p1 <- GroupHeatmap(matPC,
                   features = c("VEGFA", "VEGFB", "PGF"),
                   group.by = c("Cell_type_fine_harmony"), 
                   split.by = "group_anno", 
                   add_dot = FALSE, 
                   heatmap_palcolor = c("#50859f", "white", "#d66692")
) 

p1

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_matPC_VEGF_by_cell_type_x_group.pdf"), 
    width = 10, # The width of the plot in inches
    height = 3) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off() 
