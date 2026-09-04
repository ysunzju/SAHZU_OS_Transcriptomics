library(here) # project-root-relative paths; run scripts from repository root
##################### Myeloid cell pipeline #####################
rm(list = ls())

# Annotation ----- 
library(Seurat)
library(tidyverse)
# disabled, run from repo root: setwd(here("data/Per_cell_type_pipeline_w_new_idents/Myeloid"))
scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Myeloid/All_myeloid_reident_updated_27_05_25.rds"))
Idents(scrna) <- "RNA_snn_res.0.3"

scrna <- RenameIdents(scrna, "0" = "Macro_APOE+", 
                      "1" = "Macro_CCL4+", 
                      "2" = "Mono_FCN1+",
                      "3" = "cDC2_CD1C+", 
                      "4" = "Macro_ISGs+", 
                      "5" = "Macro_Ki67+",
                      "6" = "Stromal_doub", 
                      "7" = "Macro_MTs+",
                      "8" =  "Macro_APOE+", 
                      "9" = "Neu", 
                      "10" = "Promyelocytes", 
                      "11" = "cDC1_CLEC9A+", 
                      "12" = "Mature_DC")
DimPlot(scrna, label = T)
scrna <- subset(scrna, ident = "Stromal_doub", invert = T)
scrna$Cell_type_fine_harmony <- Idents(scrna)

DimPlot(scrna, label = T, group.by = "Cell_type_fine_harmony")
Idents(scrna) <- "Cell_type_fine_harmony"
#DimPlot(scrna, label = T)
markers <- FindAllMarkers(scrna,
                          logfc.threshold = 0.25, 
                          min.pct = 0.2,
                          min.diff.pct = 0.1,
                          only.pos = T)
## UMAP (Figure 6A) -----
celltype_colors <- mycol <- c(
  "#C65D57",  
  "#D4A373",  
  "#7D7461",  
  "#9B6A6C",  
  "#5B5F97",  
  "#DAA49A",  
  "#8D6A4F",  
  "#A4C2A5",  
  "#B07BAC",  
  "#7C9D96",  
  "#D8B4A0",  
  "#645D5C"   
)
scrna$Cell_type_fine_harmony <- Idents(scrna)
scrna$Cell_type_fine_harmony <- factor(scrna$Cell_type_fine_harmony, levels = 
                                         c("Mono_FCN1+", "Macro_CCL4+", "Macro_ISGs+", "Macro_MTs+", 
                                           "Macro_APOE+", 
                                           "Macro_Ki67+", "Promyelocytes", "Neu", "cDC1_CLEC9A+", "cDC2_CD1C+", "Mature_DC"
                                         ))
pdf("Myeloid_Dimplot_v2_01_05_25.pdf", height = 7, width = 7)
DimPlot(scrna, group.by = "Cell_type_fine_harmony", 
        label = T, raster = T, pt.size = 3, 
        cols = celltype_colors, 
        reduction = "umap", 
        raster.dpi = c(2048, 2048)) +
  theme_void() + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.title = element_blank(), 
    legend.position = "none"
  )
dev.off()

## Cell count 
dev.off()
reduc <- data.frame(Seurat::Embeddings(scrna, reduction = "umap"))
meta <- scrna@meta.data

pc12 <- cbind(reduc, meta)
pc12$idents <- Idents(scrna)

cell_num <- pc12 |>
  dplyr::group_by(idents,Cell_type_fine_harmony) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::arrange(n)

cols <- mycol[cell_num$Cell_type_fine_harmony]
pushViewport(viewport(x = unit(0.61, "npc"), y = unit(0.5, "npc"),
                      width = unit(0.2, "npc"),
                      height = unit(0.7, "npc"),
                      just = "left",
                      # yscale = extendrange(c(0,nrow(cell_num)),f = 0.05),
                      # xscale = extendrange(c(0,max(cell_num$n)),f = 0.05),
                      yscale = c(0,nrow(cell_num) + 0.5),
                      xscale = c(0,max(cell_num$n) + 0.1*max(cell_num$n))))
grid.xaxis()
# grid.yaxis(main = F)
# grid.rect()
grid.rect(x = rep(0,nrow(cell_num)),y = 1:nrow(cell_num),
          width = cell_num$n,height = unit(0.08,"npc"),
          just = "left",
          gp = gpar(fill = cols,col = NA),
          default.units = "native") 
grid.rect(gp = gpar(fill = "transparent"))
grid.text(label = "Number of cells",x = 0.5,y = unit(-2.5,"lines"))
popViewport()

pushViewport(viewport(x = unit(0.81, "npc"), y = unit(0.5, "npc"),
                      width = unit(0.2, "npc"),
                      height = unit(0.7, "npc"),
                      just = "left",
                      yscale = c(0,nrow(cell_num) + 0.5)))
# grid.rect()
grid.points(x = rep(0.1,nrow(cell_num)),y = 1:nrow(cell_num),pch = 19,
            gp = gpar(col = cols),
            size = unit(1.5, "char"))
cell_num$numbering <- c(1:nrow(cell_num))
grid.text(label = cell_num$numbering,x = 0.1,y = 1:nrow(cell_num),
          default.units = "native")
grid.text(label = cell_num$idents,x = 0.2,y = 1:nrow(cell_num),
          just = "left",
          default.units = "native")

popViewport()
# Manual saves 

## Mono, Macro, and DC signature feature plot (Figure 6B) -----
# 2 gene myeloid sig 
names(scrna@meta.data)
#scrna@meta.data[c(47:49)] <- NULL
Macro_genes <- c("C1QB", "C1QA")
Mono_genes <- c("VCAN", "FCN1")
DC_genes <- c("FCER1A", "CD1C", "BATF3", "CLEC9A", "LAMP3", "CCR7")

myeloid.genes <- list(Macro_genes, Mono_genes, DC_genes)
names(myeloid.genes) <- c("Macro", "Mono", "DC")
scrna <- AddModuleScore(scrna,
                        features = myeloid.genes,
                        name =names(myeloid.genes))
names(scrna@meta.data)[c(47:49)] <- c("Macro", "Mono", "DC")
FeaturePlot(scrna, c("Macro", "Mono", "DC"), ncol = 3)

pdf("Macro_sig.pdf", height = 5, width = 5)
SCpubr::do_FeaturePlot(scrna, "Macro", pt.size = 0.1) + 
  scale_color_gradientn(colors = c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124")) +
  theme(legend.position = "none")
dev.off()

pdf("Mono_sig.pdf", height = 5, width = 5)
SCpubr::do_FeaturePlot(scrna, "Mono", pt.size = 0.1) + 
  scale_color_gradientn(colors = c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124")) +
  theme(legend.position = "none")
dev.off()

pdf("DC_sig.pdf", height = 5, width = 5)
SCpubr::do_FeaturePlot(scrna, "DC", pt.size = 0.1) + 
  scale_color_gradientn(colors = c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124")) +
  theme(legend.position = "none")
dev.off()

## DotPlot of marker genes (Figure 6C) -----
DefaultAssay(scrna) <- "RNA"
levels(scrna$Cell_type_fine_harmony)
Mono_FCN1_genes <- c("FCN1", "VCAN", "CD300E") #
MAC_CCL4_genes <- c("CCL4", "KLF2", "KLF4") #
MAC_ISGs_genes <- c("IFIT1","CXCL10", "CXCL11") #
Macro_MTs_genes <- c("SPP1", "MT1X", "MT1E") # 
cDC1_CLEC9A_genes <- c("CLEC9A","XCR1","BATF3") # 
cDC2_CD1C_genes <- c("FCER1A","CD1C","CLEC10A") # 
mDC_genes <- c("LAMP3","CCR7","CCL19") # 
MAC_Ki67_genes <- c("MKI67","TOP2A","STMN1") # 
Promyelocyte_genes <- c("MPO","MYB","PCLAF") # 
MAC_APOE_genes <- c("APOC1","APOE","LIPA") 
Neu_genes <- c("S100A8","S100A9","CSF3R") # 

levels(scrna) <- levels(scrna$Cell_type_fine_harmony)
features <- list("Mono_FCN1" = Mono_FCN1_genes,                  
                 "MAC_CCL4" = MAC_CCL4_genes,                 
                 "MAC_ISGs" = MAC_ISGs_genes,                  
                 "Macro_MTs" = Macro_MTs_genes,                 
                 "MAC_APOE" = MAC_APOE_genes, 
                 "MAC_Ki67" = MAC_Ki67_genes, 
                 "Promyelocyte" = Promyelocyte_genes, 
                 "Neu" = Neu_genes, 
                 "cDC1_CLEC9A" = cDC1_CLEC9A_genes, 
                 "cDC2_CD1C" = cDC2_CD1C_genes, 
                 "mature DC" = mDC_genes
)

p <- DotPlot(object = scrna, features = features)

p1 <- ggplot(p$data, aes(x = features.plot, y = id)) +   
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +   
  facet_grid(facets = ~feature.groups,  switch = "x", scales = "free_x", space = "free_x") +    
  scale_radius(breaks = c(25, 50, 75, 100), range = c(0,6)) +   
  theme_classic() + 
  scale_color_gradient2(low = "#50859f", mid = "white", high = "#d66692") +  
  theme(axis.text.x = element_text(angle = 90,face = 1, size = 12, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, face = 1, color = "black"), 
        panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype="solid"),
        legend.text = element_text(size = 8, face = 1),         
        legend.title = element_text(size = 10, face = 1),        
        legend.position = 'top', strip.placement = "outside", strip.text.x = element_blank(),        
        axis.title = element_blank()) + guides(colour = guide_colourbar(title.vjust = 0.9, title.hjust = 0)) + 
  labs(size = "Percent Expressed", color = "Average Expression") 
p1

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_", "Myeloid_feature_bubbleplot.pdf"), 
    width = 9, # The width of the plot in inches
    height = 7) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()

### Myeloid signature heatmap (Figure 6D)----- 
#https://www.sciencedirect.com/science/article/pii/S0092867421000106
#https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00593-1
#https://www.cell.com/cancer-cell/fulltext/S1535-6108(24)00440-9#mmc4
library(msigdbr)
hs_df = msigdbr(species = "Homo sapiens") %>% as.data.frame()
hs_react = msigdbr(species = "Homo sapiens",
                   category = "C2",
                   subcategory = "CP:REACTOME") %>% as.data.frame() %>% 
  dplyr::select(gs_cat, gs_subcat, gs_name, gene_symbol)

hs_react = hs_react %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  dplyr::rename("term" = "gs_name", "gene" = "gene_symbol")
rm(hs_df)
genes <- list("M1" = c("IL23A", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD86", 
                       "IL1A", "IL1B", "IL6", "CCL5", "IRF5", "IRF1", 
                       "CD40", "IDO1", "KYNU", "CCR7"),
              "M2" = c("IL4R", "CCL4", "CCL13", "CCL20", "CCL17", "CCL18", "CCL22", "CCL24",
                       "LYVE1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "EGF", "CTSA", "CTSB",
                       "CTSC", "CTSD", "TGFB1", "TGFB2", "TGFB3", "MMP14", "MMP19", "MMP9",
                       "CLEC7A", "WNT7B", "FASLG", "TNFSF12", "TNFSF8", "CD276", "VTCN1",
                       "MSR1", "FN1", "IRF4"),
              "Angiogenic" = c("CCND2", "CCNE1", "CD44", "CXCR4", "E2F3", "EDN1", "EZH2", "FGF18", "FGFR1", "FYN",
                               "HEY1", "ITGAV", "JAG1", "JAG2", "MMP9", "NOTCH1", "PDGFA", "PTK2", "SPP1", "STC1",
                               "TNFAIP6", "TYMP", "VAV2", "VCAN", "VEGFA"),
              "Phagocytic" = c("MRC1", "CD163", "MERTK", "C1QB"), 
              "MHCII" <- c("HLA-DRA", "HLA-DRB5", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DQB1-AS1", "HLA-DQA2", "HLA-DQB2",
                           "HLA-DOB", "HLA-DMB", "HLA-DMA", "HLA-DOA", "HLA-DPA1", "HLA-DPB1", "CD74"), 
              "T_Cell_Attraction" <- c("CCL21", "CCL17", "CCL2","CXCL9", "CXCL10", "CXCL11", "CXCL12","CCL3", "CCL4", "CCL5", "CXCL16"), 
              "T_Cell_Repression" <- c("IDO1", "CD274", "PDCD1LG2", "TNFSF10","HLA-E", "HLA-G", "VTCN1", "IL10","TGFB1", "TGFB2", "PTGS2", "PTGES","LGALS9", "CCL22", "CD80", "CD86"), 
              "CD28_CO_STIMULATION" = hs_react$gene[which(hs_react$term == "REACTOME_CD28_CO_STIMULATION")], 
              "JAK_STAT_post_IL12" = hs_react$gene[which(hs_react$term == "REACTOME_GENE_AND_PROTEIN_EXPRESSION_BY_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_12_STIMULATION")])

names(scrna@meta.data)
scrna <- AddModuleScore(scrna,
                        features =genes,
                        name =names(genes))
names(scrna@meta.data)[c(47:55)] <- c("M1", "M2", "Angiogenic", "Phagocytic", "MHCII", "T_Attr", "T_Repre", "CD28_Costim", "IL12_JAK_STAT")

genelistnames <- c("M1", "M2", "Angiogenic", "Phagocytic", "MHCII", "T_Attr", "T_Repre", "CD28_Costim", "IL12_JAK_STAT")
heatmap_data <- scrna@meta.data %>%
  select(Cell_type_fine_harmony, all_of(genelistnames)) %>%
  group_by(Cell_type_fine_harmony) %>%
  summarise(across(all_of(genelistnames), mean, na.rm = TRUE)) %>%
  column_to_rownames("Cell_type_fine_harmony") %>%
  as.matrix() %>%
  t()  

heatmap_data_scaled <- t(scale(t(heatmap_data)))

p_value_matrix <- scrna@meta.data %>%
  group_by(Cell_type_fine_harmony) %>%
  summarise(across(all_of(genelistnames), ~t.test(.x)$p.value)) %>%
  column_to_rownames("Cell_type_fine_harmony") %>%
  as.matrix() %>%
  t()

library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

cellwidth = 0.7
cellheight = 0.7
cn = dim(as.matrix(heatmap_data_scaled))[2]
rn = dim(as.matrix(heatmap_data_scaled))[1]
w=cellwidth*cn
h=cellheight*rn

pdf(file = "Myeloid_sigs.pdf", width = 5, height = 15) 
Heatmap(as.matrix(heatmap_data_scaled),
        # Height of the blocks
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        #scale = F,
        #rect_gp = gpar(col = "white", lwd = 1.5),
        #border_g = gpar(col = ,lty = 1,lwd = 1.2),
        # Dent formatting 
        column_dend_height = unit(1.5, "cm"), 
        row_dend_width = unit(1.5, "cm"),
        #column_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
        #row_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
        # Gap numbers and style 
        #row_split = 6, column_split = c("1", "1", "2", "2"),
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        # Text formats 
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        show_heatmap_legend = T, 
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = 'right',
        show_column_names = T,
        show_row_names = T,
        border = T, 
        col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#50859f","white","#d66692")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if( heatmap_data_scaled[i, j] > 1) {
            grid.text("+++", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] > 0.8 &  heatmap_data_scaled[i, j] <= 1) {
            grid.text("++", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= 0.8 &  heatmap_data_scaled[i, j] > 0.2) {
            grid.text("+", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= 0.2 &  heatmap_data_scaled[i, j] > 0) {
            grid.text("+/-", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] == 0) {
            grid.text("-", x, y, gp = gpar(fontsize = 8))
          }
        }
)
dev.off()

## Cell percentages (Figure 6E) ----- 
scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Myeloid/Myeloid_final_as_at_24_07.rds"))
Idents(scrna) <- "group"
scrna <- RenameIdents(scrna, "A" = "H_Pre", 
                      "B" = "H_Post", 
                      "C" = "L_Pre", 
                      "D" = "L_Post")
scrna$group_anno <- Idents(scrna)
scrna$Cell_type_fine_harmony <- droplevels(scrna$Cell_type_fine_harmony)
Cellratio <- prop.table(table(scrna$Cell_type_fine_harmony, scrna$sample), margin = 2)  
Cellratio <- data.frame(Cellratio)

cellper <- reshape2::dcast(Cellratio, Var2~Var1, value.var="Freq") 
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]

#write.csv()
#to_save <- t(cellper)
write.csv(to_save, file='Myeloid_cells_prop_260806.csv', quote = F)

meta <- scrna@meta.data
colnames(meta)
meta <- meta[,c(56,4)] # group and sample
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
#cellper$group <- NA
meta <- unique(meta)
cellper <- merge(cellper, meta, by = "sample")
colnames(cellper)
colnames(cellper)[13] <- "group" # Group anno to group 
cellper <- as.data.frame(cellper)

pplist =list()
scrna_groups = unique(levels(scrna$Cell_type_fine_harmony))

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

mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
#fill_color <- zzm60colors
#names(fill_color) <- unique(levels(scrna$cluster))

for(group_ in scrna_groups){
  cellper_ = cellper %>% dplyr::select(one_of(c('sample', 'group', group_))) 
  colnames(cellper_)=c('sample','group','percent') 
  cellper_$percent =as.numeric(cellper_$percent) 
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper = quantile(percent,0.75),
                                                      lower = quantile(percent,0.25),
                                                      mean = mean(percent),
                                                      median = median(percent),
                                                      lower_lim = range(percent)[1], 
                                                      upper_lim = range(percent)[2]) 
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
                          paired = T, 
                          test = "t.test",
                          test.args = "two.sided", 
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  #my_comparisons1 <-list(c("L_Pre","L_Post"))
  pp1 = pp1 + geom_signif(comparisons = list(c("L_Pre","L_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))}, 
                          y_position = p_position, 
                          test = "t.test",
                          test.args = "two.sided", 
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("H_Pre","L_Pre")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position1, 
                          test = "t.test", 
                          test.args = "two.sided", 
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black") 
  #pp1 = pp1 + geom_signif(comparisons = list(c("H_Post","L_Post")), 
  #                        map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
  #                        y_position = p_position2, 
  #                        test = "t.test", 
  #                        test.args = "two.sided", 
  #                        textsize = 4, tip_length = 0,
  #                        parse = TRUE, 
  #                        color = "black")
  
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
                 pplist[[4]],
                 pplist[[5]],
                 pplist[[6]],
                 pplist[[7]],
                 pplist[[8]],
                 pplist[[9]],
                 pplist[[10]],
                 pplist[[11]], 
                 ncol = 6 #, nrow = 2
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_myeloid_cell_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 12, # The width of the plot in inches
    height = 6) # The height of the plot in inches
print(pps, newpage = FALSE)
dev.off()

## Monocle 3 trajectory (Fig 6F,G, Figure S6B,C) -----
### Macro ----- 
library(monocle3)
Idents(scrna) <- "Cell_type_fine_harmony"
unique(Idents(scrna))
scrna_subset <- subset(scrna, idents = c("Mono_FCN1+", "Macro_ISGs+", "Macro_CCL4+", "Macro_MTs+", 
                                         "Macro_APOE+"))
data <- GetAssayData(scrna_subset, assay = 'RNA', layer = 'counts')
cell_metadata <- scrna_subset@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(expression_dat = data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds) # Preprocessing data 
cds <- align_cds(cds, alignment_group = "sample") # Remove batch effects
cds <- reduce_dimension(cds) # Dimension reduction
cds <- cluster_cells(cds) # Cluster cells 
cds <- learn_graph(cds) # Predict trajectory 
#cds <- order_cells(cds) # Manual selection 

get_earliest_principal_node <- function(cds, time_bin="Mono_FCN1+"){
  cell_ids <- which(colData(cds)[, "Cell_type_fine_harmony"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

mycol <- celltype_colors
names(mycol) <- levels(scrna$Cell_type_fine_harmony)
mycol <- mycol[levels(scrna$Cell_type_fine_harmony)]

p1 <- plot_cells(cds,
                 color_cells_by = "Cell_type_fine_harmony",
                 label_cell_groups=F,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size = 1.5, 
                 group_label_size = 2,
                 trajectory_graph_color = "#EAEAEA", 
                 trajectory_graph_segment_size = 0.5, 
                 cell_size = 0.5) + 
  scale_color_manual(values = mycol[c("Mono_FCN1+", "Macro_ISGs+", "Macro_CCL4+", "Macro_MTs+", 
                                      "Macro_APOE+")]) + 
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.line.y = element_blank(), 
        panel.border = element_blank())
p1

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_cell_type_plots_for_myeloid.pdf"), 
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()

p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5, 
                 trajectory_graph_color = "#CCCCCC", 
                 trajectory_graph_segment_size = 0.5)
p2 <- p1 + theme(legend.position = "none",
                 axis.ticks.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.title.x = element_blank(), 
                 axis.line.x = element_blank(), 
                 axis.line.y = element_blank(), 
                 panel.border = element_blank()) +
  scale_color_gradientn(colours = c("#fef9ef",  "#f9dbb4",  "#f1b782",  "#e67c54",  "#d13f2d",  "#a2292d",  "#7c2124"))
p2

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_Myeloid.pdf"), 
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(p2, newpage = FALSE)
dev.off()

#### Graph test  ----- 
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))

library(ClusterGVis)
mat <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)
cluster.num = 3
exp <- mat
hclust_matrix <- exp; rm(exp); rm(mat)
km <- stats::kmeans(x = hclust_matrix, centers = cluster.num, nstart = 10)
od.res <- data.frame(od = match(names(km$cluster), rownames(hclust_matrix)), 
                     id = as.numeric(km$cluster), check.names = FALSE)
cl.info <- data.frame(table(od.res$id), check.names = FALSE)
m <- hclust_matrix[od.res$od, ]
wide.r <- m %>% data.frame(check.names = FALSE) %>% dplyr::mutate(gene = rownames(.), 
                                                                  cluster = od.res$id) %>% dplyr::arrange(cluster)
df <- reshape2::melt(wide.r, id.vars = c("cluster", "gene"), 
                     variable.name = "cell_type", value.name = "norm_value")
df$cluster_name <- paste("cluster ", df$cluster, sep = "")
cltn <- table(wide.r$cluster)
df <- purrr::map_df(unique(df$cluster_name), function(x) {
  tmp <- df %>% dplyr::filter(cluster_name == x)
  cn = as.numeric(unlist(strsplit(as.character(x), 
                                  split = "cluster "))[2])
  tmp %>% dplyr::mutate(cluster_name = paste(cluster_name, 
                                             " (", cltn[cn], ")", sep = ""))
})
df$cluster_name <- factor(df$cluster_name, levels = paste("cluster ", 
                                                          1:nrow(cl.info), " (", cl.info$Freq, ")", sep = ""))
wide <- wide.r
wide$cluster <- paste("C", wide$cluster, sep = "")
cluster.list <- split(wide$gene, wide$cluster)
ck <- list(wide.res = wide.r, long.res = df, cluster.list = cluster.list, 
           type = "kmeans", geneMode = "none", geneType = "none")

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "_ClusterVis_Macro_genes_Moran_0.25", "_v1.pdf"), 
    width = 3, # The width of the plot in inches
    height = 6) # The height of the plot in inches
visCluster(object = ck,
           plot.type = "heatmap",
           add.sampleanno = F,
           ht.col.list = list(col_range = c(-2, 0, 2), col_color =  c("#e9f4f6", "#8f9fab", "#1e2235"))#, 
           #markGenes = sample(rownames(mat),30,replace = F
)
dev.off()

#### GO enrichment ----
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(cols4all)
hs_df = msigdbr(species = "Homo sapiens") %>% as.data.frame()
hs_C5 = msigdbr(species = "Homo sapiens", category = "C5", subcategory = NULL) %>% as.data.frame() %>% dplyr::select(gs_cat, gs_subcat, gs_name, gene_symbol)
hs_C5 = hs_C5 %>% dplyr::select(gs_name, gene_symbol) %>% dplyr::rename("term" = "gs_name", "gene" = "gene_symbol")
rm(hs_df)
organism = 'hsa'    
OrgDb = 'org.Hs.eg.db'
org.Hs.eg.db <- org.Hs.eg.db

# C1 - Late
gene <- ck$cluster.list$C1
enricher_res_C1 <- enricher(gene, TERM2GENE = hs_C5, pvalueCutoff = 1,
                            pAdjustMethod = "BH")
gsea_results_C1 <- enricher_res_C1@result
gsea_results_C1$pval_log <- -log10(gsea_results_C1$pvalue)
#gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]

go_select <- gsea_results_C1 %>% 
  filter(pvalue < 0.05) %>% 
  arrange(pvalue) %>% 
  head(50)  %>% 
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

go_select <- go_select[which(go_select$ID %in% c("GOBP_DEFENSE_RESPONSE_TO_VIRUS", "GOBP_LEUKOCYTE_PROLIFERATION", "GOMF_CYTOKINE_ACTIVITY", "GOMF_CHEMOKINE_ACTIVITY", 
                                                 "GOBP_MONOCYTE_CHEMOTAXIS", "GOBP_HUMORAL_IMMUNE_RESPONSE")), ]
go_select$pval_log <- -log10(go_select$pvalue)
scale_factor <- max(go_select$pval_log) / max(go_select$RichFactor)

mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11),
                 plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5))
mytheme2 <- mytheme + theme(axis.text.y = element_blank()) 

p1 <- ggplot(data = go_select, aes(x = -log10(pvalue), y = Description, fill = Count)) +
  #scale_fill_continuous_c4a_seq('pu_bu') +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9, fill = "#8f9fab") +
  labs(x = 'RichFactor', y = '') +
  geom_text(aes(x = 0.03, 
                label = Description),
            hjust = 0)+ 
  theme_classic() + mytheme2
p1

# C2 - Early
gene <- ck$cluster.list$C2
enricher_res_C2 <- enricher(gene, TERM2GENE = hs_C5, pvalueCutoff = 1,
                            pAdjustMethod = "BH")
gsea_results_C2 <- enricher_res_C2@result
gsea_results_C2$pval_log <- -log10(gsea_results_C2$pvalue)

# C3 - Term
gene <- ck$cluster.list$C3
enricher_res_C3 <- enricher(gene, TERM2GENE = hs_C5, pvalueCutoff = 1,
                            pAdjustMethod = "BH",)
gsea_results_C3 <- enricher_res_C3@result
gsea_results_C3$pval_log <- -log10(gsea_results_C3$pvalue)

GO_terms <- unique(c(gsea_results_C2$Description, gsea_results_C1$Description, gsea_results_C3$Description))

mat_for_heatmap <- as.data.frame(matrix(ncol = 3, nrow = length(GO_terms)))
rownames(mat_for_heatmap) <- GO_terms
colnames(mat_for_heatmap) <- c("Early", "Late", "Terminal")
mat_for_heatmap$Early[which(rownames(mat_for_heatmap) %in% gsea_results_C2$Description)] <- gsea_results_C2[rownames(mat_for_heatmap)[which(rownames(mat_for_heatmap) %in% gsea_results_C2$Description)], "pval_log"]
mat_for_heatmap$Late[which(rownames(mat_for_heatmap) %in% gsea_results_C1$Description)] <- gsea_results_C1[rownames(mat_for_heatmap)[which(rownames(mat_for_heatmap) %in% gsea_results_C1$Description)], "pval_log"]
mat_for_heatmap$Terminal[which(rownames(mat_for_heatmap) %in% gsea_results_C3$Description)] <- gsea_results_C3[rownames(mat_for_heatmap)[which(rownames(mat_for_heatmap) %in% gsea_results_C3$Description)], "pval_log"]
mat_for_heatmap[is.na(mat_for_heatmap)] <- 0 
mat_for_heatmap <- as.matrix(mat_for_heatmap)

library(ComplexHeatmap)
Heatmap(mat_for_heatmap[, 1:3], show_row_names = F, cluster_rows = F, cluster_columns = F)
heatmap(mat_for_heatmap)

#GO_to_keep <- unique(c("GOBP_TAXIS", "GOBP_HUMORAL_IMMUNE_RESPONSE", "GOBP_NEGATIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY", "GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",
#                       "GOBP_MAINTENANCE_OF_LOCATION", "GOBP_DEFENSE_RESPONSE_TO_VIRUS", "GOBP_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE", 
#                       "GOBP_NEGATIVE_REGULATION_OF_T_CELL_PROLIFERATION", "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX", "GOBP_REGULATION_OF_T_CELL_ACTIVATION", 
#                       "GOBP_MYELOID_LEUKOCYTE_DIFFERENTIATION", "GOBP_REGULATION_OF_LEUKOCYTE_MEDIATED_IMMUNITY"
#))

#GO_to_keep <- unique(c("GOBP_CELL_CHEMOTAXIS", "GOBP_HUMORAL_IMMUNE_RESPONSE", "GOMF_CYTOKINE_ACTIVITY", 
#                       "GOBP_LEUKOCYTE_AGGREGATION", "GOBP_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY", 
#                       "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR", "GOBP_INTERLEUKIN_2_PRODUCTION", 
#                       "GOBP_ALPHA_BETA_T_CELL_ACTIVATION", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
#                       "GOBP_LEUKOCYTE_CELL_CELL_ADHESION", "GOBP_NEGATIVE_REGULATION_OF_LYMPHOCYTE_ACTIVATION", 
#                       "GOBP_NEGATIVE_REGULATION_OF_DEFENSE_RESPONSE",
#                       "GOBP_NEGATIVE_REGULATION_OF_T_CELL_PROLIFERATION", "GOBP_MYELOID_LEUKOCYTE_ACTIVATION", 
#                       "GOBP_REGULATION_OF_LIPID_TRANSPORT", "GOBP_REGULATION_OF_LIPID_TRANSPORT", 
#                       "GOBP_INTERLEUKIN_10_PRODUCTION", "GOBP_MYELOID_LEUKOCYTE_DIFFERENTIATION", 
#                       "GOBP_MITOTIC_NUCLEAR_DIVISION", "GOBP_RECOMBINATIONAL_REPAIR"
#))

mat_for_heatmap <- mat_for_heatmap[which(rownames(mat_for_heatmap) %in% GO_to_keep), ] 
mat_for_heatmap <- mat_for_heatmap[order(mat_for_heatmap[, 1], decreasing = T), ]
mat_for_heatmap <- mat_for_heatmap[order(mat_for_heatmap[, 2], decreasing = T), ]
mat_for_heatmap <- mat_for_heatmap[order(mat_for_heatmap[, 3], decreasing = T), ]
mat_for_heatmap <- mat_for_heatmap[GO_to_keep, ]
Heatmap(mat_for_heatmap[, 1:3], show_row_names = T, cluster_rows = F, cluster_columns = F)

library(reshape2)
ht <- reshape2::melt(mat_for_heatmap,value.name="Exp",na.rm = F)
colnames(ht)[1:2] <- c("gene","sample")
head(ht, 12)

# Change to factor
ht$gene <- factor(ht$gene,
                  levels = rev(unique(ht$gene)),
                  ordered = T)
ht$sample <- factor(ht$sample,
                    levels = unique(ht$sample),
                    ordered = T)

cols <- colorRampPalette(c("#fefefe","#ec66a0", rep("#4c2467", 3)))(200)
p <- ggplot(data = ht,
            aes(x = sample, y= gene))+
  geom_tile(aes(fill = Exp, height = 1, width = 1)) + 
  theme(axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 7), 
        axis.title.y = element_blank())+
  labs(x = NULL)+
  scale_fill_gradientn(colours = cols) +
  guides(fill=guide_colorbar(barheight = 4))+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 30, 
                                   size = 7, 
                                   hjust = 1, 
                                   vjust = 1.2), 
        axis.text.y = element_text(size = 7),
        axis.text.y.right = element_text(size = 7), 
        legend.position = "none", 
        panel.border = element_rect(fill=NA, color="black", linewidth = 1, linetype="solid")
  )
p

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Mat_HM_GO_enrichment_of_dif_phase_momac", ".pdf"), 
    width = 4.6, # The width of the plot in inches
    height = 8) # The height of the plot in inches
print(p, newpage = FALSE)
dev.off()

#### Cell density plots  ----- 
library(dplyr)
library(tidyr)
pseudotime <- pseudotime(cds) %>% as.data.frame()
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "pseudotime"
celltype <- cds@colData$group %>% as.data.frame()
celltype$cell <- colnames(cds)
colnames(celltype)[1] <- "clusters"
merge <-merge(pseudotime, celltype, by = 'cell')
merge <- merge[order(merge$pseudotime), ]
#to_keep <- unique(merge$clusters)[!unique(merge$clusters) %in% c("Mono_S100A8", "MoMac_IL1B", "Macro_FABP5", "Macro_FOLR2")]
#to_keep <- unique(merge$clusters)[unique(merge$clusters) %in% c("Mono_S100A8", "MoMac_IL1B", "Macro_FABP5")]

library(ggridges)
#class(merge$pseudotime)
merge$clusters <- factor(merge$clusters)
merge_na_omit <- na.omit(merge)
class(merge_na_omit$pseudotime) <- "numeric"
ggplot(merge_na_omit, aes(x = pseudotime, y = clusters#, fill=clusters
)) + 
  geom_density_ridges() +
  scale_y_discrete(position = 'right') +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8))+
  scale_x_continuous(position = "bottom")

p1 <- ggplot(merge_na_omit, aes(x = pseudotime, y = clusters, fill = clusters)) + 
  geom_boxplot(outlier.size = 0, outlier.shape = NA) +   
  #scale_fill_manual(values = mycol[c("Macro_FOLR2", "Macro_CCT3", "Macro_FABP5", "Macro_KLF2", 
  #                                   "MoMac_CLEC5A", "Macro_PCNA", "MoMac_IL1B", "Macro_ITGB2-", "Macro_IFIT1", 
  #                                   "Macro_FN1", "Macro_MT1G", "Mono_S100A8", "Macro_LYVE1")]) + 
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'), 
        axis.text.x =element_text(size = 12, color = 'black')) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.background = element_rect(fill = "#EAEAEA"), 
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)
  ) 
p1

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_boxplot_selected_Myeloids.pdf"), 
    width = 3, # The width of the plot in inches
    height = 7) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()
          
#### Plot genes -----
#### Pseudotime but 4-group partition
peuGene_exp_H_BF <- peuGene_exp[,grep(pattern="^H_BF",colnames(peuGene_exp))] 
pseudotime_H_BF <- pseudotime[colnames(peuGene_exp_H_BF),]
peuGene_exp_L_BF <- peuGene_exp[,grep(pattern="^L_BF",colnames(peuGene_exp))]
pseudotime_L_BF <- pseudotime[colnames(peuGene_exp_L_BF),]
peuGene_exp_H_AF <- peuGene_exp[,grep(pattern="^H_AF",colnames(peuGene_exp))] 
pseudotime_H_AF <- pseudotime[colnames(peuGene_exp_H_AF),]
peuGene_exp_L_AF <- peuGene_exp[,grep(pattern="^L_AF",colnames(peuGene_exp))]
pseudotime_L_AF <- pseudotime[colnames(peuGene_exp_L_AF),]

scaleFUN <- function(x) sprintf("%.2f", x)
pplist = list()
for (i in c("IL1B", "IL1RN", "CXCL2", "CXCL8", "CCL3", "FOS", "ATF3", "NFKBIA", 
            "NFKBIZ", "CD83", "CD300E", "SLAMF9", "FCN1", "S100A8", "S100A9", "BCL2A1", "MCL1", 
            "NR4A1", "DUSP1", "NR4A2", "NR4A3", 
            "CD14", "CD163", "MRC1", "MARCO", "MSR1", "MS4A4A", "IDO1", "TREM2", "CCL18", 
            "HLA-DRA", "CD74", "FCGR2A", "CYBB", "CXCL10", "CXCL11", "ISG15", "DAB2", 
            "VSIG4", "STAB1", "GPNMB", "SLCO2B1", "CSF1R", "FABP5", "APOC1", "FOLR2", "TOP2A", 
            "TK1", "CENPF", "UBE2C")) {
  j <- peuGene_exp_H_BF[i,]
  j <- as.data.frame(t(j))
  j$pseudotime <- pseudotime_H_BF
  j$group <- 'H_Pre'
  k <- peuGene_exp_H_AF[i,]
  k <- as.data.frame(t(k))
  k$pseudotime <- pseudotime_H_AF
  k$group <- 'H_Post'
  l <- peuGene_exp_L_BF[i,]
  l <- as.data.frame(t(l))
  l$pseudotime <- pseudotime_L_BF
  l$group <- 'L_Pre'
  m <- peuGene_exp_L_AF[i,]
  m <- as.data.frame(t(m))
  m$pseudotime <- pseudotime_L_AF
  m$group <- 'L_Post'
  
  peu_trend <- rbind(j, k, l, m)
  p <- ggplot(peu_trend, aes_string(x = "pseudotime", y = i, color = "group"))+
    geom_smooth(aes(fill=group)) +  
    xlab('pseudotime') +
    ylab('Relative expression') +
    #ggtitle(i) +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(color = 'black',size = 12),
          axis.title = element_text(color = 'black',size = 14))+
    scale_color_manual(name=NULL, values = c("#D9A0B3", "#54426D", "#6B798E", "#0f5688")) +
    scale_fill_manual(name=NULL, values = c("#D9A0B3", "#54426D", "#6B798E","#0f5688")) #+ 
  scale_x_continuous(limits = c(0, 35)) 
  
  pp1 <- p + annotate("text", x = range(p$data$pseudotime)[2]*0.8, y = range(p$data[i])[2]*0.8, size = 5,
                      label = i,
                      fontface="italic",
                      colour="black")  + 
    scale_y_continuous(labels=scaleFUN) +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
          plot.title = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(), 
          strip.background = element_blank(), 
          strip.text.y = element_text(size = 8, angle = 0, hjust = 0), 
          legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8), 
          legend.position = "none") 
  pplist[[i]] <- pp1 
}

library(cowplot)
pp3 <- plot_grid(pplist[[1]], pplist[[2]], pplist[[3]], pplist[[4]], pplist[[5]], pplist[[6]], 
                 pplist[[7]], pplist[[8]], pplist[[10]], pplist[[11]], pplist[[13]], pplist[[14]], 
                 pplist[[15]], pplist[[16]], pplist[[18]], pplist[[19]], pplist[[22]], 
                 pplist[[23]], pplist[[25]], pplist[[27]], pplist[[29]], pplist[[32]], pplist[[39]], 
                 pplist[[40]], pplist[[41]], pplist[[42]], pplist[[43]], pplist[[44]], pplist[[45]], 
                 pplist[[47]], pplist[[48]], pplist[[49]], ncol = 8)

pp3
time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "_Pseudotime_gene_trend_Myeloid.pdf"), 
    width = 24, # The width of the plot in inches
    height = 12) # The height of the plot in inches
print(pp3, newpage = FALSE)
dev.off()

## cDC2 traj (Figure 6H) -----
scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Myeloid/Myeloid_final_as_at_24_07.rds"))

library(monocle3)
unique(Idents(scrna)) # Selected subsets only 
scrna_sbuset <- subset(scrna, idents = c("cDC2_CD1C+"))
#scrna_sbuset <- subset(scrna, idents = c("cDC1_CLEC9A+", "cDC2_CD1C+", "Mature_DC"))
#"Mono_FCN1+", "Macro_CCL4+", "Macro_ISGs+", "Macro_MTs+", #"Macro_LYVE1+", 
#"Macro_APOE+", 
#"Macro_Ki67+", "Promyelocytes", "Neu", "cDC1_CLEC9A+", "cDC2_CD1C+", "Mature_DC"

data <- GetAssayData(scrna_sbuset, assay = 'RNA', layer = 'counts')
cell_metadata <- scrna_sbuset@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(expression_dat = data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds) # Preprocessing data 
cds <- align_cds(cds, alignment_group = "sample") # Remove batch effects
cds <- reduce_dimension(cds) # Dimension reduction
cds <- cluster_cells(cds) # Cluster cells 
cds <- learn_graph(cds) # Predict trajectory 

time_bin = "C"
cell_ids <- which(colData(cds)[, "group"] == time_bin)
closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                            (which.max(table(closest_vertex[cell_ids,]))))]

#cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
#cds <- order_cells(cds, root_pr_nodes = "Y_70")
cds <- order_cells(cds, root_pr_nodes = c("Y_2","Y_76"))
plot_cells(cds, color_cells_by = "pseudotime", 
           label_branch_points = F,
           label_roots = TRUE,
           label_leaves = F)

#cds <- order_cells(cds, root_pr_nodes='Y_21') # Set root nodes  
#cds <- order_cells(cds)
plot_cells(cds, label_cell_groups = T)
plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by="treatment_status")
plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by="Cell_type_fine_harmony", label_branch_points = F)
plot_cells(cds, color_cells_by = "pseudotime")

p1 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5, 
                 cell_size = 1,
                 trajectory_graph_color = "#CCCCCC", 
                 trajectory_graph_segment_size = 1)
p2 <- p1 + theme(legend.position = "none",
                 axis.ticks.x = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(), 
                 axis.title.y = element_blank(), 
                 axis.title.x = element_blank(), 
                 axis.line.x = element_blank(), 
                 axis.line.y = element_blank(), 
                 panel.border = element_blank()) +
  scale_color_gradientn(colours = c("#e9f4f6", "#8f9fab", "#1e2235"))
p2

time <- gsub(" ", "_", Sys.time()) 
time <- gsub("-", "_", time) 
time <- gsub(":", "_", time) 
pdf(file = paste0(time, "Pseudotime_cDC2.pdf"), 
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(p2, newpage = FALSE)
dev.off()

#### Graph test  ----- 
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 8)
genes <- row.names(subset(modulated_genes, q_value < 0.5 & morans_I > 0.2)) # Relaxed criteria 

#install.packages("ClusterGVis")
library(ClusterGVis)

#### Heatplot of key genes (w complext heatmap) ------ 
library(dplyr)
library(tidyr)
pseudotime <- pseudotime(cds) %>% as.data.frame()
pseudotime$cell <- rownames(pseudotime)
colnames(pseudotime)[1] <- "pseudotime"
celltype <- cds@colData$group %>% as.data.frame()
celltype$cell <- colnames(cds)
colnames(celltype)[1] <- "clusters"
merge <-merge(pseudotime, celltype, by = 'cell')

#devtools::install_github("junjunlab/ClusterGVis", force = T)
exp <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)

exp <- data.frame(t(exprs(cds)))
exp <- exp[, c("FOS", 
               "JUN", 
               "ATF3",
               "CD83", 
               "NR4A1", 
               "NR4A2", 
               "NR4A3",
               "IL1B",
               "BCL2A1"
)
]

heatmap_data_scaled <- data.frame(scale(exp))

library(ComplexHeatmap)
library(circlize)
heatmap_data_scaled$cell <- rownames(heatmap_data_scaled)
merge_mat <- merge(merge, heatmap_data_scaled, by = "cell")
merge_mat <- merge_mat[order(merge_mat$pseudotime), ]
merge_mat_for_plot <- as.data.frame(t(merge_mat)[c(4:12), ])
for (i in 1:ncol(merge_mat_for_plot)) {
  merge_mat_for_plot[, i] <- as.numeric(merge_mat_for_plot[, i])
}

library(circlize)
col_fun = colorRamp2(c(0, 7.5, 15), c("#e9f4f6", "#8f9fab", "#1e2235"))

ht_list = HeatmapAnnotation(
  Pseudotime = anno_barplot(merge_mat$pseudotime, 
                            gp = gpar(fill = col_fun(merge_mat$pseudotime), col = col_fun(merge_mat$pseudotime)), 
                            height = unit(2.1, "cm")), 
  Group = merge_mat$clusters, 
  col = list(Group = c("A" = "#54426D", "B" = "#D9A0B3", "C" = "#0f5688", D = "#6B798E"))
)

cellheight = 0.7
rn = dim(as.matrix(merge_mat_for_plot))[1]
h=cellheight*rn
pdf(file = #"DC_pseudotime_gene_HM.pdf", 
      "DC_pseudotime_gene_HM_260226.pdf",
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
Heatmap(as.matrix(#merge_mat[c("ATF3", "BCL2A1", "CD83", "FOS", "IL1B", "JUN", "NR4A1", "NR4A2", "NR4A3"), ]
  merge_mat_for_plot),
  # Height of the blocks
  #width = unit(w, "cm"),
  height = unit(h, "cm"),
  #scale = F,
  #rect_gp = gpar(col = "white", lwd = 1.5),
  #border_g = gpar(col = ,lty = 1,lwd = 1.2),
  # Dent formatting 
  #column_dend_height = unit(1.5, "cm"), 
  #row_dend_width = unit(1.5, "cm"),
  #column_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
  #row_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
  # Gap numbers and style 
  #row_split = 6, column_split = c("1", "1", "2", "2"),
  #row_gap = unit(2, "mm"),
  #column_gap = unit(2, "mm"),
  # Text formats 
  row_title = NULL,column_title = NULL,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  show_heatmap_legend = T, 
  cluster_rows = F,
  cluster_columns = F,
  row_names_side = 'right',
  show_column_names = F,
  show_row_names = T,
  border = T, 
  top_annotation = ht_list, 
  col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#50859f","white", "#d66692"))
)
dev.off()

## Code related to Figure S6A ----- 
## Ti and Pi -----
library(reshape2)
Idents(scrna) <- "reportname"
scrna <- RenameIdents(scrna,
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
scrna$TNR_val <- Idents(scrna)

# Pi 
Idents(scrna) <- "treatment_status"
scrna_pre <- subset(scrna, idents = "Pre")
Cellratio <- prop.table(table(scrna_pre$Cell_type_fine_harmony, scrna_pre$sample), margin = 2)  
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio, Var2~Var1, value.var="Freq") 
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]
cellper <- na.omit(cellper)

meta <- scrna_pre@meta.data
colnames(meta)
meta <- meta[,c(4, 47)] #group and TNR
meta <- as.data.frame(meta)
cellper$sample <- rownames(cellper)
#cellper$group <- NA
meta <- unique(meta) 
cellper <- merge(cellper, meta, by = "sample")
colnames(cellper)
#colnames(cellper)[18] <- "group"
cellper <- as.data.frame(cellper)
#cellper$patient <- str_replace(cellper$sample, "_BF_", "")

# Calc 
library(ggplot2)
library(dplyr)
cellper$TNR_val <- as.numeric(as.character(cellper$TNR_val))

#for (i in c(2:length(colnames(cellper)))) {
#  class(cellper[, i]) <- "numeric"
#} # Change to numeric 
#standardize = function(x){
#  z <- (x - mean(x)) / sd(x)
#  return( z)
#} # Define stardardise function
#cellper$TNR_val_norm <- apply(cellper["TNR_val"], 2, standardize) #  Create TNR_val_norm col
#for (i in as.vector(unique(scrna$Cell_type_fine_harmony))) {
#  cellper[i] <- apply(cellper[i], 2, standardize)
#} # Normalise cell percentages 

# Create the matrix for plot
matrix_for_Plots <- matrix(data = rep(NA, 3*length(unique(scrna$Cell_type_fine_harmony))), 
                           ncol = 3, nrow = length(unique(scrna$Cell_type_fine_harmony)), dimnames = list(c(as.vector(unique(scrna$Cell_type_fine_harmony))), c("Cell_type", "Pi", "pval"))) %>% as.data.frame()
matrix_for_Plots$Cell_type <- rownames(matrix_for_Plots)

for (i in as.vector(unique(scrna$Cell_type_fine_harmony))) {
  m <- lm(unlist(cellper["TNR_val"]) ~ unlist(cellper[i]))
  b <- as.numeric(round(unname(coef(m)[2]), digits = 2))
  r2 = as.numeric(round(summary(m)$r.squared, digits = 2))
  Pi <- b/(abs(b)) * r2
  pval <- summary(m)$coefficients
  pval <- as.numeric(round(pval[2, 4], digits = 3))
  matrix_for_Plots[i, "Pi"] <- Pi
  matrix_for_Plots[i, "pval"] <- pval 
  matrix_for_Plots$log10pval <- -log10(matrix_for_Plots$pval)
}

matrix_for_Plots <- matrix_for_Plots %>%
  arrange(Pi) %>%
  mutate(Cell_type = factor(Cell_type,levels = Cell_type)) 
matrix_for_Plots <- na.omit(matrix_for_Plots)

# Plot
#mycol <- c("#92699e", "#a5a9b0", "#696a6c", "#d4c2db", "#b95055", "#6d6fa0", "#c49abc", "#d8a0c0",
#"#63a3b8", "#9a70a8", "#4490c4", "#a44e89", "#927c9a", "#de8b36", "#9f8d89", "#408444")

mycol <- celltype_colors
names(mycol) <- levels(scrna$Cell_type_fine_harmony)
mycol <- mycol[levels(matrix_for_Plots$Cell_type)]

#bg_colors = c("white", "#F6F9E4")
#gradient_grob <- rasterGrob(colorRampPalette(bg_colors)(256), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
p <- ggplot(matrix_for_Plots, aes(x = Cell_type, y = Pi)) +  
  #annotation_custom(gradient_grob,xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_segment(aes(x = Cell_type, xend = Cell_type, y = 0, yend = Pi, color = Cell_type),                 
               linetype = "solid", size = 1, color = mycol) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour="gray80") +  
  geom_point(aes(color = Cell_type, size = log10pval), color = mycol
  ) + 
  scale_size_continuous(range=c(3,15)) +
  #geom_text(aes(label = Ti), size = 3) +  
  geom_text(aes(label = Cell_type, y = 0), 
            hjust = ifelse(matrix_for_Plots$Pi >= 0, 1, 0), 
            #vjust = -0.5, 
            angle = 90, 
            fontface = 'italic',
            color = ifelse(matrix_for_Plots$pval >= 0.05, "black", 'red'),             
            size = 4) +    
  theme_classic() +
  theme(axis.text = element_text(size =10),axis.title = element_text(size = 10),legend.text = element_text(size =10),
        legend.title = element_text(size =10),plot.title = element_text(size =10,face ='plain'),legend.position ='none')+
  labs(title = "Predicative Index", y = "Predicative Index") + 
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18, hjust = 0.5),
        legend.position = "none") + 
  NoLegend() + 
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid")) 
p

# Ti 
# Needs cell per change 
library(reshape2)
# Calc post
Idents(scrna) <- "treatment_status"
scrna_post <- subset(scrna, idents = "Post")
Cellratio_post <- prop.table(table(scrna_post$Cell_type_fine_harmony, scrna_post$sample), margin =2) 
Cellratio_post <- data.frame(Cellratio_post)

cellper_post <- dcast(Cellratio_post, Var2~Var1, value.var="Freq") 
rownames(cellper_post)<- cellper_post[,1]
cellper_post <- cellper_post[,-1]
cellper_post <- na.omit(cellper_post)

# Calc pre
#Idents(scrna) <- "treatment_status"
#scrna_post <- subset(scrna, idents = "Pre")
Cellratio_pre <- prop.table(table(scrna_pre$Cell_type_fine_harmony, scrna_pre$sample), margin =2) 
Cellratio_pre <- data.frame(Cellratio_pre)

cellper_pre <- dcast(Cellratio_pre, Var2~Var1, value.var="Freq") 
rownames(cellper_pre)<- cellper_pre[,1]
cellper_pre <- cellper_pre[,-1]
cellper_pre <- na.omit(cellper_pre)
#cellper_pre <- log2(cellper_pre + 1)

cellper_delta <- cellper_post - cellper_pre

meta <- scrna_post@meta.data
colnames(meta)
meta <- meta[,c(4,47)] #group and TNR
meta <- as.data.frame(meta)
cellper_delta$sample <- rownames(cellper_delta)
#cellper_post$group <- NA
meta <- unique(meta)
cellper_delta <- merge(cellper_delta, meta, by = "sample")
colnames(cellper_delta)
#colnames(cellper_post)[18] <- "group"
cellper_delta <- as.data.frame(cellper_delta)
#cellper_delta$patient <- str_replace(cellper_delta$sample, "_AF_", "")

# Calc 
library(ggplot2)
library(dplyr)
cellper_delta$TNR_val <- as.numeric(as.character(cellper_delta$TNR_val))

#for (i in c(2:length(colnames(cellper_delta)))) {
#  class(cellper_delta[, i]) <- "numeric"
#} # Change to numeric 
#standardize = function(x){
#  z <- (x - mean(x)) / sd(x)
#  return( z)
#} # Define stardardise function
#cellper_delta$TNR_val_norm <- apply(cellper_delta["TNR_val"], 2, standardize) #  Create TNR_val_norm col
#for (i in as.vector(unique(scrna$Cell_type_fine_harmony))) {
#  cellper_delta[i] <- apply(cellper_delta[i], 2, standardize)
#} # Normalise cell percentages 

# Create the matrix for plot
matrix_for_Plots_delta <- matrix(data = rep(NA, 3*length(unique(scrna$Cell_type_fine_harmony))), 
                                 ncol = 3, nrow = length(unique(scrna$Cell_type_fine_harmony)), dimnames = list(c(as.vector(unique(scrna$Cell_type_fine_harmony))), c("Cell_type", "Ti", "pval"))) %>% as.data.frame()
matrix_for_Plots_delta$Cell_type <- rownames(matrix_for_Plots_delta)

for (i in as.vector(unique(scrna$Cell_type_fine_harmony))) {
  m <- lm(unlist(cellper_delta["TNR_val"]) ~ unlist(cellper_delta[i]))
  b <- as.numeric(round(unname(coef(m)[2]), digits = 2))
  r2 = as.numeric(round(summary(m)$r.squared, digits = 2))
  Ti <- b/(abs(b)) * r2
  pval <- summary(m)$coefficients
  pval <- as.numeric(round(pval[2, 4], digits = 3))
  matrix_for_Plots_delta[i, "Ti"] <- Ti
  matrix_for_Plots_delta[i, "pval"] <- pval 
  matrix_for_Plots_delta$log10pval <- -log10(matrix_for_Plots_delta$pval)
}

matrix_for_Plots_delta <- matrix_for_Plots_delta %>%
  arrange(Ti) %>%
  mutate(Cell_type = factor(Cell_type,levels = Cell_type)) 
matrix_for_Plots_delta <- na.omit(matrix_for_Plots_delta)
#c4a_gui()

mycol1 <- as.data.frame(mycol); mycol1$celltype <- rownames(mycol1); 
mycol1 <- mycol1[match(matrix_for_Plots_delta$Cell_type, mycol1$celltype), ]       
mycol <- mycol1$mycol

#bg_colors = c("white", "#E7CAD3")
#gradient_grob <- rasterGrob(colorRampPalette(bg_colors)(256), width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
p1 <- ggplot(matrix_for_Plots_delta, aes(x = Cell_type, y = Ti)) +  
  #annotation_custom(gradient_grob,xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_segment(aes(x = Cell_type, xend = Cell_type, y = 0, yend = Ti, color = Cell_type),                 
               linetype = "solid", size = 1, color = mycol) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour="gray80") +  
  geom_point(aes(color = Cell_type, size = log10pval), color = mycol
  ) + 
  scale_size_continuous(range=c(3,15)) +
  #geom_text(aes(label = Ti), size = 3) +  
  geom_text(aes(label = Cell_type, y = 0), 
            hjust = ifelse(matrix_for_Plots_delta$Ti >= 0, 1, 0), 
            #vjust = -0.5, 
            angle = 90, 
            fontface = 'italic',
            color = ifelse(matrix_for_Plots_delta$pval >= 0.05, "black", 'red'),             
            size = 4) + 
  theme_classic() +
  theme(axis.text = element_text(size =10),axis.title = element_text(size = 10),legend.text = element_text(size =10),
        legend.title = element_text(size =10),plot.title = element_text(size =10,face ='plain'),legend.position ='none')+
  labs(title = "Therapeutic Index", y = "Therapeutic Index") + 
  theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
        plot.title = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18, hjust = 0.5),
        legend.position = "none") + 
  NoLegend() + 
  theme(panel.border = element_rect(fill = NA, color="black", linewidth = 1, linetype = "solid")) 
scaleFUN <- function(x) sprintf("%.2f", x)
p1 <- p1 + scale_y_continuous(labels=scaleFUN)


library(cowplot)
p2 <-plot_grid(p, p1, ncol = 1)
p2

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Pi_Ti_all_cells_for_fig.pdf"), 
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
print(p2, newpage = FALSE)
dev.off()

rm(list = ls())

## InfoPlots -----
Idents(scrna) <- "group"
scrna <- RenameIdents(scrna, "A" = "Pre", 
                      "B" = "Post", 
                      "C" = "Pre", 
                      "D" = "Post")
scrna$treatment_status
Prog_calc <- function (object) {
  Idents(object) <- "treatment_status"
  ##### Pi #####
  object_pre <- subset(object, idents = "Pre")
  object_pre$Cell_type_fine_harmony <- droplevels(object_pre$Cell_type_fine_harmony)
  Cellratio <- prop.table(table(object_pre$Cell_type_fine_harmony, object_pre$sample), margin = 2)  
  Cellratio <- data.frame(Cellratio)
  
  cellper <- reshape2::dcast(Cellratio, Var2~Var1, value.var="Freq") 
  rownames(cellper)<- cellper[,1]
  cellper <- cellper[,-1]
  cellper <- na.omit(cellper)
  
  meta <- object_pre@meta.data
  colnames(meta)
  meta <- meta[,c("sample", "TNR_val")] 
  meta <- as.data.frame(meta)
  cellper$sample <- rownames(cellper)
  #cellper$group <- NA
  meta <- unique(meta) 
  cellper <- merge(cellper, meta, by = "sample")
  #colnames(cellper)
  #colnames(cellper)[18] <- "group"
  cellper <- as.data.frame(cellper)
  #cellper$patient <- str_replace(cellper$sample, "_BF_", "")
  
  ##### Calc #####
  library(ggplot2)
  library(dplyr)
  cellper$TNR_val <- as.numeric(as.character(cellper$TNR_val))
  
  for (i in c(2:length(colnames(cellper)))) {
    class(cellper[, i]) <- "numeric"
  } # Change to numeric 
  standardize = function(x){
    z <- (x - mean(x)) / sd(x)
    return( z)
  } # Define stardardise function
  
  cellper$TNR_val_norm <- apply(cellper["TNR_val"], 2, standardize) #  Create TNR_val_norm col
  for (i in as.vector(unique(object$Cell_type_fine_harmony))) {
    cellper[i] <- apply(cellper[i], 2, standardize)
  } # Normalise cell percentages 
  
  ##### Create the matrix for plot#####
  matrix_for_Plots <- matrix(data = rep(NA, 3*length(unique(object$Cell_type_fine_harmony))), 
                             ncol = 3, nrow = length(unique(object$Cell_type_fine_harmony)), dimnames = list(c(as.vector(unique(object$Cell_type_fine_harmony))), c("Cell_type", "Pi", "pval"))) %>% as.data.frame()
  matrix_for_Plots$Cell_type <- rownames(matrix_for_Plots)
  
  for (i in as.vector(unique(object$Cell_type_fine_harmony))) {
    m <- lm(unlist(cellper["TNR_val_norm"]) ~ unlist(cellper[i]))
    b <- as.numeric(round(unname(coef(m)[2]), digits = 2))
    r2 = as.numeric(round(summary(m)$r.squared, digits = 2))
    Pi <- b/(abs(b)) * r2
    pval <- summary(m)$coefficients
    pval <- as.numeric(round(pval[2, 4], digits = 3))
    matrix_for_Plots[i, "Pi"] <- Pi
    matrix_for_Plots[i, "pval"] <- pval 
    matrix_for_Plots$log10pval <- -log10(matrix_for_Plots$pval)
  }
  
  Pi_temp <- matrix_for_Plots
  rm(matrix_for_Plots)
  
  ##### Ti #####
  ### Needs cell per change 
  # Calc post
  Idents(object) <- "treatment_status"
  object_post <- subset(object, idents = "Post")
  object_post$Cell_type_fine_harmony <- droplevels(object_post$Cell_type_fine_harmony)
  Cellratio_post <- prop.table(table(object_post$Cell_type_fine_harmony, object_post$sample), margin =2) 
  Cellratio_post <- data.frame(Cellratio_post)
  
  cellper_post <- reshape2::dcast(Cellratio_post, Var2~Var1, value.var="Freq") 
  rownames(cellper_post)<- cellper_post[,1]
  cellper_post <- cellper_post[,-1]
  cellper_post <- na.omit(cellper_post)
  
  # Calc pre
  #Idents(object) <- "treatment_status"
  #object_post <- subset(object, idents = "Pre")
  Cellratio_pre <- prop.table(table(object_pre$Cell_type_fine_harmony, object_pre$sample), margin =2) 
  Cellratio_pre <- data.frame(Cellratio_pre)
  
  cellper_pre <- reshape2::dcast(Cellratio_pre, Var2~Var1, value.var="Freq") 
  rownames(cellper_pre)<- cellper_pre[,1]
  cellper_pre <- cellper_pre[,-1]
  cellper_pre <- na.omit(cellper_pre)
  #cellper_pre <- log2(cellper_pre + 1)
  
  cellper_delta <- cellper_post - cellper_pre
  
  meta <- object_post@meta.data
  colnames(meta)
  meta <- meta[,c("sample","TNR_val")]#group和orig.ident信息
  meta <- as.data.frame(meta)
  cellper_delta$sample <- rownames(cellper_delta)
  #cellper_post$group <- NA
  meta <- unique(meta)
  cellper_delta <- merge(cellper_delta, meta, by = "sample")
  colnames(cellper_delta)
  #colnames(cellper_post)[18] <- "group"
  cellper_delta <- as.data.frame(cellper_delta)
  #cellper_delta$patient <- str_replace(cellper_delta$sample, "_AF_", "")
  
  ##### Calc #####
  library(ggplot2)
  library(dplyr)
  cellper_delta$TNR_val <- as.numeric(as.character(cellper_delta$TNR_val))
  #cellper_post[, as.vector(unique(object$Cell_type_fine))]
  
  for (i in c(2:length(colnames(cellper_delta)))) {
    class(cellper_delta[, i]) <- "numeric"
  } # Change to numeric 
  
  standardize = function(x){
    z <- (x - mean(x)) / sd(x)
    return( z)
  } # Define stardardise function
  
  cellper_delta$TNR_val_norm <- apply(cellper_delta["TNR_val"], 2, standardize) #  Create TNR_val_norm col
  for (i in as.vector(unique(object$Cell_type_fine_harmony))) {
    cellper_delta[i] <- apply(cellper_delta[i], 2, standardize)
  } # Normalise cell percentages 
  
  ##### Create the matrix for plot#####
  matrix_for_Plots_delta <- matrix(data = rep(NA, 3*length(unique(object$Cell_type_fine_harmony))), 
                                   ncol = 3, nrow = length(unique(object$Cell_type_fine_harmony)), dimnames = list(c(as.vector(unique(object$Cell_type_fine_harmony))), c("Cell_type", "Ti", "pval"))) %>% as.data.frame()
  matrix_for_Plots_delta$Cell_type <- rownames(matrix_for_Plots_delta)
  
  for (i in as.vector(unique(object$Cell_type_fine_harmony))) {
    m <- lm(unlist(cellper_delta["TNR_val_norm"]) ~ unlist(cellper_delta[i]))
    b <- as.numeric(round(unname(coef(m)[2]), digits = 2))
    r2 = as.numeric(round(summary(m)$r.squared, digits = 2))
    Ti <- b/(abs(b)) * r2
    pval <- summary(m)$coefficients
    pval <- as.numeric(round(pval[2, 4], digits = 3))
    matrix_for_Plots_delta[i, "Ti"] <- Ti
    matrix_for_Plots_delta[i, "pval"] <- pval 
    matrix_for_Plots_delta$log10pval <- -log10(matrix_for_Plots_delta$pval)
  }
  
  matrix_for_Plots_delta
  Ti_temp <- matrix_for_Plots_delta
  colnames(Ti_temp) <- paste0(colnames(Ti_temp), "_for_Ti")
  Ti_temp$Cell_type <- rownames(Ti_temp)
  Pi_Ti_Matrix <- merge(Pi_temp, Ti_temp, by = "Cell_type")
  return(Pi_Ti_Matrix)
} # Remember to add "TNR_val" before calculation
Pi_Ti_Matrix <- Prog_calc(scrna)
Pi_Ti_Matrix$group_info <- ifelse(Pi_Ti_Matrix$Ti_for_Ti > 0.2, "Overwhelmingly favourable", 
                                  ifelse(Pi_Ti_Matrix$Ti_for_Ti < -0.2, "Overwhelmingly unfavourable", 
                                         ifelse(Pi_Ti_Matrix$Ti_for_Ti <= 0.2 & Pi_Ti_Matrix$Ti_for_Ti > 0.05, "Favourable",
                                                ifelse(Pi_Ti_Matrix$Ti_for_Ti < -0.05 & Pi_Ti_Matrix$Ti_for_Ti >= -0.2, "Unfavourable", "Indeterminate"))))
## Info plot -----
Matrix_for_plot_all_cells <- data.frame(Expression = scrna$CytoTRACE2_Relative, celltype = scrna$Cell_type_fine_harmony)
#Matrix_for_plot_all_cells$celltype <- factor(Matrix_for_plot_all_cells$celltype, levels = rev(c("GMP", "Macro_MKI67", "Macro_PCNA", "Macro_CCT3", 
#                                                                                                "Macro_LYVE1", "Macro_FN1", "Macro_IFIT1", "Macro_MT1G", "Macro_FABP5", "Macro_FOLR2", "Macro_KLF2", "Macro_ITGB2-",
#                                                                                                "Mono_S100A8", "MoMac_IL1B", "MoMac_CLEC5A", "MoMac_CXCL10", 
#                                                                                                "cDC1", "cDC2_FCGR2B", "cDC2_NR4A3", "cDC2_CD69", "mDC", 
#                                                                                                "Neu")))
#Matrix_for_plot_all_cells$celltype <- factor(Matrix_for_plot_all_cells$celltype, levels = c("", ))
RidgePlot <- ggplot(Matrix_for_plot_all_cells, aes(x = Expression, y = celltype)) +  
  ggridges::geom_density_ridges(aes(fill = celltype), alpha = 0.8, show.legend = FALSE, 
                                rel_min_height = 0.01, 
                                #rel_min_height = 0.01, #尾部修剪，数值越大修剪程度越高
                                scale = 1, 
                                quantile_lines = F,
                                #quantiles = 0.5, 
                                color = "white") +
  scale_x_continuous(limits = c(0.05, 0.75),
                     breaks = seq(0.05, 0.75, by = 0.2)) +
  #scale_fill_viridis_c(name = "CytoTRACE2_Score", option = "C")+
  coord_cartesian(clip = "off") +
  #scale_fill_manual(values = mycol2) +
  geom_vline(xintercept = c(0.05, 0.75),
             size = 0.5,
             color = 'grey',
             lty = 'dashed') + 
  scale_fill_manual(values = mycol) +
  theme_classic() + 
  #ggridges::theme_ridges(grid = F) +
  theme(panel.grid = element_blank(), #移除背景网格线
        panel.background = element_rect(fill = NULL, colour = "black", linewidth = 1), 
        axis.text.y = element_text(size = 8, vjust = -0.3), 
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1), 
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8), 
        #axis.text.y.left = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()
  ) 

RidgePlot 

Matrix_for_plot_all_cells$group_info <- NA
for (i in unique(Matrix_for_plot_all_cells$celltype)) {
  Matrix_for_plot_all_cells$group_info[which(Matrix_for_plot_all_cells$celltype == i)] <- Pi_Ti_Matrix$group_info[which(Pi_Ti_Matrix$Cell_type == i)]
}
Matrix_for_plot_all_cells$group_info <- factor(Matrix_for_plot_all_cells$group_info, levels = c("Overwhelmingly favourable", "Favourable", "Indeterminate", "Unfavourable", "Overwhelmingly unfavourable"), ordered = T)
HL_Group_info <- ggplot(Matrix_for_plot_all_cells, aes(x = 1, y = celltype)) +  
  geom_point(aes(fill = group_info), size = 4, shape = 21, 
             stroke = NA) +  
  scale_x_discrete(breaks = c(1)) +  
  scale_fill_manual(name = "", values = c('#4F6F46','#A9BE7B', '#B2B6B6', '#BA5B49', '#7C191E'),) +  
  theme_void() +  
  theme(legend.position = "top",   
        axis.text.y = element_text(size = 8, vjust = -0.3), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
HL_Group_info

### Log2fc for each -----
log2fc_calc_post_v_pre <- function(object) {
  object$Cell_type_fine_harmony <- droplevels(object$Cell_type_fine_harmony)
  Cellratio <- prop.table(table(object$Cell_type_fine_harmony, object$treatment_status), margin = 2)  
  Cellratio <- data.frame(Cellratio)
  
  cellper <- reshape2::dcast(Cellratio, Var2~Var1, value.var="Freq") 
  rownames(cellper)<- cellper[,1]
  cellper <- cellper[,-1]
  cellper <- na.omit(cellper) %>% t() %>% as.data.frame
  
  log2(cellper$Post / cellper$Pre)
  cellper$log2fc <- log2(cellper$Post / cellper$Pre)
  
  return(cellper)
}
cellper <- log2fc_calc_post_v_pre(scrna)
cellper$celltype <- rownames(cellper)
cellper$celltype  <- factor(cellper$celltype , levels = levels(scrna$Cell_type_fine_harmony))

log2fc_plot <- ggplot(cellper, aes(x = celltype, y = log2fc)) +  
  geom_col(fill = ifelse(cellper$log2fc > 0, "#BBA1CB", "#C9CFC1"), width = 0.5, color = "black", linewidth = 1) +  
  scale_y_continuous(limits = c(-6, 6),                      
                     breaks = c(-6, -4, -2, 0, 2, 4, 6),                     
                     expand = c(0, 0)) +  
  labs(y = "") +  coord_flip() +  theme_classic() +
  theme(legend.position = "none",        
        axis.title.x = element_text(size = 12),        
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),        
        axis.title.y = element_blank(),        
        axis.ticks.y = element_blank(),        
        axis.text.y = element_blank(),        
        axis.line.y = element_blank(),        
        plot.margin = unit(c(0, 1, 0, 0), "cm")) +  
  ylab(bquote(Log[2] ~ italic('Fold Change')))
log2fc_plot

## Updated Milo Log2fc plot -----
### Use for all cells partitioned by pre and post status ----- 
library(miloR)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
library(ggbeeswarm)
library(scater)
library(scales)
library(forcats)
library(data.table)
library(stringr)
library(dplyr)

# All cells with treatment status 
scrna_sce <- as.SingleCellExperiment(scrna)
scrna_sce_milo <- miloR::Milo(scrna_sce)  
scrna_sce_milo <- miloR::buildGraph(scrna_sce_milo, k = 30, d = 50)
scrna_sce_milo <- makeNhoods(scrna_sce_milo, 
                             prop = 0.2, 
                             k = 30, 
                             d = 50, 
                             refined = TRUE)
scrna_sce_milo <- countCells(scrna_sce_milo, meta.data = data.frame(colData(scrna_sce_milo)), 
                             sample="sample") 

traj_design <- data.frame(colData(scrna_sce_milo))[,c("sample", "treatment_status")]
traj_design$sample <- as.factor(traj_design$sample)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample
scrna_sce_milo <- calcNhoodDistance(scrna_sce_milo, d = 50)

da_results <- testNhoods(scrna_sce_milo, 
                         design = ~ treatment_status, 
                         design.df = traj_design)
scrna_sce_milo <- buildNhoodGraph(scrna_sce_milo)
da_results <- annotateNhoods(scrna_sce_milo, da_results, coldata_col = "Cell_type_fine_harmony")
da_results$Cell_type_fine_harmony <- factor(da_results$Cell_type_fine_harmony, levels = levels(scrna$Cell_type_fine_harmony))
p4 <- plotDAbeeswarm(da_results, group.by = "Cell_type_fine_harmony") +
  scale_color_gradient2(low="#50859f", 
                        mid="darkgrey",
                        high="#d66692",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black')) 
p5 <- p4 + #+ scale_y_continuous(limits = c(-6, 6),                      
  #                  breaks = c(-6, -4, -2, 0, 2, 4, 6),                     
  #                  expand = c(0, 0)) +  
  labs(y = "") +  #coord_flip() + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme(legend.position = "none",        
        axis.title.x = element_text(size = 12),        
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),        
        axis.title.y = element_blank(),        
        axis.ticks.y = element_blank(),        
        axis.text.y = element_blank(),        
        axis.line.y = element_blank(),        
        plot.margin = unit(c(0, 1, 0, 0), "cm")) +  
  ylab("Treatment status")
#ylab(bquote(Log[2] ~ italic('Fold Change')))
p5

saveRDS(da_results, "da_results_post_vs_pre_all_myeloid.rds")

### Use for all cells partitioned by efficacy status ----- 
scrna$efficacy <- scrna$orig.ident
scrna$efficacy <- factor(scrna$efficacy, levels = c("L", "H"))
scrna_sce <- as.SingleCellExperiment(scrna)
scrna_sce_milo <- miloR::Milo(scrna_sce)  
scrna_sce_milo <- miloR::buildGraph(scrna_sce_milo, k = 30, d = 50) 
scrna_sce_milo <- makeNhoods(scrna_sce_milo, 
                             prop = 0.2,  
                             k = 30,  
                             d = 50,   
                             refined = TRUE)
scrna_sce_milo <- countCells(scrna_sce_milo, meta.data = data.frame(colData(scrna_sce_milo)), 
                             sample="sample")   
  
traj_design <- data.frame(colData(scrna_sce_milo))[,c("sample", "efficacy")]  
traj_design$sample <- as.factor(traj_design$sample)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample
  
scrna_sce_milo <- calcNhoodDistance(scrna_sce_milo, d = 50)

da_results <- testNhoods(scrna_sce_milo, 
                         design = ~ efficacy, 
                         design.df = traj_design)  
scrna_sce_milo <- buildNhoodGraph(scrna_sce_milo)  
da_results <- annotateNhoods(scrna_sce_milo, da_results, coldata_col = "Cell_type_fine_harmony")
da_results$Cell_type_fine_harmony <- factor(da_results$Cell_type_fine_harmony, levels = levels(scrna$Cell_type_fine_harmony))
p6 <- plotDAbeeswarm(da_results, group.by = "Cell_type_fine_harmony") +
  scale_color_gradient2(low="#50859f", 
                        mid="darkgrey",
                        high="#d66692",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black')) 
p7 <- p6 + #+ scale_y_continuous(limits = c(-6, 6),                      
  #                  breaks = c(-6, -4, -2, 0, 2, 4, 6),                     
  #                  expand = c(0, 0)) +  
  labs(y = "") +  #coord_flip() + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme(legend.position = "none",        
        axis.title.x = element_text(size = 12),        
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),        
        axis.title.y = element_blank(),        
        axis.ticks.y = element_blank(),        
        axis.text.y = element_blank(),        
        axis.line.y = element_blank(),        
        plot.margin = unit(c(0, 1, 0, 0), "cm")) +  
  ylab("Response")
#ylab(bquote(Log[2] ~ italic('Fold Change')))
p7

saveRDS(da_results, "da_results_post_vs_pre_all_myeloid_by_HL.rds")

### Use for pre-treatment (baseline) cells partitioned by status ----- 
scrna$efficacy <- scrna$orig.ident
scrna$efficacy <- factor(scrna$efficacy, levels = c("L", "H"))
Idents(scrna) <- "treatment_status"
scrna_subset <- subset(scrna, idents = "Pre")
scrna_subset$sample <- droplevels(scrna_subset$sample)
scrna_subset$Cell_type_fine_harmony <- droplevels(scrna_subset$Cell_type_fine_harmony)
scrna_sce <- as.SingleCellExperiment(scrna_subset)#; rm(scrna_subset)
scrna_sce_milo <- miloR::Milo(scrna_sce)  
scrna_sce_milo <- miloR::buildGraph(scrna_sce_milo, k = 30, d = 50) 
  
scrna_sce_milo <- makeNhoods(scrna_sce_milo, 
                             prop = 0.2,  
                             k = 30,  
                             d = 50,   
                             refined = TRUE)
scrna_sce_milo <- countCells(scrna_sce_milo, meta.data = data.frame(colData(scrna_sce_milo)), 
                             sample="sample")   
  
traj_design <- data.frame(colData(scrna_sce_milo))[,c("sample", "efficacy")]  
traj_design$sample <- as.factor(traj_design$sample)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample
  
scrna_sce_milo <- calcNhoodDistance(scrna_sce_milo, d = 50)
  
da_results <- testNhoods(scrna_sce_milo, 
                         design = ~ efficacy, 
                         design.df = traj_design)  
scrna_sce_milo <- buildNhoodGraph(scrna_sce_milo)  
da_results <- annotateNhoods(scrna_sce_milo, da_results, coldata_col = "Cell_type_fine_harmony")
da_results$Cell_type_fine_harmony <- factor(da_results$Cell_type_fine_harmony, levels = levels(scrna_subset$Cell_type_fine_harmony))

##### Due to lack of datapoint for promye use below for plot (delete dummy datapoint when actually plotting in AI) ----- 
da.res <- da_results
group.by = "Cell_type_fine_harmony"
alpha = 0.1
da.res <- mutate(da.res, group_by = da.res[, group.by])
#Create a dummy datapoint for the promyelocyte row 
dim(da.res)
newrow <- da.res[8178, ]
rownames(newrow) <- "8179"
newrow$Cell_type_fine_harmony <- "Promyelocytes"
newrow$group_by <- "Promyelocytes"
da.res <- rbind(da.res, newrow)

beeswarm_pos <- ggplot_build(da.res %>% mutate(is_signif = ifelse(SpatialFDR < 
                                                                    alpha, 1, 0)) %>% arrange(group_by) %>% ggplot(aes(group_by, 
                                                                                                                       logFC)) + geom_quasirandom())
pos_x <- beeswarm_pos$data[[1]]$x
pos_y <- beeswarm_pos$data[[1]]$y
n_groups <- unique(da.res$group_by) %>% length()
#n_groups <- n_groups

p8 <- da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 
                                           1, 0)) %>% mutate(logFC_color = ifelse(is_signif == 1, 
                                                                                  logFC, NA)) %>% arrange(group_by) %>% mutate(Nhood = factor(Nhood, 
                                                                                                                                              levels = unique(Nhood))) %>% mutate(pos_x = pos_x, pos_y = pos_y) %>% 
  ggplot(aes(pos_x, pos_y, color = logFC_color)) + scale_color_gradient2() + 
  guides(color = "none") + xlab(group.by) + ylab("Log Fold Change") + 
  scale_x_continuous(breaks = seq(1, n_groups), labels = setNames(levels(da.res$group_by), 
                                                                  seq(1, n_groups))) + geom_point() + coord_flip() + 
  theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0)) +
  scale_color_gradient2(low="#50859f", 
                        mid="darkgrey",
                        high="#d66692",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10) +
  theme(axis.text = element_text(colour = 'black')) 

p9 <- p8 + #+ scale_y_continuous(limits = c(-6, 6),                      
  #                  breaks = c(-6, -4, -2, 0, 2, 4, 6),                     
  #                  expand = c(0, 0)) +  
  labs(y = "") +  #coord_flip() + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme(legend.position = "none",        
        axis.title.x = element_text(size = 12),        
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),        
        axis.title.y = element_blank(),        
        axis.ticks.y = element_blank(),        
        axis.text.y = element_blank(),        
        axis.line.y = element_blank(),        
        plot.margin = unit(c(0, 1, 0, 0), "cm")) +  
  ylab("Efficacy (Pre)")
#ylab(bquote(Log[2] ~ italic('Fold Change')))
p9

saveRDS(da_results, "da_results_pre_myeloid_by_HL.rds")

### Use for post-treatment (dynamic) cells partitioned by status ----- 
scrna$efficacy <- scrna$orig.ident
scrna$efficacy <- factor(scrna$efficacy, levels = c("L", "H"))
Idents(scrna) <- "treatment_status"
scrna_subset <- subset(scrna, idents = "Post")
scrna_subset$sample <- droplevels(scrna_subset$sample)
scrna_subset$Cell_type_fine_harmony <- droplevels(scrna_subset$Cell_type_fine_harmony)
scrna_sce <- as.SingleCellExperiment(scrna_subset)
scrna_sce_milo <- miloR::Milo(scrna_sce)  
scrna_sce_milo <- miloR::buildGraph(scrna_sce_milo, k = 30, d = 50) 
  
  
scrna_sce_milo <- makeNhoods(scrna_sce_milo, 
                             prop = 0.2,  
                             k = 30,  
                             d = 50,   
                             refined = TRUE)
scrna_sce_milo <- countCells(scrna_sce_milo, meta.data = data.frame(colData(scrna_sce_milo)), 
                             sample="sample")   
  
traj_design <- data.frame(colData(scrna_sce_milo))[,c("sample", "efficacy")]  
traj_design$sample <- as.factor(traj_design$sample)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample
  
scrna_sce_milo <- calcNhoodDistance(scrna_sce_milo, d = 50)

  
da_results <- testNhoods(scrna_sce_milo, 
                         design = ~ efficacy, 
                         design.df = traj_design)  
scrna_sce_milo <- buildNhoodGraph(scrna_sce_milo)  
da_results <- annotateNhoods(scrna_sce_milo, da_results, coldata_col = "Cell_type_fine_harmony")
da_results$Cell_type_fine_harmony <- factor(da_results$Cell_type_fine_harmony, levels = levels(scrna_subset$Cell_type_fine_harmony))
pp4 <- plotDAbeeswarm(da_results, group.by = "Cell_type_fine_harmony") +
  scale_color_gradient2(low="#50859f", 
                        mid="darkgrey",
                        high="#d66692",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10) +
  theme(axis.text = element_text(colour = 'black')) 
pp5 <- pp4 + #+ scale_y_continuous(limits = c(-6, 6),                      
  #                  breaks = c(-6, -4, -2, 0, 2, 4, 6),                     
  #                  expand = c(0, 0)) +  
  labs(y = "") +  #coord_flip() + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme_classic() +
  theme(legend.position = "none",        
        axis.title.x = element_text(size = 12),        
        axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1),        
        axis.title.y = element_blank(),        
        axis.ticks.y = element_blank(),        
        axis.text.y = element_blank(),        
        axis.line.y = element_blank(),        
        plot.margin = unit(c(0, 1, 0, 0), "cm")) +  
  ylab("Efficacy (Post)")
#ylab(bquote(Log[2] ~ italic('Fold Change')))
pp5

saveRDS(da_results, "da_results_post_myeloid_by_HL.rds")

## Combine plots ----- 
library(patchwork)
p1 <- HL_Group_info + p5 + p7 + p9 + pp5 + 
  plot_layout(nrow = 1, widths = c(0.1, 0.5, 0.5, 0.5, 0.5), guides = "collect") & 
  theme(legend.position = "none", legend.margin = margin(t = -0.5, unit = "cm"),
        legend.text = element_text(#angle = 0,
          face = "italic",
          size = 5#,
          #hjust= 0,
          #vjust = 0
        ))
p1

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Combined_plot_Ti_status_and_Mye_milo_R_4_part_updated_03_07_25", "_v1.pdf"), 
    width = 10, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches
print(p1, newpage = FALSE)
dev.off()

rm(scrna_sce_milo)
rm(da_results)
gc()

### Responders ----- 
Idents(scrna) <- "orig.ident"
scrna_H <- subset(scrna, idents = "H")
Cellratio <- prop.table(table(scrna$Cell_type_fine_harmony, scrna$treatment_status), margin =2) 
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio, Var2~Var1, value.var="Freq") 
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]
cellper <- as.data.frame(t(cellper))
cellper$cell_type <- rownames(cellper)
cellper$cell_type <- factor(cellper$cell_type, levels = c("Mono_FCN1+", "Macro_CCL4+", "Macro_ISGs+", "Macro_MTs+", #"Macro_LYVE1+", 
                                                          "Macro_APOE+", "Macro_Ki67+", "Promyelocytes", "Neu", "cDC1_CLEC9A+", "cDC2_CD1C+", "Mature_DC")
)
mycol <- celltype_colors <- c(
  "#C65D57",     
  "#D4A373",     
  "#7D7461",     
  "#9B6A6C",     
  "#5B5F97",     
  "#DAA49A",     
  "#8D6A4F",     
  "#A4C2A5",     
  "#B07BAC",     
  "#7C9D96",     
  "#D8B4A0",    
  "#645D5C"     
)

p2 <- ggplot(data = cellper,
             aes(x = Pre, 
                 y = Post, 
                 color = cell_type, 
                 label = cell_type
                 #group = , colour = Genus
             )) +
  #geom_smooth(se=FALSE, linetype="dashed", size = 2, color = "#A67EB7") +
  geom_point(aes(size = Post, color = cell_type)) + 
  scale_size_continuous(range = c(3, 15)) + 
  geom_text() + 
  scale_color_manual(values = mycol) + 
  geom_abline(linetype = "dashed") + 
  scale_x_continuous(limits = c(0, 0.5)) + 
  scale_y_continuous(limits = c(0, 0.5)) + 
  theme_bw(base_size = 12) + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA, linewidth = 1, color = "black"), 
        #plot.title = element_text(hjust = 0.5, face = "bold"), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) 
#geom_xspline(spline_shape = -0.4)
#scale_x_discrete(limits=c("RS","RE","VE","SE","LE","P"))+
#scale_colour_manual(values=phy.cols) 
p2

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Cell_prop_geom_point_plot_Myeloid", "_v1.pdf"), 
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(p2, newpage = FALSE)
dev.off()
rm(scrna_H)

### Non-responders -----
Idents(scrna) <- "orig.ident"
scrna_L <- subset(scrna, idents = "L")
Cellratio <- prop.table(table(scrna_L$Cell_type_fine_harmony, scrna_L$treatment_status), margin =2) 
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio, Var2~Var1, value.var="Freq") 
rownames(cellper)<- cellper[,1]
cellper <- cellper[,-1]
cellper <- as.data.frame(t(cellper))
cellper$cell_type <- rownames(cellper)
cellper$cell_type <- factor(cellper$cell_type, levels = c("Mono_FCN1+", "Macro_CCL4+", "Macro_ISGs+", "Macro_MTs+", #"Macro_LYVE1+", 
                                                          "Macro_APOE+", "Macro_Ki67+", "Promyelocytes", "Neu", "cDC1_CLEC9A+", "cDC2_CD1C+", "Mature_DC")
)
mycol <- celltype_colors <- c(
  "#C65D57",     
  "#D4A373",     
  "#7D7461",     
  "#9B6A6C",     
  "#5B5F97",     
  "#DAA49A",     
  "#8D6A4F",     
  "#A4C2A5",     
  "#B07BAC",     
  "#7C9D96",     
  "#D8B4A0",    
  "#645D5C"     
)

p2 <- ggplot(data = cellper,
             aes(x = Pre, 
                 y = Post, 
                 color = cell_type, 
                 label = cell_type
                 #group = , colour = Genus
             )) +
  #geom_smooth(se=FALSE, linetype="dashed", size = 2, color = "#A67EB7") +
  geom_point(aes(size = Post, color = cell_type)) + 
  scale_size_continuous(range = c(3, 15)) + 
  geom_text() + 
  scale_color_manual(values = mycol) + 
  geom_abline(linetype = "dashed") + 
  scale_x_continuous(limits = c(0, 0.5)) + 
  scale_y_continuous(limits = c(0, 0.5)) + 
  theme_bw(base_size = 12) + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = NA, linewidth = 1, color = "black"), 
        #plot.title = element_text(hjust = 0.5, face = "bold"), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) 
#geom_xspline(spline_shape = -0.4)
#scale_x_discrete(limits=c("RS","RE","VE","SE","LE","P"))+
#scale_colour_manual(values=phy.cols) 
p2

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Cell_prop_geom_point_plot_Myeloid_NR", "_v1.pdf"), 
    width = 5, # The width of the plot in inches
    height = 5) # The height of the plot in inches
print(p2, newpage = FALSE)
dev.off()

## Neutrophil +/- Promyelocyte subtyping (Related to Figure S6F) ----- 
Neu <- subset(scrna, idents = "Neu")
### Sigs ----- 
Neu_sigs <- data.frame(readxl::read_xlsx(here("data/Per_cell_type_pipeline_w_new_idents/Myeloid/Neu_DEGs_2024_Cell.xlsx")))
Neu_sigs_list <- list()
for (i in unique(Neu_sigs$Cluster)) {
  Neu_sigs_list[[i]] <- as.character(Neu_sigs$Gene[which(Neu_sigs$Cluster == i)])[1:100]
}

Neu <- NormalizeData(Neu) %>% FindVariableFeatures() %>% ScaleData()
Neu <- AddModuleScore(Neu, features = Neu_sigs_list, name = names(Neu_sigs_list))
names(Neu@meta.data)[c(47:56)] <- names(Neu_sigs_list)

library(dplyr)
genelistnames <- names(Neu_sigs_list)
heatmap_data <- Neu@meta.data %>%
  dplyr::select(group, all_of(genelistnames)) %>%
  group_by(group) %>%
  summarise(across(all_of(genelistnames), mean, na.rm = TRUE)) %>%
  column_to_rownames("group") %>%
  as.matrix() %>%
  t()     
   
heatmap_data_scaled <- t(scale(t(heatmap_data)))

p_value_matrix <- Neu@meta.data %>%
  group_by(group) %>%
  summarise(across(all_of(genelistnames), ~t.test(.x)$p.value)) %>%
  column_to_rownames("group") %>%
  as.matrix() %>%
  t()

library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

cellwidth = 0.7
cellheight = 0.7
cn = dim(as.matrix(heatmap_data_scaled))[2]
rn = dim(as.matrix(heatmap_data_scaled))[1]
w=cellwidth*cn
h=cellheight*rn
Heatmap(as.matrix(heatmap_data_scaled),
        # Height of the blocks
        width = unit(w, "cm"),
        height = unit(h, "cm"),
        #scale = F,
        #rect_gp = gpar(col = "white", lwd = 1.5),
        #border_g = gpar(col = ,lty = 1,lwd = 1.2),
        # Dent formatting 
        column_dend_height = unit(1.5, "cm"), 
        row_dend_width = unit(1.5, "cm"),
        #column_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
        #row_dend_gp = gpar(col = "#9C8D9B",lwd = 1.4),
        # Gap numbers and style 
        #row_split = 6, column_split = c("1", "1", "2", "2"),
        row_gap = unit(2, "mm"),
        column_gap = unit(2, "mm"),
        # Text formats 
        row_title = NULL,column_title = NULL,
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        show_heatmap_legend = T, 
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = 'right',
        show_column_names = T,
        show_row_names = T,
        border = T, 
        col = colorRamp2(breaks = c(-2, 0, 2), colors = c("#50859f","white","#d66692")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if( heatmap_data_scaled[i, j] > 1) {
            grid.text("+++", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] > 0.8 &  heatmap_data_scaled[i, j] <= 1) {
            grid.text("++", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= 0.8 &  heatmap_data_scaled[i, j] > 0.2) {
            grid.text("+", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] <= 0.2 &  heatmap_data_scaled[i, j] > 0) {
            grid.text("+/-", x, y, gp = gpar(fontsize = 8))
          } else if( heatmap_data_scaled[i, j] == 0) {
            grid.text("-", x, y, gp = gpar(fontsize = 8))
          }
        }
)
# Manuval saves 
