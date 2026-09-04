library(here) # project-root-relative paths; run scripts from repository root
##################### Bulk tube morph sig score (Figure S4B) ##################### 
library(msigdbr)

# Load data 
expr_combat <- read.table(here("data/EcoTyper_run_16_06_25/OS_182_post_batch_rem_260731.txt"))
group <- read.csv(here("data/Re-clustering_fine_cell_types_21_05_25/Metadata_260802/updated_bulk_metadta_260802.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8")

hs_df = msigdbr(species = "Homo sapiens") %>% as.data.frame()
hs_C5 = msigdbr(species = "Homo sapiens",
                category = "C5",
                subcategory = "GO:BP") %>% as.data.frame() %>% 
  dplyr::select(gs_cat, gs_subcat, gs_name, gene_symbol)

hs_C5 = hs_C5 %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  dplyr::rename("term" = "gs_name", "gene" = "gene_symbol")
rm(hs_df)

geneset_GOBP <- hs_C5
#geneset <- gson::read.gmt(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/h.all.v2025.1.Hs.symbols.gmt.txt"))
#table(geneset$term)
geneset_GOBP$term <- gsub(pattern = "GOBP_","", geneset_GOBP$term)
geneset_GOBP <- geneset_GOBP[, c(2, 1)]
colnames(geneset_GOBP) <- c("Metagene", "Cell type") 

geneSet_GOBP = split(1:dim(geneset_GOBP)[1], geneset_GOBP$`Cell type`) %>%
  lapply(function(x) {
    geneset_GOBP[unlist(x), 1]
  })
head(geneSet_GOBP, 4) 

geneSet_GOBP_sel <- geneSet_GOBP[which(names(geneSet_GOBP) %in% c(
  "TUBE_MORPHOGENESIS"
  ))]

geneset <- gson::read.gmt(here("data/Re-clustering_fine_cell_types_21_05_25/Osteoblasts_and_cancer/h.all.v2025.1.Hs.symbols.gmt.txt"))
table(geneset$term)
geneset$term <- gsub(pattern = "HALLMARK_","", geneset$term)
geneset <- geneset[, c(2, 1)]
colnames(geneset) <- c("Metagene", "Cell type") 

geneSet = split(1:dim(geneset)[1], geneset$`Cell type`) %>%
  lapply(function(x) {
    geneset[unlist(x), 1]
  })
head(geneSet, 4) 

dat <- as.matrix(expr_combat[, -1]) 
group_sel <- group
geneSet1 <- geneSet[names(geneSet)[c(1, 4, 7, 9, 11, 12, 13, 14, 17, 18, 19, 22, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 
                                     37,  42, 44, 45, 46, 49)]]
GeneSet_comb <- append(geneSet_GOBP_sel, geneSet1)

library(GSVA) 
gsvaPar <- ssgseaParam(exprData = dat, 
                       geneSets = GeneSet_comb,
                       normalize = TRUE) 
ssgsea <- gsva(gsvaPar, verbose = FALSE) 

#### Grouped boxplot of the key signatures ----- 
#dt_for_box <- as.data.frame(t(dt))
dt_for_box <- as.data.frame(t(ssgsea))
dt_for_box$Sample <- rownames(dt_for_box)
dt_for_box <- merge(dt_for_box, group_sel, by = "Sample")
dt_for_box$Strat <- ifelse(dt_for_box$Treatment == "Pre" & dt_for_box$TNR >= 90, "R_Pre", 
                           ifelse(dt_for_box$Treatment == "Pre" & dt_for_box$TNR < 90, "NR_Pre", 
                                  ifelse(dt_for_box$Treatment == "Post" & dt_for_box$TNR >= 90, "R_Post", "NR_Post")))
dt_for_box$Strat <- factor(dt_for_box$Strat, levels = c("R_Pre", "R_Post", "NR_Pre", "NR_Post"))
table(dt_for_box$Strat, dt_for_box$Treatment)

#write.table(dt_for_box, file='OS_Sigs_incl_Tube_260802.txt', quote = F, sep = "\t")

pplist <- list()
dt_for_box <- na.omit(dt_for_box) # Omit missing NAs from TNR 

for (group_ in names(geneSet_GOBP_sel)) {
  cellper_ = dt_for_box %>% dplyr::select(one_of(c('Sample', 'Strat', group_)))
  colnames(cellper_)=c('Sample','Strat','percent')
  cellper_$percent =as.numeric(cellper_$percent)
  cellper_ <- cellper_ %>% group_by(Strat) %>% mutate(upper = quantile(percent,0.75),
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
  
  mycol <- rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D"))
  pp1 = ggplot(cellper_,aes(x=Strat,y=percent)) + 
    geom_jitter(shape =21, aes(fill = Strat), width = 0.25) + 
    geom_boxplot(alpha = 0.5, size = 1, outlier.shape = NA, aes(fill = Strat)) + 
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
  compare_means(percent ~ Strat,  data = cellper_)
  #my_comparisons <-list(c("H_Pre","H_Post") 
  #c("H_Pre","L_Pre"), 
  #c("L_Pre","L_Post")
  #)
  
  pp1 = pp1 + geom_signif(comparisons = list(c("R_Pre","R_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position, 
                          #paired = T, 
                          #test = "t.test",
                          test.args = list(#paired = TRUE, 
                            alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  #my_comparisons1 <-list(c("L_Pre","L_Post"))
  pp1 = pp1 + geom_signif(comparisons = list(c("NR_Pre","NR_Post")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))}, 
                          y_position = p_position, 
                          #paired = T, 
                          #test = "t.test",
                          test.args = list(#paired = TRUE, 
                            alternative = "two.sided"),
                          textsize = 2, tip_length = 0,
                          parse = TRUE, 
                          color = "black")  
  pp1 = pp1 + geom_signif(comparisons = list(c("R_Pre","NR_Pre")), 
                          map_signif_level = function(p) {paste("italic(P) == ", sprintf("%.2g", p))},
                          y_position = p_position1, 
                          #paired = T, 
                          #test = "t.test", 
                          test.args   = list(#paired = TRUE, 
                            alternative = "two.sided"),
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

pps <- plot_grid(pplist[[1]]
)
pps

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "Tube_morph_sig", "_change_H_Pre_H_Post_L_Pre_L_Post.pdf"), 
    width = 3, # The width of the plot in inches
    height = 3) # The height of the plot in inches
print(pps, newpage = FALSE)
dev.off()

##### Violin plots for capEC genes (Figure S4C) ----- 
Endo <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Endo/EC_w_Cap_updated_251231.rds"))

Endo[["RNA"]] <- as(object = Endo[["RNA"]], Class = "Assay")
my_comparisons <- list(c("H_Pre", "H_Post"), c("H_Pre", "L_Pre"), c("L_Pre", "L_Post"), c("H_Post", "L_Post"))
levels(Endo) 

### Stalk
Endo_subset <- subset(Endo, idents = "Stalk-like")
pplist <- list()
for (i in c("ICAM1", "DNASE1L3", "SELE", "SELP")) {
  pplist[[i]] <- FeatureStatPlot(srt = Endo_subset, group.by = "group_anno",  
                                 stat.by = i, 
                                 add_box = TRUE,  
                                 palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")), 
                                 bg_palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")),
                                 #comparisons = my_comparisons,  
                                 add_trend = T, 
                                 bg.by = 'group_anno') + theme(legend.position = "none") + 
    theme(legend.position = "none", 
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.title.y = element_blank())
}

pps <- plot_grid(pplist[[1]], 
                 pplist[[2]],
                 pplist[[3]], 
                 pplist[[4]],
                 #pplist[[5]],
                 #pplist[[6]],
                 #pplist[[7]],
                 #pplist[[8]],
                 #pplist[[9]], 
                 #pplist[[10]],
                 #pplist[[11]], 
                 #pplist[[12]],
                 #pplist[[13]],
                 ncol = 4
)
pps

pdf("Stalk-like_gene_exp_260103.pdf", height = 3, width = 8)
print(pps)
dev.off()

## Trans FLThi
levels(Endo)
Endo_subset <- subset(Endo, idents = "Transition_FLT1hi")
pplist <- list()
for (i in c("FLT1", "KDR")) {
  pplist[[i]] <- FeatureStatPlot(srt = Endo_subset, group.by = "group_anno",  
                                 stat.by = i, 
                                 add_box = TRUE,  
                                 palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")), 
                                 bg_palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")),
                                 #comparisons = my_comparisons,  
                                 add_trend = T, 
                                 bg.by = 'group_anno') + theme(legend.position = "none") + 
    theme(legend.position = "none", 
          axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.title.y = element_blank())
}

pps <- plot_grid(
  pplist[[1]], 
  pplist[[2]],
  #pplist[[5]],
  #pplist[[6]],
  #pplist[[7]],
  #pplist[[8]],
  #pplist[[9]], 
  #pplist[[10]],
  #pplist[[11]], 
  #pplist[[12]],
  #pplist[[13]],
  ncol = 2
)
pps

pdf("Transition_FLT1hi_gene_exp_260103.pdf", height = 3, width = 8)
print(pps)
dev.off()

## Tip-like 
levels(Endo)
Endo_subset <- subset(Endo, idents = "Tip-like")
pplist <- list()
for (i in c("PGF", "PXDN", "CD99", "PDGFB", "APLN", "CD276", "VSIR")) {
  pplist[[i]] <- FeatureStatPlot(srt = Endo_subset, group.by = "group_anno",  
                                 stat.by = i, 
                                 add_box = TRUE,  
                                 palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")), 
                                 bg_palcolor = rev(c("#6B798E", "#0f5688", "#D9A0B3","#54426D")),
                                 #comparisons = my_comparisons,  
                                 add_trend = T, 
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
                 #pplist[[8]],
                 #pplist[[9]], 
                 #pplist[[10]],
                 #pplist[[11]], 
                 #pplist[[12]],
                 #pplist[[13]],
                 ncol = 7
)
pps

pdf("Tip_gene_exp_260103.pdf", height = 3, width = 14)
print(pps)
dev.off()

## GO enrichment of capEC subset markers (Figure S4D) ----- 
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

mar <- FindAllMarkers(scrna_cap, group.by = "Cell_type_fine_harmony", only.pos = T)

# Tip
gene <- mar$gene[mar$cluster == "Tip-like"][1:100]
enricher_res_Tip <- enricher(gene, TERM2GENE = hs_C5, pvalueCutoff = 1,
                             pAdjustMethod = "BH")
enricher_res_Tip <- enricher_res_Tip@result
enricher_res_Tip$pval_log <- -log10(enricher_res_Tip$pvalue)
#gsea_results <- gsea_results[order(gsea_results$enrichmentScore,decreasing=T),]

go_select <- enricher_res_Tip %>% 
  filter(pvalue < 0.05) %>% 
  arrange(pvalue) %>% 
  #head(50)  %>% 
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

go_select <- go_select[which(go_select$ID %in% 
                               c("GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX", 
                                 "GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT", 
                                 "HP_ABNORMALITY_OF_BLOOD_CIRCULATION", 
                                 "GOBP_LEUKOCYTE_MIGRATION",
                                 "GOBP_REGULATION_OF_T_CELL_ACTIVATION"
                               )
), ]
go_select$pval_log <- -log10(go_select$pvalue)
go_select$Cluster <- "Tip-like"
#scale_factor <- max(go_select$pval_log) / max(go_select$RichFactor)
Mat_for_plot <- go_select

# FLT1hi Tran 
gene <- mar$gene[mar$cluster == "Transition_FLT1hi"][1:100]
enricher_res_FLThi_Trans <- enricher(gene, TERM2GENE = hs_C5, pvalueCutoff = 1,
                                     pAdjustMethod = "BH")
enricher_res_FLThi_Trans <- enricher_res_FLThi_Trans@result
enricher_res_FLThi_Trans$pval_log <- -log10(enricher_res_FLThi_Trans$pvalue)
go_select <- enricher_res_FLThi_Trans %>% 
  filter(pvalue < 0.05) %>% 
  arrange(pvalue) %>% 
  #head(50)  %>% 
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

go_select <- go_select[which(go_select$ID %in% 
                               c("GOBP_CELL_JUNCTION_ASSEMBLY", 
                                 "GOBP_PLATELET_ACTIVATION", 
                                 "GOBP_REGULATION_OF_CELL_GROWTH", 
                                 "GOBP_MAINTENANCE_OF_LOCATION",
                                 "GOBP_NEGATIVE_REGULATION_OF_TRANSPORT"
                               )
), ]
go_select$pval_log <- -log10(go_select$pvalue)
go_select$Cluster <- "Transition_FLT1hi"

Mat_for_plot <- rbind(Mat_for_plot, go_select)

# FLT1lo Tran 
gene <- mar$gene[mar$cluster == "Transition_FLT1lo"][1:100] 
enricher_res_FLTlo_Trans <- enricher(gene, TERM2GENE = hs_C5, pvalueCutoff = 1,
                                     pAdjustMethod = "BH")
enricher_res_FLTlo_Trans <- enricher_res_FLTlo_Trans@result
enricher_res_FLTlo_Trans$pval_log <- -log10(enricher_res_FLTlo_Trans$pvalue)
go_select <- enricher_res_FLTlo_Trans %>% 
  filter(pvalue < 0.05) %>% 
  arrange(pvalue) %>% 
  #head(50)  %>% 
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

go_select <- go_select[which(go_select$ID %in% 
                               c("GOCC_MITOCHONDRIAL_PROTEIN_CONTAINING_COMPLEX", 
                                 "GOBP_ALTERNATIVE_MRNA_SPLICING_VIA_SPLICEOSOME", 
                                 "GOBP_PROTEIN_DNA_COMPLEX_ASSEMBLY", 
                                 "GOBP_MRNA_TRANSPORT",
                                 "GOMF_HELICASE_ACTIVITY"
                               )
), ]
go_select$pval_log <- -log10(go_select$pvalue)
go_select$Cluster <- "Transition_FLT1lo"

Mat_for_plot <- rbind(Mat_for_plot, go_select)

# Stalk-like 
gene <- mar$gene[mar$cluster == "Stalk-like"][1:100] 
enricher_res_FLTlo_Trans <- enricher(gene, TERM2GENE = hs_C5, pvalueCutoff = 1,
                                     pAdjustMethod = "BH")
enricher_res_FLTlo_Trans <- enricher_res_FLTlo_Trans@result
enricher_res_FLTlo_Trans$pval_log <- -log10(enricher_res_FLTlo_Trans$pvalue)
go_select <- enricher_res_FLTlo_Trans %>% 
  filter(pvalue < 0.05) %>% 
  arrange(pvalue) %>% 
  #head(50)  %>% 
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

go_select <- go_select[which(go_select$ID %in% 
                               c("GOBP_LEUKOCYTE_PROLIFERATION", 
                                 "GOBP_LEUKOCYTE_CELL_CELL_ADHESION", 
                                 "GOBP_RESPONSE_TO_TYPE_II_INTERFERON", 
                                 "GOBP_RESPONSE_TO_TUMOR_CELL",
                                 "GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE"
                               )
), ]
go_select$pval_log <- -log10(go_select$pvalue)
go_select$Cluster <- "Stalk-like"

Mat_for_plot <- rbind(Mat_for_plot, go_select)
Mat_for_plot$Description <- factor(Mat_for_plot$Description, levels = rev(Mat_for_plot$Description))

mycol <- rev(c("#B08D57", "#6B4E3D", "#9A9A9A","#7A1E2D"))
p <- ggplot() +
  geom_bar(data = Mat_for_plot,
           aes(x = -log10(pvalue), y = Description, fill = Cluster),
           width = 0.5, 
           stat = 'identity') +
  theme_classic() + 
  scale_x_continuous(expand = c(0,0)) +
  theme(axis.text.y = element_blank()) +
  geom_text(data = Mat_for_plot,
            aes(x = 0.1, 
                y = Description, 
                label = Description),
            size = 4.5,
            hjust = 0) +
  geom_text(data = Mat_for_plot,
            aes(x = 0.1, y = Description, label = geneID, color = Cluster),
            size = 4,
            fontface = 'italic',
            hjust = 0,
            vjust = 2.3) +
  labs(x = '-Log10P', 
       title = 'Enriched pathways of top 100 marker genes') + 
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.ticks.y = element_blank()) + 
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol)
p

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "CapEC_subset_GOBP", ".pdf"), 
    width = 6, 
    height = 12) 
print(p, newpage = FALSE)
dev.off()

## ST demonstrationo of tip and stalk signature distrbution (Figure S4E) ------ 
Tip_sig <- list(c("ADM", "ANGPT2", "ANKRD37", "APLN", "C1QTNF6", "CD93", "CLDN5", "COL4A1", "COL4A2", 
                  "COTL1", "CXCR4", "DLL4", "EDNRB", "ESM1", "FSCN1", "GPIHBP1", "HSPG2", "IGFBP3", 
                  "INHBB", "ITGA5", "JUP", "KCNE3", "KCNJ8", "KDR", "LAMA4", "LAMB1", "LAMC1", "LXN", 
                  "MARCKS", "MARCKSL1", "MCAM", "MEST", "MYH9", "MYO1B", "N4BP3", "NID2", "NOTCH4", 
                  "PDGFB", "PGF", "PLOD1", "PLXND1", "PMEPA1", "PTN", "RAMP3", "RBP1", "RGCC", "RHOC", 
                  "SMAD1", "SOX17", "SOX4", "SPARC", "TCF4", "UNC5B", "VIM"))
#scrna <- AddModuleScore(scrna, features = Tip_sig, name = "Tip", assay = "RNA") 
#names(scrna@meta.data)[47] <- "Tip"
Stalk_sig <- list(c("ACKR1", "AQP1", "C1QTNF9", "CD36", "CSRP2", "EHD4", "FBLN5", "HSPB1", "LIGP1", 
                    "IL6ST", "JAM2", "LGALS3", "LRG1", "MEOX2", "PLSCR2", "CAVIN2", "SELP", "SPINT2", 
                    "TGFBI", "TGM2", "TMEM176A", "TMEM176B", "TMEM252", "TSPAN7", "FLT1", "VWF"))
#scrna <- AddModuleScore(scrna, features = Stalk_sig, name = "Stalk", assay = "RNA") 
#names(scrna@meta.data)[48] <- "Stalk"
#scrna@meta.data[38] <- NULL

SAHST005 <- readRDS(here("data/New_ST_251011/Output/SAHST005_cloned.rds"))
names(SAHST005@meta.data)
SAHST005 <- AddModuleScore(SAHST005, features = Tip_sig, name = "Tip", assay = "SCT") 
SAHST005 <- AddModuleScore(SAHST005, features = Stalk_sig, name = "Stalk", assay = "SCT") 
names(SAHST005@meta.data)[81] <- "Tip"
names(SAHST005@meta.data)[82] <- "Stalk"

SpatialFeaturePlot(SAHST005, c("Tip", "Stalk"))

#### Custom Plot 
# Tip
SAHST005@reductions$spatial = SAHST005@reductions$umap
SAHST005@reductions$spatial@key = 'spatial_'
SAHST005@reductions$spatial@cell.embeddings = as.matrix(SAHST005@images$image@coordinates[,c(3,2)])
SAHST005@reductions$spatial@cell.embeddings[,2] = -SAHST005@reductions$spatial@cell.embeddings[,2]
colnames(SAHST005@reductions$spatial@cell.embeddings) = c('spatial_1','spatial_2')
data_mat <- data.frame(SAHST005$Tip, #SAHST005$S_phase, SAHST002_plot$Os_min, 
                       SAHST005@reductions$spatial@cell.embeddings)

# Clone DimPlot 
# color_scale_val <- range(ST_merge$ECM)
color_scale_val <- range(SAHST005$Tip)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST005.Tip)) +
  geom_point(shape = 16, size = 1) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026"), 
                        #colours = c("lightblue","lightyellow","red"),  
                        limits = color_scale_val) + 
  #scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
  #                                  "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
  #                                  "#d73027", "#a50026"), 
  #colours = c("lightblue","lightyellow","red"),  
  #                      limits = color_scale_val) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 1.7) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf('SAHST005_Tip_Spatial_mapping.pdf', height = 3, width = 3)
print(p1)
dev.off()

# Stalk 
data_mat <- data.frame(SAHST005$Stalk, #SAHST005$S_phase, SAHST002_plot$Os_min, 
                       SAHST005@reductions$spatial@cell.embeddings)
color_scale_val <- range(SAHST005$Stalk)
p1 <- ggplot(data_mat, aes(x = spatial_1, y = spatial_2, colour = SAHST005.Stalk)) +
  geom_point(shape = 16, size = 1) +
  scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
                                    "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
                                    "#d73027", "#a50026"), 
                        #colours = c("lightblue","lightyellow","red"),  
                        limits = color_scale_val) + 
  #scale_color_gradientn(colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
  #                                  "#e0f3f8", "#ffffbf", "#fee090", "#f46d43",
  #                                  "#d73027", "#a50026"), 
  #colours = c("lightblue","lightyellow","red"),  
  #                      limits = color_scale_val) + 
  #scale_color_viridis(option = "H", direction = 1, limits = color_scale_val) + 
  coord_fixed(ratio = 1.7) +
  theme(panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.text = element_blank(),
        axis.ticks=element_blank())
p1

pdf('SAHST005_Stalk_Spatial_mapping.pdf', height = 3, width = 3)
print(p1)
dev.off()

## TNK and Endo subset co-occurence (Fig 4J) ------ 
library(Seurat) 
library(reshape2) 
library(Hmisc) 
library(reshape2) 
library(ggplot2) 
library(tidyverse) 

scRNA <- readRDS(here("data/CellPhoneDB_pipeline_03_05_25/temp_all_cells_prior_to_reint.rds"))
scrna <- readRDS(here("data/Re-clustering_fine_cell_types_21_05_25/Endo/EC_w_Cap_updated_251231.rds"))

scRNA$Cell_type_fine_harmony_keep <- scRNA$Cell_type_fine_harmony
scRNA$Cell_type_fine_harmony <- scRNA$cluster

levels(scrna) 
scrna$Cell_type_fine_harmony <- factor(scrna$Cell_type_fine_harmony, 
                                       levels = c("artEC", "Tip-like", "Transition_FLT1hi", "Transition_FLT1lo", 
                                                  "Stalk-like", "venEC", "LEC",  "EndoMT-I", "EndoMT-II", "Cycling_ECs"))

Cellratio_1 <- prop.table(table(scRNA$Cell_type_fine_harmony, scRNA$sample), margin = 2) %>% data.frame()
cellper_i <- reshape2::dcast(Cellratio_1, Var2~Var1, value.var="Freq")  
rownames(cellper_i)<- cellper_i[,1] 
cellper_i <- cellper_i[,-1] 

Cellratio_2 <- prop.table(table(scrna$Cell_type_fine_harmony, scrna$sample), margin = 2)  %>% data.frame()
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
plot_df$label <- paste0(sprintf("%.2f", plot_df$value.x) 
)
plot_df$label <- ifelse(plot_df$value.x > 0.3 | plot_df$value.x < -0.3, plot_df$label, "")

#plot_df <- plot_df[plot_df$Var1 %in% unique(scRNA$Cell_type_fine_harmony), ]
plot_df <- plot_df[plot_df$Var1 %in% "TandNK", ]
plot_df <- plot_df[plot_df$Var2 %in% unique(scrna$Cell_type_fine_harmony), ]

p <- ggplot() +
  geom_tile(data = plot_df, aes(x = Var2, y = Var1, fill = value.x), color = NA) +
  geom_text(data = plot_df, aes(x = Var2, y = Var1, label = label),
            colour = ifelse(plot_df$value.y < 0.05, "red", "black"), 
            size = 4) +
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
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )

time <- gsub(" ", "_", Sys.time())
time <- gsub("-", "_", time)
time <- gsub(":", "_", time)
pdf(file = paste0(time, "EC_subsets_all_major_celltype", "_co-occur_plot.pdf"), 
    width = 10, 
    height = 6) 
print(p, newpage = FALSE)
dev.off()
