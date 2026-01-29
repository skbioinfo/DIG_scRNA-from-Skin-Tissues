.libPaths("/home/jinylim/R/x86_64-pc-linux-gnu-library/4.2")
homelib = "/home/jinylim/R/x86_64-pc-linux-gnu-library/4.2"
library(ggplot2, lib.loc = homelib)
# library(Seurat, lib.loc = homelib)
library(dplyr)
library(cowplot)
library(ggrepel)

source("/clinicalInfo.R")
meta <- read.csv("/Seuratmeta.csv")
meta_ct2 <-read.csv( "~/ingleR.csv")

all(meta$uid == meta_ct2$X)
meta$celltype1_old <- meta$celltype1
meta$celltype1 <- meta_ct2$SingleR_label
meta$max.score <- meta_ct2$max.score
meta$delta.next <- meta_ct2$delta.next

meta$UMAP1 <- meta$umap_1
meta$UMAP2 <- meta$umap_2

#---celltype1------
meta$celltype1 <- dplyr::recode(
  meta[['celltype2']],
  ## DC
  "DC1" = "DC",
  "DC2" = "DC",
  "MigDC" = "DC",
  "mDC" = "DC",
  "pDC" = "DC",
  "moDC_1" = "DC",
  "moDC_3" = "DC",
  
  ## Fibroblast
  "F1" = "Fibroblast",
  "F2" = "Fibroblast",
  "F3" = "Fibroblast",
  
  ## T cells
  "Tc" = "T cell",
  "Th" = "T cell",
  "Treg" = "T cell",
  "Tc17_Th17" = "T cell",
  
  ## KC
  "Differentiated_KC"  = "Keratinocyte",
  "Differentiated_KC*" = "Keratinocyte",
  "Undifferentiated_KC" = "Keratinocyte",
  "Proliferating_KC" = "Keratinocyte",
  
  ## Macrophage / Mono
  "Macro_1" = "Macrophage",
  "Macro_2" = "Macrophage",
  "Mono_mac" = "Macrophage",
  "Inf_mac"  = "Macrophage",
  
  ## Langerhans
  "LC_1" = "Langerhans cell",
  "LC_2" = "Langerhans cell",
  "LC_3" = "Langerhans cell",
  "LC_4" = "Langerhans cell",
  
  ## Endothelial-like
  "LE1" = "Endothelial cell",
  "LE2" = "Endothelial cell",
  "VE1" = "Endothelial cell",
  "VE" = "Endothelial cell",
  "VE2" = "Endothelial cell",
  "VE3" = "Endothelial cell",
  
  ## Pericyte
  "Pericyte_1" = "Pericyte",
  "Pericyte_2" = "Pericyte",
  "Schwann_1" = "Schwann",
  
  ## ILC
  "ILC1_3" = "ILC",
  "ILC1_NK" = "ILC",
  "ILC2" = "ILC",
  
  ## 그대로 유지
  .default = meta[['celltype2']]
)


meta$celltype1 <- dplyr::recode(
  meta[['celltype2']],
  ## DC
  "DC1" = "DC",
  "DC2" = "DC",
  "MigDC" = "DC",
  "mDC" = "DC",
  "pDC" = "DC",
  "moDC_1" = "DC",
  "moDC_3" = "DC",
  
  ## Fibroblast
  "F1" = "Fibroblast",
  "F2" = "Fibroblast",
  "F3" = "Fibroblast",
  
  ## T cells
  "Tc" = "T cell",
  "Th" = "T cell",
  "Treg" = "T cell",
  "Tc17_Th17" = "T cell",
  
  ## KC
  "Differentiated_KC"  = "Keratinocyte",
  "Differentiated_KC*" = "Keratinocyte",
  "Undifferentiated_KC" = "Keratinocyte",
  "Proliferating_KC" = "Keratinocyte",
  
  ## Macrophage / Mono
  "Macro_1" = "Macrophage",
  "Macro_2" = "Macrophage",
  "Mono_mac" = "Macrophage",
  "Inf_mac"  = "Macrophage",
  
  ## Langerhans
  "LC_1" = "Langerhans cell",
  "LC_2" = "Langerhans cell",
  "LC_3" = "Langerhans cell",
  "LC_4" = "Langerhans cell",
  
  ## Endothelial-like
  "LE1" = "Endothelial cell",
  "LE2" = "Endothelial cell",
  "VE1" = "Endothelial cell",
  "VE" = "Endothelial cell",
  "VE2" = "Endothelial cell",
  "VE3" = "Endothelial cell",
  
  ## Pericyte
  "Pericyte_1" = "Pericyte",
  "Pericyte_2" = "Pericyte",
  "Schwann_1" = "Schwann",
  
  ## ILC
  "ILC1_3" = "ILC",
  "ILC1_NK" = "ILC",
  "ILC2" = "ILC",
  
  ## 그대로 유지
  .default = meta[['celltype2']]
)
## R — immune / non-immune 순서로 factor 정렬

celltype1_order <- c(
  ## immune
  "T cell",
  "NK",
  "Plasma",
  "ILC",
  "Macrophage",
  "DC",
  "Langerhans cell",
  "Mast_cell",
  
  ## non-immune
  "Keratinocyte",
  "Fibroblast",
  "Endothelial cell",
  "Pericyte",
  "Schwann",
  "Melanocyte"
)

# 예: pred$labels_coarse 또는 meta$celltype
meta$celltype1 <- factor(meta$celltype1, levels = celltype1_order)
preferred_within <- list(
  "DC" = c("DC1","DC2","MigDC","moDC_1","moDC_3","mDC","pDC"),
  "T cell" = c("Tc", "Th", "Tc17_Th17", "Treg"),
  "ILC" = c("ILC1_3","ILC1_NK","ILC2"),
  "Macrophage" = c("Macro_1","Macro_2","Mono_mac","Inf_mac"),
  "Langerhans cell" = c("LC_1","LC_2","LC_3","LC_4"),
  "Keratinocyte" = c("Differentiated_KC","Differentiated_KC*","Undifferentiated_KC","Proliferating_KC"),
  "Fibroblast" = c("F1","F2","F3"),
  "Endothelial cell" = c("VE1","VE2","VE3","LE1","LE2"),
  "Pericyte" = c("Pericyte_1","Pericyte_2"),
  "Schwann" = c("Schwann_1")
)

celltype2_order <- unlist(lapply(levels(meta$celltype1), function(ct1){
  ct2s <- unique(as.character(meta$celltype2)[meta$celltype1 == ct1])
  pref <- preferred_within[[ct1]]
  if (!is.null(pref)) {
    c(intersect(pref, ct2s), sort(setdiff(ct2s, pref)))
  } else {
    sort(ct2s)
  }
}))

meta$celltype2 <- factor(meta$celltype2, levels = celltype2_order)

names(col_celltype1) <- celltype1_order
col_celltype2 <- c(col_celltype2, col_celltype1)
names(col_celltype2) <- celltype2_order

col_celltype1 <- c(
  ## Immune (blue–purple–red 계열)
  "T cell"            = "#4C72B0",  # blue
  "NK"                = "#6BAED6",  # light blue
  "Plasma"            = "#9ECAE1",
  "ILC"               = "#E8C4F2",  # very light blue
  
  "Macrophage"        = "#C44E52",  # red
  "DC"                = "#DD8452",  # orange-red
  "Langerhans cell"   = "#E24A33",  # strong red-orange
  "Mast_cell"         = "#8C1D18",  # dark red
  
  ## Non-immune (green–orange–brown 계열)
  "Keratinocyte"      = "#F28E2B",  # orange
  "Fibroblast"        = "#FFBE7D",  # light orange
  
  "Endothelial cell"  = "#59A14F",  # green
  "Pericyte"          = "#8CD17D",  # light green
  
  "Schwann"           = "#B6992D",  # olive
  "Melanocyte"        = "#1B9E77"   # dark teal
)

## R — celltype2 색상 (celltype1 계열 유지)

col_celltype2 <- c(
  ## T cell (blue)
  "Treg"        = "#3B5B92",
  "Th"          = "#4C72B0",
  "Tc"          = "#6C8ED5",
  "Tc17_Th17"   = "#9BB7E5",
  
  ## NK / ILC (light blue–cyan)
  "NK"          = "#6BAED6",
  "ILC1_3"      = "#9ECAE1",
  "ILC1_NK"     = "#BFDCEC",
  "ILC2"        = "#D6EAF8",
  
  ## DC (orange-red)
  "DC1"         = "#D3722C",
  "DC2"         = "#DD8452",
  "MigDC"       = "#E6A15C",
  "moDC_1"      = "#F1B87A",
  "moDC_3"      = "#F6C89F",
  
  ## Macrophage (red)
  "Macro_1"     = "#B03A2E",
  "Macro_2"     = "#C44E52",
  "Mono_mac"    = "#D98880",
  "Inf_mac"     = "#E6B0AA",
  
  ## Langerhans (dark red)
  "LC_1"        = "#8C1D18",
  "LC_2"        = "#A93226",
  "LC_3"        = "#C0392B",
  "LC_4"        = "#E74C3C",
  
  ## Mast
  "Mast_cell"   = "#641E16",
  
  ## Keratinocyte (orange)
  "Undifferentiated_KC" = "#F5B041",
  "Proliferating_KC"   = "#F8C471",
  "Differentiated_KC"  = "#F28E2B",
  "Differentiated_KC*" = "#FAD7A0",
  
  ## Fibroblast (light orange/pink)
  "F1"          = "#F5B7B1",
  "F2"          = "#FF9896",
  "F3"          = "#FADBD8",
  
  ## Endothelial (green)
  "VE1"         = "#2E8B57",
  "VE2"         = "#59A14F",
  "VE3"         = "#86BC86",
  "LE1"         = "#9FD9A3",
  "LE2"         = "#C7E9C0",
  
  ## Pericyte (light green)
  "Pericyte_1"  = "#7FBF7B",
  "Pericyte_2"  = "#A6DBA0",
  
  ## Schwann / Melanocyte
  "Schwann_1"   = "#B6992D",
  "Melanocyte" = "#1B9E77"
)

## factor level 순서에 맞게 정렬 (권장)
col_celltype2 <- col_celltype2[levels(meta$celltype2)]

## 확인 
#--------------
table(pred$labels_coarse)

meta$celltype2 <- meta_ct2$celltype2_2

dir <- '/mnt/ultrastarTrois1/2025_Spatial_JY/AnalFigure/Basic3_singleR-aggre'
meta$condition <- factor(meta$condition, levels = c("Healthy", "AD", "PSO", "HS"))
meta$celltype1 <- factor(meta$celltype1, levels = names(col_celltype1))
meta$celltype2 <- factor(meta$celltype2, levels = names(col_celltype2))

celltype1_loc <- meta %>%
  group_by(celltype1) %>%
  summarise(
    UMAP1 = mean(UMAP1, na.rm = TRUE),
    UMAP2 = mean(UMAP2, na.rm = TRUE)
  )

celltype2_loc <- meta %>%
  group_by(celltype2) %>%
  summarise(
    UMAP1 = mean(UMAP1, na.rm = TRUE),
    UMAP2 = mean(UMAP2, na.rm = TRUE)
  )

SingleR_celltype2_loc <- meta %>%
  group_by(SingleR_celltype2) %>%
  summarise(
    UMAP1 = mean(UMAP1, na.rm = TRUE),
    UMAP2 = mean(UMAP2, na.rm = TRUE)
  )

#--- umap-----------
all(meta$X == meta_singleR$X)
meta$SingleR_celltype2 <- meta_singleR$SingleR_celltype2


#----celltype-----------
meta -> meta_old
meta <- meta_old %>% dplyr::filter(max.score  > 0.25)
p_umap1 <- ggplot(meta, aes(UMAP1, UMAP2, color = celltype1_old)) +
  geom_point(size = 0.25) +
  # geom_text_repel(
  #   data = celltype2_loc,
  #   aes(label = celltype2),
  #   color = "black",
  #   size = 3,
  #   fontface = "bold",
  #   box.padding = 0.3,
  #   point.padding = 0.2,
  #   max.overlaps = Inf
  # ) +
  theme_bw()+
  theme(
    legend.background = element_rect(fill = "white", color = NA)) +
  scale_color_manual(values = col_celltype1) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(plot = p_umap1 + facet_wrap(~condition), 
       file.path(dir, "Umap_celltype1_filtered0.25_condition.png"),width = 7, height = 5)


p_umap1 <- ggplot(meta, aes(UMAP1, UMAP2, color = max.score )) +
  geom_point(size = 0.25) +
  theme_bw()+
  theme(
    legend.background = element_rect(fill = "white", color = NA)) +
  theme(legend.position = "right") 
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(plot = p_umap1, 
       file.path(dir, "Umap_prediction.score.max.png"),width = 8, height = 5)

# meta_singleR <- read.csv("~/2025_Spatial_JY/data/xenium_meta_singleR_1.csv") 


#----Umap by condition--------
pdf(file.path(dir, "Umap_celltype1_sample.pdf"),width = 5, height = 5)
conds <- sort(unique(meta$condition))
for (cond in conds) {
  df <- meta %>% filter(condition == cond)
  p_umap1 <- ggplot(df, aes(UMAP1, UMAP2, color = celltype1)) +
    geom_point(size = 0.25) +
    theme_bw() +
    scale_color_manual(values = col_celltype1) +
    theme(legend.position = "none") +
    ggtitle(cond)
  print(p_umap1)  
}
dev.off()
#--------


#--- Barplot by condition---------
ggplot(meta, aes(x = slide_ID, fill = celltype1))+ 
  geom_bar(position = 'fill')+
  scale_fill_manual(values =col_celltype1) +
  facet_grid(~condition)+ theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())
ggsave( file.path(dir, "BarFilled_celltype1.pdf"),width = 5, height = 3.5)

ggplot(meta %>% dplyr::filter(condition %in% c("Healthy", "PSO")),
       aes(x = slide_ID, fill = celltype2))+ 
  geom_bar(position = 'fill')+
  scale_fill_manual(values =col_celltype2) +
  facet_grid(~condition)+ theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())
ggsave( file.path(dir, "BarFilled_celltype2_PSO.pdf"),width = 6, height = 3.5)

ggplot(meta, aes(x = slide_ID, fill = celltype2))+ 
  geom_bar(position = 'fill')+
  scale_fill_manual(values =col_celltype2) +
  facet_grid(~condition)+ theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_blank())
ggsave( file.path(dir, "BarFilled_celltype1.pdf"),width = 5, height = 3.5)

#----spatial by condition--------
pdf(file.path(dir, "Spatialmap_celltype12_sample.pdf"),width = 4, height = 5)
sampleid <- sort(unique(meta$sample_id))
for (sample_i in sampleid) {
  df <- meta %>% filter(sample_id == sample_i)
  p_stmap1 <- ggplot(df, aes(spatial_x, spatial_y, color = celltype1)) +
    geom_point(size = 0.25) +
    theme_bw() +
    scale_color_manual(values = col_celltype1) +
    theme(legend.position = "none") +
    ggtitle(sample_i)
  print(p_stmap1)  
  
  p_stmap <- ggplot(df, aes(spatial_x, spatial_y, color = celltype2)) +
    geom_point(size = 0.25) +
    theme_bw() +
    scale_color_manual(values = col_celltype2) +
    theme(legend.position = "none") +
    ggtitle(sample_i)
  print(p_stmap)  
}
dev.off()
#--------
genes_dotplot <- c(
  "MPZ","SOX10","PLP1","PMP22","S100B", #mela
  "ACTA2","RGS5","MCAM","CSPG4","NOTCH3", #per
  "PECAM1","VWF","PLVAP","EMCN","KDR", #EN
  "LYVE1","PROX1","PDPN","CCL21","FLT4", #En
  "COL1A1","VCAN","PDGFRA","CXCL12","THBS2", #Fib
  "KRT17","KRT6B","DMKN","DSG1","COL17A1", #Kera
  "TPSB2","TPSAB1","CTSG","MS4A2","KIT", #MAs
  "CXCL8","CD83","FCER1A","CLEC10A","CD1C", #Mac
  "CAPN3","TYR","DCT","MLANA","PMEL", 
  "KLRD1","NKG7","GNLY","PRF1",   # ILC/NK
  "PTPRC","CD3E","IL7R","CCL5","LTB" #T
)
genes_use <- intersect(genes_dotplot, rownames(qry))
qry <-  AddMetaData(qry, metadata = meta %>% dplyr::select("celltype1",
                                                            "celltype2"))
qry <-qry[,colnames(qry)[!is.na(qry$celltype1)]]

dp <- DotPlot(qry, features = genes_use, group.by = 'celltype1')
dp_df <- dp$data
dp_df$id <- factor(dp_df$id)
dp_df$id <- factor(dp_df$id, levels = rev(levels(dp_df$id )))

# 3) ggplot으로 직접 dotplot 다시 그리기
p <- ggplot(dp_df, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  ) +
  scale_color_distiller(palette = "RdBu", direction = -1) +
  labs(x = NULL, y = NULL, color = "Avg. expression (scaled)", size = "% cells")
p

