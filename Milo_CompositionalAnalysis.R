## R (miloR) — xenium.sce에서 Healthy / PSO만 subset 후 DA
library(SingleCellExperiment)
library(miloR)
library(scater)
library(dplyr)
source("/clinicalInfo.R")
## 0) Assign the columns #338910
condition_col <- "condition"   # e.g., "condition"
sample_col    <- "sample_id"   # e.g., "sample_id" or "Sample"

## 1) subset (Healthy, PSO) #64555

library(anndataR)
singleR <- read_h5ad("/singleR.h5ad")
singleR.sce <- singleR$as_SingleCellExperiment()
common_cells <- intersect(colnames(xenium.sce), colnames(singleR.sce))
idx_x <- match(common_cells, colnames(xenium.sce))
idx_s <- match(common_cells, colnames(singleR.sce))

# singleR label로 덮어쓰기
xenium.sce$celltype2_2 <- as.character(xenium.sce$celltype2)
singleR.sce$labels_merged<- as.character(singleR.sce$labels_merged)
xenium.sce$celltype2_2[idx_x] <- singleR.sce$labels_merged[idx_s]

xen_sub <- xenium.sce[, xenium.sce[[condition_col]] %in% c("Healthy", "PSO")]
xen_sub[[condition_col]] <- droplevels(factor(xen_sub[[condition_col]]))
xen_sub[[sample_col]]    <- droplevels(factor(xen_sub[[sample_col]]))

## 2) (없으면) logcounts / PCA / UMAP
if (!"logcounts" %in% assayNames(xen_sub)) {
  logcounts(xen_sub) <- log1p(counts(xen_sub))
}
if (!"PCA" %in% reducedDimNames(xen_sub)) {
  xen_sub <- runPCA(xen_sub, ncomponents = 30, exprs_values = "logcounts")
}
if (!"UMAP" %in% reducedDimNames(xen_sub)) {
  xen_sub <- runUMAP(xen_sub, dimred = "PCA")
}

## 3) Milo object
xen_milo <- Milo(xen_sub)
reducedDim(xen_milo, "PCA")  <- reducedDim(xen_sub, "X_pca_harmony")
reducedDim(xen_milo, "UMAP") <- reducedDim(xen_sub, "X_umap")

## 4) Graph + neighbourhoods
xen_milo <- buildGraph(xen_milo, k = 30, d = 30)
xen_milo <- makeNhoods(xen_milo, prop = 0.1, k = 30, d = 30, refined = TRUE)

## 5) Count cells per nhood per sample
meta_df <- as.data.frame(colData(xen_milo))
xen_milo <- countCells(xen_milo, meta.data = meta_df, samples = sample_col)

## 6) Design matrix (sample-level)
design_df <- meta_df %>%
  dplyr::select(all_of(c(sample_col, condition_col))) %>%
  distinct()
rownames(design_df) <- design_df[[sample_col]]
design_df <- design_df[colnames(nhoodCounts(xen_milo)), , drop = FALSE]

## 7) DA test
xen_milo <- calcNhoodDistance(xen_milo, d = 30)

da_res <- testNhoods(
  xen_milo,
  design = stats::as.formula(paste0("~ ", condition_col)),
  design.df = design_df
)
da_res %>% arrange(SpatialFDR) %>% head

## 8) Plot (optional)
xen_milo <- buildNhoodGraph(xen_milo)
da_results <- annotateNhoods(xen_milo, da_res, coldata_col = "celltype2_2")
head(da_results[order(da_results$SpatialFDR), ])

outdir <- "MiloR"
png(file.path(outdir, "milo_nhood_DA.png"),
    width = 1800, height = 1600, res = 200)
plotNhoodGraphDA(xen_milo, da_res, alpha = 0.05)
dev.off()



umap_pl <- plotReducedDim(xen_milo, dimred = "UMAP", colour_by="condition", text_by = "celltype2_2", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")
nh_graph_pl <- plotNhoodGraphDA(xen_milo, da_results, layout="UMAP",alpha=0.1) 

png(file.path(outdir, "UMAP_milo_nhood_DA2_2.png"), width = 2900, height = 1600, res = 200)
umap_pl + nh_graph_pl +plot_layout(guides="collect")
dev.off()

keep_lvls <- intersect(celltype2_2_order, unique(da_results$celltype2_2))
da_results$celltype2_2 <- factor(da_results$celltype2_2, levels = rev(keep_lvls))

png(file.path(outdir, "DAbeeswarm_milo_nhood_DA2.png"), width = 1500, height = 2000, res = 200)
plotDAbeeswarm(da_results, group.by = "celltype2_2")
dev.off()

save(
  da_results,
  xen_milo,
  file =  "/MiloR_celltype2_analized_250106.Rdata"
)

saveRDS(xenium.sce,
        file =  "Xenium_singleRT_sce.Rds")

#------
celltype2_2_order <- c(
  "B cell",
  "Plasma cell",
  "CD4_Tn",
  "CD4_Tm",
  "CD8_Tm",
  "CD8_Tex",
  "Treg",
  "NK",
  "ILC",
  "Macrophage",
  "Neutrophil",
  "mDC",
  "cDC2",
  "cDC1",
  "pDC",
  "Ductal epithelial cell",
  "Basal KC",
  "Basal KC2",
  "Suprabasal KC",
  "Granular KC",
  "Proliferative cell",
  "Fibroblast_Inflamed",
  "Fibroblast_MMP1+",
  "Fibroblast–Plasma",
  "Vascular EC",
  "Perivascular EC",
  "Lymphatic EC",
  "Schwann cell–like",
  "Adipocyte"
)

