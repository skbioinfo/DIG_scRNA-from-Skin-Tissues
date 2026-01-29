# A function to perform data pre-processing and clustering

#' dimreduc_cluster_seurat: Perform data normalization, scaling, variable feature selection, and clustering in a single automated pipeline.
#'
#' Description: Subset a list of Seurat objects or a single Seurat object given a sample ID pattern to match.
#'
#' @param seurat A seurat object.
#' @param out_dir Project directory. This function will create the subfolders "PCA" and "UMAP" to store the necessary plots
#' @param object_name Name of the Seurat object (for file naming purposes). This could be written into the function but significantly increases runtime.
#' @param n_pcs Number of PCs used in dimensionality reduction
silhouette_seurat <- function(seurat, out_dir, object_name, n_pcs) {
    ## Check out_dir and create necessary subfolders ----
    if(!endsWith(out_dir, "/")) {
        out_dir <- paste0(out_dir, "/")
    }

    # Create subdirectory
    if(!dir.exists(paste0(out_dir, "silhouette"))) {
        dir.create(paste0(out_dir, "silhouette"))
    }

    ## Calculate average silhouette width (SilScore) across all resolutions in seurat ----
    # Pull cluster res from seurat object
    res_idents <- grep("RNA_snn_res", colnames(seurat@meta.data), value = TRUE)
    cluster_res <- substr(res_idents, 13, nchar(res_idents)) %>% as.numeric()  

    # Calculate average silhouette width
    SilScores <- data.frame(Res = numeric(), AvgSilWidth = numeric())
    dist.matrix <- dist(x = Embeddings(object = seurat[["pca"]])[, 1:n_pcs])
    for (res in cluster_res) {
        seurat <- SetIdent(seurat, value = paste0("RNA_snn_res.", res))
        clusters <- Idents(seurat)
        sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
        temp_sil_df <- data.frame(Res = res, AvgSilWidth = mean(sil[,3]))
        SilScores <- rbind(SilScores, temp_sil_df)
    }

    ## Generate line plot of avg sil width and save to pdf ----
    # Plot Silhouette scores as a function of leiden resolution
    SilPlot <- ggplot(SilScores, aes(x = Res, y = AvgSilWidth)) + 
    geom_line(size = 1) +
    geom_point(size = 3) +
    scale_x_continuous(name = "Leiden Resolution", breaks = cluster_res) +
    ggtitle(paste0(object_name, " Silhouette Plot")) +
    theme_bw()

    # Save plot to PDF
    pdf(paste0(out_dir , "/silhouette/", object_name, "_silscore.pdf"), width = 8, height = 6)
    print(SilPlot)
    dev.off()

    ## Generate output ----
    # Save best resolution
    best_res <- SilScores$Res[which.max(SilScores$AvgSilWidth)] %>% format(, nsmall = 2)
    return(list(SilScores = SilScores, best_res = best_res))
}
