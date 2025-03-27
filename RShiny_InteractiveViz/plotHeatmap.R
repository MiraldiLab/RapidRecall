plotHeatmap <- function(counts, meta, genes, celltypes) {
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(circlize)
    library(matrixStats)
    library(dplyr)
    library(ggplot2)

    # Ensure colnames are correct
    heat_col_tfmrna <- colorRamp2(c(-2, 0, 2), c('dodgerblue','white','red')) ## Colors for heatmap
    counts <- as.matrix(counts)
    colnames(counts) <- gsub("\\.", "/", colnames(counts))

    # Z-score counts
    counts <- t(scale(t(counts)))

    # Filter meta and counts by selected celltypes
    index <- which(meta$celltype %in% celltypes)
    meta <- meta[index, ]
    counts <- counts[, index]

    # Define celltype_timepoint
    celltype_timepoint <- c("Naive_Resting", "Naive_2", "Naive_5", "Naive_15",
                            "TCM_Resting", "TCM_2", "TCM_5", "TCM_15",
                            "Th1_Resting", "Th1_2", "Th1_5", "Th1_15",
                            "Th2_Resting", "Th2_2", "Th2_5", "Th2_15",
                            "Th17_Resting", "Th17_2", "Th17_5", "Th17_15",
                            "MHCII_Resting", "MHCII_2", "MHCII_5", "MHCII_15",
                            "CTL_Resting", "CTL_2", "CTL_5", "CTL_15",
                            "Treg_Resting", "Treg_2", "Treg_5", "Treg_15",
                            "TCM/TEM_2", "TCM/TEM_5", "TCM/TEM_15")
    celltype_timepoint <- celltype_timepoint[which(celltype_timepoint %in% meta$celltype_condition)]

    # Create Annotation
    annotation_matrix <- as.matrix(unique(rownames(meta)))
    rownames(annotation_matrix) <- annotation_matrix[, 1]
    index <- match(rownames(annotation_matrix), rownames(meta))
    annotation_matrix[, 1] <- (meta[, "celltype"])[index]
    annotation_matrix <- cbind(annotation_matrix, meta$celltype_condition[index])
    annotation_matrix <- cbind(annotation_matrix, meta$celltype_donor[index])
    annotation_matrix <- cbind(annotation_matrix, meta$Condition[index])
    index <- which(celltype_timepoint %in% annotation_matrix[, 2])
    celltype_timepoint <- celltype_timepoint[index]
    index <- which(annotation_matrix[, 2] %in% celltype_timepoint)
    annotation_matrix <- annotation_matrix[index, ]

    # Define colors
    heatmap_colors <- c("#E41A1C", "#4A72A6", "#48A462", "#7E6E85", "#D16948", "#E1C62F", "#B75F49", "#EC83BA", "#999999")
    heatmap_colors <- setNames(heatmap_colors, c("Naive", "TCM", "Th1", "Th2", "Th17", "MHCII", "CTL", "Treg", "TCM/TEM"))
    heatmap_colors <- heatmap_colors[which(names(heatmap_colors) %in% celltypes)]
    condition_colors <- c("#ebebeb", "#C0C0C0", "#A9A9A9", "#808080")
    names(condition_colors) <- c("Resting", "2", "5", "15")

    # Check if genes are present in counts
    missing_genes <- setdiff(genes, rownames(counts))
    if (length(missing_genes) > 0) {
        warning("The following genes are missing from the counts matrix: ", paste(missing_genes, collapse = ", "))
        genes <- setdiff(genes, missing_genes)
    }

    # Create HeatmapAnnotation
    ha <- HeatmapAnnotation(condition = annotation_matrix[, 4], celltype = annotation_matrix[, 1],
                            col = list(condition = condition_colors, celltype = heatmap_colors),
                            simple_anno_size = unit(1, "cm"),
                            annotation_legend_param = list(celltype = list(at = celltypes)))

    #splits <- unit(c(0, 0, 0, 1.5, 0, 0, 0, 1.5, 0, 0, 0, 1.5, 0, 0, 0, 1.5, 0, 0, 0, 1.5, 0, 0, 0, 1.5, 0, 0, 0, 1.5, 0, 0, 0, 1.5, 0, 0, 0), "mm")
    splits <- unit(rep(c(0,0,0,1.5), length(celltypes)), "mm")
    if("TCM/TEM" %in% celltypes){
        splits <- splits[-(length(splits)-2)]
    }
    # Create and draw the heatmap
    ht <- Heatmap(counts[genes, ], name = "Z-score", show_row_names = TRUE, show_row_dend = FALSE,
                  show_column_dend = FALSE, border = TRUE, show_column_names = FALSE,
                  column_split = factor(annotation_matrix[, 2], levels = celltype_timepoint),
                  col = heat_col_tfmrna, column_gap = splits, row_gap = unit(0, "mm"),
                  cluster_columns = FALSE, top_annotation = ha, row_names_gp = gpar(fontsize = 16),
                  use_raster = FALSE, column_title = NULL, row_title = NULL)
    return(ht)
}

