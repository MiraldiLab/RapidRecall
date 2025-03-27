plotPattern <- function(cluster, meta, counts, out){
    meta$celltype[which(meta$celltype %in% c("Th1","Th2","Th17"))] <- "TEM"
    meta$celltype_condition <- paste(meta$celltype, meta$Condition, sep = "_")
    celltype_condition <- c("Naive_Resting","Naive_2","Naive_5","Naive_15","TCM_Resting","TCM_2","TCM_5","TCM_15","TEM_Resting","TEM_2","TEM_5","TEM_15")
    celltype <- c("Naive","Naive","Naive","Naive","TCM","TCM","TCM","TCM","TEM","TEM","TEM","TEM")
    condition <- c("Resting","2","5","15","Resting","2","5","15","Resting","2","5","15")
    means <- c()
    sds <- c()
    for(i in unique(celltype_condition)){
        curr_average <- rowMeans(counts[cluster, which(meta$celltype_condition == i)])
        curr_sds <- apply(counts[cluster, which(meta$celltype_condition == i)], 1, sd, na.rm=TRUE)
        means <- c(means, mean(curr_average))
        sds <- c(sds, mean(curr_sds))
    }

    df <- data.frame(Counts = means, Std = sds, Celltype = celltype, Condition = condition)
    df$Condition <- factor(df$Condition, levels = c("Resting","2","5","15"))
    df$Celltype <- factor(df$Celltype, levels = c("Naive","TCM","TEM"))
    cols <- c("#E41A1C","#4A72A6","#48A462")
    pdf(out, width = 10, height = 7)
    print(ggplot(data = df, aes(x = Condition, y = Counts, group = Celltype)) + 
        geom_line(aes(color = Celltype), linewidth = 3.2) + 
        #geom_errorbar(aes(ymin = Counts - Std, ymax = means + Std, color = celltype), width = 0.1) +
        geom_point(aes(color = Celltype)) + 
        theme_classic(base_size = 25) + 
        scale_color_manual(values = cols))
    dev.off()
}


