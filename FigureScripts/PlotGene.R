plotGene <- function(counts, meta, gene, celltypes, cols, outdir){
	    .e <- environment()
    library(ggplot2)
        counts <- as.matrix(counts)
        counts <- counts[gene, which(meta$celltype %in% celltypes)]
	    meta <- meta[which(meta$celltype %in% celltypes), ]
	    means <- rep(0, length(unique(meta$celltype_condition)))
	        sds <- rep(0, length(unique(meta$celltype_condition)))
	        celltype <- rep("", length(unique(meta$celltype_condition)))
		    timepoint <- rep("", length(unique(meta$celltype_condition)))
		    index <- 1
		        for(i in unique(meta$celltype_condition)){
				        means[index] <- mean(counts[which(meta$celltype_condition == i)])
		            sds[index] <- sd(counts[which(meta$celltype_condition == i)])
			            celltype[index] <- meta$celltype[which(meta$celltype_condition == i)[1]]
			            timepoint[index] <- meta$Condition[which(meta$celltype_condition == i)[1]]
				            index <- index + 1
				        }
		        sds[is.na(sds)] <- 0
			    df <- data.frame(means = means, sds = sds, celltype = celltype, timepoint = timepoint)
			    df$celltype <- factor(df$celltype, levels = celltypes)
			        df$timepoint <- factor(df$timepoint, levels = c("Resting", "2","5","15"))
			        pdf(paste(outdir, "/", gene, ".pdf", sep = ""), width = 10, height = 7)
				    print(ggplot(data = df, aes(x = timepoint, y = means, group = celltype)) + 
					          geom_line(aes(color = celltype, linetype = celltype), size = 2) + 
						          geom_point(aes(color = celltype)) + 
							          geom_errorbar(aes(ymin = means - sds, ymax = means + sds, color = celltype), width = 0.1) + 
								          theme_classic(base_size = 25) + scale_color_manual(values = cols) +
									          scale_linetype_manual(values = c("solid","solid","solid","solid","solid","solid","solid","solid","solid")))
				    dev.off()
}
