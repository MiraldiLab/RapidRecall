
plotPattern <- function(cluster, meta, counts, out){
	    #meta$celltype[which(meta$celltype %in% c("Th1","Th2","Th17"))] <- "TEM"
	    meta$celltype_condition <- paste(meta$celltype, meta$Condition, sep = "_")
    celltype_condition <- c("Naive_Resting","Naive_2","Naive_5","Naive_15","TCM_Resting","TCM_2","TCM_5","TCM_15","TCM_5","TCM_Proliferation_15","Th1_Resting","Th1_2","Th1_5","Th1_15",
			                                "Th2_Resting","Th2_2","Th2_5","Th2_15","Th17_Resting","Th17_2","Th17_5","Th17_15","TEM_act_Resting","TEM_act_2","TEM_act_5","TEM_act_15",
							                            "TEM_act2_Resting","TEM_act2_2","TEM_act2_5","TEM_act2_15","TCM/TEM_2","TCM/TEM_5","TCM/TEM_15","Treg_Resting","Treg_2","Treg_5","Treg_15")
        celltype <- c("Naive","Naive","Naive","Naive","TCM","TCM","TCM","TCM","TCM_Proliferation","TCM_Proliferation","Th1","Th1","Th1","Th1","Th2","Th2","Th2","Th2","Th17","Th17","Th17","Th17",
		                      "TEM_act","TEM_act","TEM_act","TEM_act","TEM_act2","TEM_act2","TEM_act2","TEM_act2","TCM/TEM","TCM/TEM","TCM/TEM","Treg","Treg","Treg","Treg")
        condition <- c("Resting","2","5","15","Resting","2","5","15","5","15","Resting","2","5","15","Resting","2","5","15","Resting","2","5","15",
		                           "Resting","2","5","15","Resting","2","5","15","2","5","15","Resting","2","5","15")
	    means <- c()
	    sds <- c()
	        for(i in celltype_condition){
			        if(length(which(meta$celltype_condition == i)) > 1){
					            curr_average <- rowMeans(counts[cluster, which(meta$celltype_condition == i)])
	                curr_sds <- apply(counts[cluster, which(meta$celltype_condition == i)], 1, sd, na.rm=TRUE)
			            means <- c(means, mean(curr_average))
			            sds <- c(sds, mean(curr_sds))
				            }
	            else{
			                curr_average <- counts[cluster, which(meta$celltype_condition == i)]
		                curr_sds <- 0
				            means <- c(means, mean(curr_average))
				            sds <- c(sds, mean(curr_sds))
					            }
		        }

	        df <- data.frame(Counts = means, Std = sds, Celltype = celltype, Condition = condition)
		    df$Condition <- factor(df$Condition, levels = c("Resting","2","5","15"))
		    #df$Celltype <- factor(df$Celltype, levels = c("Naive","TCM","TCM_Proliferation","Th1","Th2","Th17","TEM_act","TEM_act2","TCM/TEM","Treg"))
		    df$Celltype <- factor(df$Celltype, levels = rev(c("TCM_Proliferation","TCM","TCM/TEM","Naive","Th1","Th2","Th17","TEM_act","TEM_act2","Treg")))
		        cols <- rev(c("red","#4A72A6","#E1C62F","#666666","#666666","#666666","#666666","#666666","#666666","#666666"))
		        pdf(out, width = 5, height = 4, compress = F)
			    print(ggplot(data = df, aes(x = Condition, y = Counts, group = Celltype)) + 
				          geom_line(aes(color = Celltype, linetype = Celltype, size = Celltype, alpha = Celltype),show.legend = FALSE) + 
					          #geom_errorbar(aes(ymin = Counts - Std, ymax = means + Std, color = celltype), width = 0.1) +
					          geom_point(aes(color = Celltype),show.legend = FALSE, size=2.5) + 
						          theme_classic(base_size = 25) +
							          scale_linetype_manual(values = rev(c(2,1,1,3,1,1,1,1,1,1))) +
								          scale_size_manual(values = rev(c(2,2,2,1,1,1,1,1,1,1))) + 
									          scale_alpha_manual(values = rev(c(1,1,1,0.4,0.4,0.4,0.4,0.4,0.4,0.4))) +
										          scale_color_manual(values = cols)) 
			    dev.off()
}
