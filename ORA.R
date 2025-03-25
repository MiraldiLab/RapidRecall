
performOraHypergeometric <- function(pathway, queryGenes, backgroundGenes = NULL, minSize = NULL, maxSize = NULL) {
	    
	    # if (is.null(minSize)){
	    #     minSize = 5
	    # }
	    
	    # if (is.null(backgroundGenes)){
	    #     backgroundGenes <- unique(unlist(pathway, use.names = T))
	    # }
	    
	    if (is.null(minSize) || minSize < 5) {
		            warning('minSize is less 5 or NULL. Minimun Geneset/Pathway Size will Default to 5')
        minSize <- 5
	    }
    if (is.null(backgroundGenes)){
	            backgroundGenes <- unique(unlist(pathway, use.names = T))
        }

        if (is.null(maxSize) || maxSize > length(backgroundGenes) - 1) {
		        maxSize <- length(backgroundGenes) - 1
	    }
        
        pathwaySizes <- sapply(pathway, length)
	    toKeep <- which(minSize <= sapply(pathway, length) & sapply(pathway, length) <= maxSize)
	    
	    if (length(toKeep) ==  0){
		            return(data.frame())
	        }
	        pathway <- pathway[toKeep]
	        pathwaySizes <- pathwaySizes[toKeep]
		    
		    queryGenesFiltered <- intersect(queryGenes, backgroundGenes)
		    N <- length(backgroundGenes)  # Total number of background genes. This could be all genes in the pathway. 
		        n <- length(queryGenesFiltered)  # Number of genes in the gene set/set of interest
		        pvalues <- c()
			    
			    for (ix in seq_along(pathway)) {
				            K <- length(pathway[[ix]])  # Number of genes in the current pathway
			        overlap <- length(intersect(queryGenesFiltered, pathway[[ix]]))  # Overlap between gene set and pathway
				        p_value <- phyper(overlap - 1, K, N - K, n, lower.tail = FALSE)
				        pvalues <- c(pvalues, p_value)
					    }
			    
			    padj <- p.adjust(pvalues, method = "BH") # Multiple testing correction
			        oraRes <- data.frame(Pathway = names(pathway), pAdjusted = padj)

			        return(oraRes)
}
