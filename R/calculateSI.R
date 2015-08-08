#------------------------------
# SPLICING INDEX CALCULATION
#------------------------------

calculateSI <- function(geneEffect, exonEffect, genenames, probenamesGene) {

	if(ncol(geneEffect) != ncol(exonEffect)) stop("Matrices do not match.")

	# Sort the probenames and the exonEffect matrix.
	temp <- sort.int(genenames, index.return=TRUE)
	sortedGenenames <- temp$x
	geneEffect <- geneEffect[temp$ix,]
	temp2 <- sort.int(probenamesGene, index.return=TRUE)
	sortedProbenamesGene <- temp2$x
	exonEffect <- exonEffect[temp2$ix,]

	# Pick the first instance of each genename.
	loc <- match(sortedGenenames, sortedProbenamesGene)
	loc[length(loc) + 1] <- length(sortedProbenamesGene) + 1 

	# Transpose the matrix to calculate the difference between a vector and a matrix.
	exonEffect <- t(exonEffect)

    ## Probe-level splicing index values.
    SI <- matrix(nrow=nrow(exonEffect), ncol=ncol(exonEffect))
    cumMatch <- NULL
    pb <- txtProgressBar(min=0, max=length(loc)-1, style=3)
    for(i in 1:(length(loc)-1) ) {
      setTxtProgressBar(pb, i)
      whMatch <- loc[i]:(loc[i+1]-1)
      cumMatch <- c(cumMatch, whMatch)
      if(length(whMatch)==0) {
        warning("Gene ",i," (",sortedGenenames[i],") is not found in the exonEffect matrix.")
        next
      }
      SI[,whMatch] <- exonEffect[,whMatch] - geneEffect[i,]
    }
    close(pb)
	
	if(length(cumMatch) != ncol(exonEffect)) {
		warning("Some genes in the exonEffect matrix were not found in the geneEffect matrix.")
	}
	
	return(list(SI=t(SI), ix=temp2$ix))

}

