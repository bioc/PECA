`PECA_df` <-
function(df=NULL, id=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t", type="median", paired=FALSE) {
	# Read dataframe
	probeintensities <- df
	probenamesGene <- subset(probeintensities,select=id)
	probeintensities <- subset(probeintensities,select=c(samplenames1,samplenames2))
	probeintensities <- as.matrix(probeintensities)
	PECA(probenamesGene=probenamesGene, probeintensities=probeintensities, samplenames1=samplenames1, samplenames2=samplenames2, normalize=normalize, test=test, type=type, paired=paired)
}

`PECA_tsv` <-
function(file=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t", type="median", paired=FALSE) {
	# Read tsv-file
	message("Reading data.")
	flush.console()
	probeintensities <- read.csv(file, sep="\t")
	probenamesGene <- probeintensities[,1]
	probeintensities <- subset(probeintensities,select=c(samplenames1,samplenames2))
	probeintensities <- as.matrix(probeintensities)
	PECA(probenamesGene=probenamesGene, probeintensities=probeintensities, samplenames1=samplenames1, samplenames2=samplenames2, normalize=normalize, test=test, type=type, paired=paired)
}

`PECA_AffyBatch` <-
function(affy=NULL, normalize=FALSE, test="t", type="median", paired=FALSE) {
	# Read AffyBatch
	data <- affy
	probenamesGene <- probeNames(data)
	probeintensities <- pm(data)
	l <- length(sampleNames(data))
	sn1 <- sampleNames(data)[1:(l/2)]
	sn2 <- sampleNames(data)[(l/2+1):l]
	rm(data)
	gc()
	PECA(probenamesGene=probenamesGene, probeintensities=probeintensities, samplenames1=sn1, samplenames2=sn2, normalize=normalize, test=test, type=type, paired=paired)
}

`PECA_CEL` <-
function(samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t", type="median", paired=FALSE) {
	# Read affymetrix CEL-files
	message("Reading data")
	flush.console()
	data <- ReadAffy(filenames=c(samplenames1,samplenames2))
	probenamesGene <- probeNames(data)
	probeintensities <- pm(data)
	rm(data)
	gc()
	PECA(probenamesGene=probenamesGene, probeintensities=probeintensities, samplenames1=samplenames1, samplenames2=samplenames2, normalize=normalize, test=test, type=type, paired=paired)
}

# PECA function called by different wrappers
`PECA` <-
function(probenamesGene=NULL, probeintensities=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t", type="median", paired=FALSE) {

# Normalization
if (normalize==TRUE | normalize=="quantile") {
	message("Performing quantile normalization")
	flush.console()
	probeintensities <- normalize.quantiles(probeintensities)
}
if (normalize=="median") {
	message("Performing median normalization")
	flush.console()
	probeintensities <- normalizeMedianValues(probeintensities)
}

# Log transformation
message("Performing log-transformation")
flush.console()
probeintensities[probeintensities == 0] <- 1
probeintensities <- log2(probeintensities)
#probeintensities[is.infinite(probeintensities)] <- 0
colnames(probeintensities) <- c(samplenames1,samplenames2)

# PECA slr and t-statistic
message("Calculating low-level statistics")
flush.console()
if (test == "t") {
    if (paired) {
		probeSLR <- matrix(nrow=nrow(probeintensities), ncol=length(samplenames1)) 
		for(i in 1:length(samplenames1)) probeSLR[,i] <- probeintensities[,samplenames1[i]] - probeintensities[,samplenames2[i]]
        t <- rowttests(probeSLR, tstatOnly=TRUE)
		df.total <- length(samplenames1)-1
    }
    else {
        labels <- factor(c(rep(1,length(samplenames1)), rep(2,length(samplenames2))))    
		t <- rowttests(as.matrix(cbind(probeintensities[,samplenames1], probeintensities[,samplenames2])), fac=labels, tstatOnly=TRUE)
		df.total <- length(samplenames1)+length(samplenames2)-2
    }
    probeSLR <- t$dm
    t <- t$statistic
}
if (test == "modt") {
    if (paired) {
		probeSLR <- matrix(nrow=nrow(probeintensities), ncol=length(samplenames1)) 
		for(i in 1:length(samplenames1)) probeSLR[,i] <- probeintensities[,samplenames1[i]] - probeintensities[,samplenames2[i]]
		# probeSLR[is.infinite(probeSLR)] <- NA
		fit <- lmFit(probeSLR)
		fit <- eBayes(fit)
		probeSLR <- fit$coefficients
		t <- fit$t
    }
    else {
        design <- cbind(G1=1,G1vsG2=c(rep(1,length(samplenames1)), rep(0,length(samplenames2))))
		probeSLR <- as.matrix(cbind(probeintensities[,samplenames1], probeintensities[,samplenames2]))
		# probeSLR[is.infinite(probeSLR)] <- NA
        fit <- lmFit(probeSLR, design)
		fit <- eBayes(fit)
		probeSLR <- fit$coefficients[,2]
		t <- fit$t[,2]
    }
	df.total<- fit$df.residual[1] + fit$df.prior
	rm(fit)
	gc()
}
if (test == "rots") {
  #SV_ROTS<-c(0,0,0,1,1,1)
  #ROTS.out <- ROTS(data=probeintensities,groups=SV_ROTS,seed=1,B=1000,K=10000)
  #p <- ROTS.out$Pvalue1
}

# Aggregating t-statistics
message("Aggregating statistics")
flush.console()
gene.n <- tapply(t, probenamesGene, function(x) sum(!is.na(x)))  
if (type=="median") {
	geneSLR <- tapply(probeSLR, probenamesGene, median, na.rm=TRUE)
	t <- tapply(t, probenamesGene, median, na.rm=TRUE)
}
if (type=="tukey") {
	geneSLR <- tapply(probeSLR, probenamesGene, tukey)
	t <- tapply(t, probenamesGene, tukey)
}

# P-values
gene.p <- 2 * pt(abs(t), df=df.total, lower.tail=FALSE)
gene.p2 <- gene.p
if (type=="median") {
	gene.p2 <- pbeta(gene.p, gene.n/2 + 0.5, gene.n - (gene.n/2 + 0.5) + 1)
}
if (type=="tukey") {
	message("Simulating distributions")
	flush.console()
	distributions <- generateDistributions(max=max(gene.n), k=10000)
	gene.p2 <- mapply(psim, p=gene.p, list=distributions[gene.n])
}
gene.p.fdr <- p.adjust(gene.p2, method="fdr")

# Cleanup
rm(probeintensities)
rm(probeSLR)
gc()

# Return a table containing slr, t-statistic and p-value
result <- data.frame(cbind(slr=geneSLR, t=t, score=gene.p, p=gene.p2, p.fdr=gene.p.fdr))
message("Done")
return(result)

}
