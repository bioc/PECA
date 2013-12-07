`PECA_tsv` <-
function(file=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t", type="median", paired=FALSE) {
	# Read tsv-file
	message("Reading data.")
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
	message("Reading data.")
	data <- ReadAffy(filenames=c(samplenames1,samplenames2))
	probenamesGene <- probeNames(data)
	probeintensities <- pm(data)
	rm(data)
	gc()
	PECA(probenamesGene=probenamesGene, probeintensities=probeintensities, samplenames1=samplenames1, samplenames2=samplenames2, normalize=normalize, test=test, type=type, paired=paired)
}

`PECA` <-
function(probenamesGene=NULL, probeintensities=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t", type="median", paired=FALSE) {

if (normalize) {
	message("Performing quantile-normalization.")
	probeintensities <- normalize.quantiles(probeintensities)
}

message("Performing log-transformation.")
probeintensities <- log2(probeintensities)
probeintensities[is.infinite(probeintensities)] <- 0
colnames(probeintensities) <- c(samplenames1,samplenames2)

# PECA slr and t-statistic
message("Calculating probe-level statistics.")
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

message("Calculating gene-level statistics.")  
n.gene <- tapply(t, probenamesGene, function(x) sum(!is.na(x)))  
if (type=="median") {
	geneSLR <- tapply(probeSLR, probenamesGene, median, na.rm=TRUE)
	t <- tapply(t, probenamesGene, median, na.rm=TRUE)
}
if (type=="tukey") {
	geneSLR <- tapply(probeSLR, probenamesGene, tukey)
	t <- tapply(t, probenamesGene, tukey)
}
p.gene <- 2 * pt(abs(t), df=df.total, lower.tail=FALSE)
p.beta.gene <- pbeta(p.gene, n.gene/2 + 0.5, n.gene - (n.gene/2 + 0.5) + 1)

## Cleanup
rm(probeintensities)
rm(probeSLR)
gc()

## Return a table containing slr, t-statistic and p-value
result <- cbind(slr=geneSLR, t=t, score=2*pt(abs(t), df=df.total, lower.tail=FALSE), p=p.beta.gene)
message("Done.")
return(result)

}
