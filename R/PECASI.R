`PECASI` <- function(path, dataFolder, chipType, cdfTag=NULL, samplenames1, samplenames2, test="t") {

setwd(path)

if(length(samplenames1) != length(samplenames2)) print("The samplenames are not paired.")

# Annotation data
message("Reading annotations")
flush.console()
cdf <- AffymetrixCdfFile$byChipType(chipType, tags=cdfTag) 

# Expression data and quantile normalization
message("Reading data")
flush.console()
cs <- AffymetrixCelSet$byName(dataFolder, cdf=cdf)
	
message("Performing quantile normalization")
flush.console()
qn <- QuantileNormalization(cs, typesToUpdate="pm")
csN <- process(qn)

# The cell information
cells <- getCellIndices(cdf, units=1:nbrOfUnits(cdf))
cells <- flattenCellIndices(cells)
probenamesGene <- cells$unitNames
probenamesExon <- cells$groupNames

# Probe intensity matrix (log-transformed)
message("Performing log-transformation")
flush.console()
probeintensities <- getData(csN, indices=cells$cell, fields=c("intensities"))$intensities
probeintensities <- log2(probeintensities)
colnames(probeintensities) <- cs$names
rm(cells)
gc()

# Probe-level signal log-ratio (slr) for each sample pair
probeSLR<- matrix(nrow=nrow(probeintensities), ncol=length(samplenames1)) 
for(i in 1:length(samplenames1)) probeSLR[,i]<- probeintensities[,samplenames1[i]] - probeintensities[,samplenames2[i]]
rm(probeintensities)
gc()

# PECA signal log-ratio (slr) for each sample pair
geneSLR<- tapply(probeSLR[,1], factor(probenamesGene), median, na.rm=TRUE)
for(i in 2:length(samplenames1)) geneSLR<- cbind(geneSLR, tapply(probeSLR[,i], probenamesGene, median, na.rm=TRUE))

# Probe-level SI-values and sorting of the data
message("Calculating splicing index")
flush.console()
probeSI <- calculateSI(geneEffect=geneSLR, exonEffect=probeSLR, genenames=row.names(geneSLR), probenamesGene=as.character(probenamesGene))
probenamesExon <- probenamesExon[probeSI$ix]
probeSI <- probeSI$SI

message("Calculating low-level statistics")
flush.console()
if(test == "t") {
    # PECA SI and t-statistic
    t <- rowttests(probeSI, tstatOnly=TRUE)
    probeSI <- t$dm
    t <- t$statistic
    df.total <- length(samplenames1)-1
}

if(test == "modt") {
    # PECA SI and modified t-statistic
    fit <- lmFit(probeSI)
    fit <- eBayes(fit)
    probeSI <- fit$coef[,1]
    t <- fit$t[,1]
    df.total<- fit$df.residual[1] + fit$df.prior
	rm(fit)
}

# Aggregating t-statistics
message("Aggregating statistics")
flush.console()
exonSI <- tapply(probeSI, probenamesExon, median, na.rm=TRUE)
t <- tapply(t, probenamesExon, median, na.rm=TRUE)
p <- 2*pt(abs(t), df=df.total, lower.tail=FALSE)
rm(probeSI)
gc()

# Return a table containing si, t-statistic and p-value
result <- data.frame(cbind(si=exonSI, t=t, p=p))
return(result)

}


