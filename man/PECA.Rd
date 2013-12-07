\name{PECA}
\alias{PECA}
\alias{PECA_AffyBatch}
\alias{PECA_CEL}
\alias{PECA_tsv}
\title{PECA differential gene expression}
\description{
  Calculates the PECA ordinary or modified t-statistic to determine differential expression between two groups of samples in Affymetrix gene expression studies or peptide-based proteomic studies.
}

\details{
\code{PECA} determines differential gene expression using directly the probe-level measurements from Affymetrix gene expression microarrays or proteomic datasets. An expression change between two groups of samples is first calculated for each probe/peptide on the array. The gene/protein-level expression changes are then defined as medians over the probe-level changes. For more details about the probe-level expression change averaging (PECA) procedure, see Elo et al. (2005), Laajala et al. (2009) and Suomi et al.

\code{PECA} calculates the probe-level expression changes using the ordinary or modified t-statistic. The ordinary t-statistic is calculated using the function \code{rowttests} in the Bioconductor \code{genefilter} package. The modified t-statistic is calculated using the linear modeling approach in the Bioconductor \code{limma} package. Both paired and unpaired tests are supported.

The significance of an expression change is determined based on the analytical p-value of the gene-level test statistic. Unadjusted p-values are reported along with the corresponding p-values looked up from beta ditribution. The quality control and filtering of the data (e.g. based on low intensity or probe specificity) is left to the user.
}

\usage{
## Read AffyBatch object
PECA_AffyBatch(affy=NULL, normalize=FALSE, test="t", type="median", paired=FALSE)

## Read CEL-files
PECA_CEL(samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t",
	type="median", paired=FALSE)

## Read tab separated text file	
PECA_tsv(file=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t",
	type="median", paired=FALSE)
}
\arguments{
  \item{affy}{AffyBatch object}
  \item{normalize}{a logical indicating whether quantile normalization is applied}
  \item{test}{a character string indicating whether the ordinary ("\code{t}") or modified ("\code{modt}") t-test is performed.}
  \item{type}{a character string indicating whether ("\code{median}") or ("\code{tukey}") is used when calculating gene/protein values}
  \item{paired}{a logical indicating whether a paired test is performed}
  \item{file}{filename of tab separated data}
  \item{samplenames1}{a character vector containing the names of the .CEL-files/columns in the first group}
  \item{samplenames2}{a character vector containing the names of the .CEL-files/columns in the second group}
}

\examples{
## Load spike-in example
library(SpikeIn)
data(SpikeIn95)

## Subset the data
data <- SpikeIn95[,13:20]

## Run the test
results <- PECA_AffyBatch(affy=data)

## Show portion of the results
head(results)
}

\references{ 
T. Suomi, O. Nevalainen and L.L. Elo:
Evaluation of probe-level expression change averaging in microarray and proteomics data.
2013

L.L. Elo, L. Lahti, H. Skottman, M. Kylaniemi, R. Lahesmaa and T. Aittokallio:
Integrating probe-level expression changes across generations of Affymetrix arrays. 
Nucleic Acids Research 33(22), e193, 2005.

E. Laajala, T. Aittokallio, R. Lahesmaa and L.L. Elo: 
Probe-level estimation improves the detection of differential splicing in Affymetrix exon array studies.
Genome Biology 10(7), R77, 2009.

H. Bengtsson, K. Simpson, J. Bullard and K. Hansen:
aroma.affymetrix: A generic framework in R for analyzing small to very large Affymetrix data sets in bounded memory.
Tech Report \#745, Department of Statistics, University of California, Berkeley, 2008.
}
\keyword{math}
