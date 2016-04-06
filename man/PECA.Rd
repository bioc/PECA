\name{PECA}

\alias{PECA}
\alias{PECA_AffyBatch}
\alias{PECA_CEL}
\alias{PECA_tsv}
\alias{PECA_df}

\title{PECA differential gene expression}

\description{
  Calculates the PECA ordinary or modified t-statistic to determine differential expression between two groups of samples in Affymetrix gene expression studies or peptide-based proteomic studies.
}

\usage{
## Read AffyBatch object
PECA_AffyBatch(affy=NULL, normalize=FALSE, test="t", type="median", paired=FALSE, progress=FALSE)

## Read CEL-files
PECA_CEL(samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t",
   type="median", paired=FALSE, progress=FALSE)

## Read tab separated text file	
PECA_tsv(file=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE,
   test="t", type="median", paired=FALSE, progress=FALSE)

## Read dataframe	
PECA_df(df=NULL, id=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE,
   test="t", type="median", paired=FALSE, progress=FALSE)
}
\arguments{
  \item{affy}{AffyBatch object.}
  \item{normalize}{A character string indicating if ("\code{quantile}") or ("\code{median}") normalization is performed.}
  \item{test}{A character string indicating whether the ordinary ("\code{t}") or modified ("\code{modt}") t-test is performed.}
  \item{type}{A character string indicating whether ("\code{median}") or ("\code{tukey}") is used when calculating gene/protein values.}
  \item{paired}{A logical indicating whether a paired test is performed.}
  \item{file}{Filename of tab separated data.}
  \item{samplenames1}{A character vector containing the names of the .CEL-files/columns in the first group.}
  \item{samplenames2}{A character vector containing the names of the .CEL-files/columns in the second group.}
  \item{df}{Dataframe to be used as an input.}
  \item{id}{Column name of dataframe used for aggregating results.}
  \item{progress}{A logical indicating whether a progress bar is shown.}
}

\value{
  \code{PECADE} returns a matrix which rows correspond to the genes under analysis and columns indicate the corresponding signal log-ratio (slr), t-statistic, p-value and FDR adjusted p-value.
}

\details{
\code{PECA} determines differential gene expression using directly the probe-level measurements from Affymetrix gene expression microarrays or proteomic datasets. An expression change between two groups of samples is first calculated for each probe/peptide on the array. The gene/protein-level expression changes are then defined as medians over the probe-level changes. For more details about the probe-level expression change averaging (PECA) procedure, see Elo et al. (2005), Laajala et al. (2009) and Suomi et al.

\code{PECA} calculates the probe-level expression changes using the ordinary or modified t-statistic. The ordinary t-statistic is calculated using the function \code{rowttests} in the Bioconductor \code{genefilter} package. The modified t-statistic is calculated using the linear modeling approach in the Bioconductor \code{limma} package. Both paired and unpaired tests are supported.

The significance of an expression change is determined based on the analytical p-value of the gene-level test statistic. Unadjusted p-values are reported along with the corresponding p-values looked up from beta ditribution. The quality control and filtering of the data (e.g. based on low intensity or probe specificity) is left to the user.
}

\examples{
## Generate example data frame
df <- data.frame(id=c(rep("a",10),rep("b",10),rep("c",10)))
df$A1 <- rnorm(30, mean=50, sd=5)
df$A2 <- rnorm(30, mean=48, sd=5)
df$A3 <- rnorm(30, mean=50, sd=5)
df$B1 <- rnorm(30, mean=52, sd=5)
df$B2 <- rnorm(30, mean=53, sd=5)
df$B3 <- rnorm(30, mean=51, sd=5)

## Run the test
group1 <- c("A1","A2","A3")
group2 <- c("B1","B2","B3")
results <- PECA_df(df, group1, group2, id=id)
}

\references{ 
T. Suomi, G.L. Corthals, O. Nevalainen and L.L. Elo:
Using peptide-level proteomics data for detecting differentially expressed proteins.
2015

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
