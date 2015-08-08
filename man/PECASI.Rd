\name{PECASI}

\alias{PECASI}

\title{PECA splicing index}

\description{
  Calculates the PECA splicing index to determine differentially spliced exons between two groups of samples in Affymetrix exon array studies.
}

\usage{
PECASI(path, dataFolder, chipType, cdfTag=NULL, samplenames1, samplenames2, test="t")
}
\arguments{
  \item{path}{A character string specifying the path of the working directory containing the expression and annotation data.}
  \item{dataFolder}{A character string specifying the name of the directory containing the raw expression data (.CEL-files).}
  \item{chipType}{A character string specifying the microarray (chip) type.}
  \item{cdfTag}{A character string indicating an optional suffix added to the name of the particular chip definition file (CDF).}
  \item{samplenames1}{A character vector containing the names of the .CEL-files in the first group without the extension .CEL.}
  \item{samplenames2}{A character vector containing the names of the .CEL-files in the second group without the extension .CEL. The paired samples are assumed to be in the same order in both of the vectors \code{samplenames1} and \code{samplenames2}.}
  \item{test}{A character string indicating whether the ordinary ("\code{t}") or modified ("\code{modt}") t-test is performed.}
}

\value{
  \code{PECASI} returns a matrix which rows correspond to the exons under analysis and columns indicate the corresponding splicing index (si), t-statistic, p-value and FDR adjusted p-value.
}

\details{
\code{PECASI} determines differential alternative splicing using directly the probe-level measurements from Affymetrix exon microarrays. Differential splicing between two groups of samples is first calculated for each probe on the array. The exon-level differential splicing is then defined as the median over the probe-level differences. For more details about the probe-level expression change averaging (PECA) procedure, see Elo et al. (2005), Elo et al. (2006) and Laajala et al.

The current implementation of \code{PECASI} calculates the probe-level differential splicing using the ordinary or modified t-statistic over splicing index values. The ordinary t-statistic is calculated using the function \code{rowttests} in the Bioconductor \code{genefilter} package. The modified t-statistic is calculated using the linear modeling approach in the Bioconductor \code{limma} package. The samples are assumed to be paired. For more details about the PECA splicing index procedure, see Laajala et al.

\code{PECASI} uses the \code{aroma.affymetrix} package to normalize and extract the probe-level data from the .CEL-files (Bengtsson et al. 2008). Therefore, it is important that the naming and structure of the data files follow exactly the rules specified in the \code{aroma.affymetrix} package.

The raw expression data (.CEL-files) need to be in the directory \code{rawData/<dataFolder>/<chipType>}, where \code{rawData} is a directory under the current working directory specified by the \code{path}, \code{dataFolder} is the name of the dataset given by the user, and \code{chipType} indicates the type of the microarray used in the experiment.

In addition to the expression data, a chip definition file (CDF) is required. The CDF-file(s) for a particular microarray type \code{chipType} need to be in the directory \code{annotationData/chipTypes/<chipType>}, where \code{annotationData} is a directory under the current working directory specified by the \code{path}. Besides the CDF-files provided by Affymetrix, various custom CDF-files are available for a particular microarray type. The different versions can be separated by adding a suffix \code{cdfTag} to the name of the CDF-file: \code{<chipType>,<cdfTag>.cdf}

The quality control and filtering of the data (e.g. based on low intensity or probe specificity) is left to the user.
}

\references{ 
L.L. Elo, L. Lahti, H. Skottman, M. Kylaniemi, R. Lahesmaa and T. Aittokallio:
Integrating probe-level expression changes across generations of Affymetrix arrays. 
Nucleic Acids Research 33(22), e193, 2005.

L.L. Elo, M. Katajamaa, R. Lund, M. Oresic, R. Lahesmaa and T. Aittokallio:
Improving identification of differentially expressed genes by integrative analysis of Affymetrix and Illumina arrays. 
OMICS A Journal of Integrative Biology 10(3), 369--380, 2006.

E. Laajala, T. Aittokallio, R. Lahesmaa and L.L. Elo: 
Probe-level estimation improves the detection of differential splicing in Affymetrix exon array studies.
Genome Biology 10(7), R77, 2009.

H. Bengtsson, K. Simpson, J. Bullard and K. Hansen:
aroma.affymetrix: A generic framework in R for analyzing small to very large Affymetrix data sets in bounded memory.
Tech Report \#745, Department of Statistics, University of California, Berkeley, 2008.
}

\seealso{\code{\link{PECA}}}

\keyword{math}
