\name{COR_Fstd}
\alias{COR_Fstd}
\title{Function to calculate pairwise Fst and correlation coefficient between Fst and geographic distance.}
\description{This function calculates pairwise genetic differentication (Fst), as well as the correlation        coefficient between genetic differentiation (Fst) and geographic distance.}
\usage{COR_Fstd(x, d, ncode)}
\arguments{
  \item{x}{A genetic data files, here the format is genepop format.}
  \item{d}{Pairwise geographic distance for calculating correlations. If no actual geographic distance, d=FALSE. The geographic distance should be a matrix that corresponds to the pairwise genetic differentiation matrix.}
  \item{ncode}{The code type of the genotype.}
}
\details{ This function returns lists of pairwise Fst matrix, and a matrix of theorical geographical distance if d is FALSE, as well as their correlation coefficient.}
\value{
  \item{pwFst}{The pairwise genetic differentiation (Fst) matrix.}
  \item{COR_Fstd}{Person correlation coefficient between genetic differentiation (Fst) and geographic distance.}
}
\references{
Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5(1), 184-186.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24 (11), 1403-1405.

Weir, B. S.; Cockerham, C. Clark (1984). "Estimating F-Statistics for the Analysis of Population Structure". Evolution. 38 (6): 1358.

Evolution and the Genetics of Populations Volume 2: the Theory of Gene Frequencies, pg 294-295, S. Wright, Univ. of Chicago Press, Chicago, 1969.}
\author{ qinxinghu@gmail.com}
\seealso{COR_detaDd}
\examples{
infile <- system.file('extdata', 'Island.gen',package='HierDpart')
locinfile <- 'Island.gen' ## file in user's directory not in R's extdata directory
check <- file.copy(infile,locinfile,overwrite=TRUE)
d2=COR_Fstd(infile,d=FALSE,ncode=3)
print(d2)
}
\keyword{Correlation between genetic differentiation and geographic distance}

