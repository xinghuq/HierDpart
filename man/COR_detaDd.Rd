\name{COR_detaDd}
\alias{COR_detaDd}
\title{Function to calculate pairwise detaD and correlation between genetic differentiation (deta D) and geographic distance}
\description{This function calculates new pairwise genetic differentication, deta D (Oscar E, et al, 2018), as well as the correlation between genetic differentiation (deta D) and geographic distance.}
\usage{COR_detaDd(f, d, ncode)}
\arguments{
  \item{f}{A genetic data files, the format is genepop format.}
  \item{d}{Pairwise geographic distance for calculating correlations. If no actual geographic distance, d=FALSE. If TRUE, the geographic distance should be a matrix that corresponds to the pairwise genetic differentiation matrix.}
  \item{ncode}{The coding type of your genetic data}
}
\details{This function returns to lists of pairwise deta D matrix, and a matrix of theorical geographical distance if d is FALSE, as well as their correlation coefficient.}
\value{
  \item{PairwiseDetaD}{The pairwise genetic differentiation-detaD matrix.}
  \item{Dgeo}{The pairwise geographical distance matrix.}
  \item{CordetaDd}{Person correlation coefficient between genetic differentiation (deta D) and geographic distance.}
}
\references{
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{qinxinghu@gmail.com}
\seealso{COR_Fst}
\examples{
infile <- system.file('extdata', 'Island.gen',package='HierDpart')
locinfile <- 'Island.gen' ## file in user's directory not in R's extdata directory
check <- file.copy(infile,locinfile,overwrite=TRUE)
d1=COR_detaDd(infile,d=FALSE,ncode=3)
print(d1)}
\keyword{Correlation between genetic differentiation and geographic distance}

