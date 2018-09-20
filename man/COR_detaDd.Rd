\name{COR_detaDd}
\alias{COR_detaDd}
\title{Function to calculate pairwise detaD and correlation between genetic differentiation (deta D) and geographic distance}
\description{This function calculates new pairwise genetic differentication, deta D (Oscar E, et al, 2018), as well as the correlation between genetic differentiation (deta D) and geographic distance.}
\usage{COR_detaDd(f, d, ncode)}
\arguments{
  \item{f}{A genetic data files, here the format is genepop format.}
  \item{d}{Pairwise geographic distance for calculating correlations. If no actual geographic distance, d=FALSE. If TRUE,the geographic distance should be a matrix that corresponds to the pairwise genetic differentiation matrix.}
  \item{ncode}{The code type of the genotype}
}
\details{This function returns lists of pairwise detaD matrix, and a matrix of theorical geographical distance if d is FALSE, as well as their correlation coefficient.}
\value{
  \item {PairwiseDetaD}{The pairwise genetic differentiation-detaD matrix}
  \item {Dgeo}{The pairwise geographical distance matrix}
  \item {CordetaDd}{Person correlation coefficient between genetic differentiation (deta D) and geographic distance}
}
\references{
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{qinxinghu@gmail.com}
\seealso{COR_Fst}
\examples{
 d1=COR_detaDd("Island.gen",d=FALSE,ncode=3)
 print(d1)}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
