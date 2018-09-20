\name{CorDdPlot}
\alias{CorDdPlot}
\title{ Plot the relation between genetic differentiation (deta D) and geographic distance.}
\description{This function plots the linear regression of genetic differentiation (deta D) against geographic distance.}
\usage{CorDdPlot(x, d, ncode)}
\arguments{
  \item{x}{The genetic file, here is genepop format file.}
  \item{d}{Geographic distance.If D is TRUE, you should input the matrix of geographic distance; if no realistic distance,d should be set as FALSE, then the default distance will be generated according to the order of the population sequence.}
  \item{ncode}{ The type of your coding allele}
}
\details{}
\value{
  \item{plot}{Return to a linear plot for the relation between genetic differentiation (deta D) and geographic distance.}
  \item{lm}{ Return to details of linear regression between genetic differentiation (deta D) and geographic distance.}
}
\references{
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{qinxinghu@gmail.com}
\note{}
\seealso{}
\examples{
CorDdPlot("Island.gen",d=FALSE,ncode=3)}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
