\name{DetaDcorplot}
\alias{DetaDcorplot}
\title{Plot correlation matrix for pairwise deta D.}
\description{This function plots the pairwise genetic differentiation (deta D).}
\usage{DetaDcorplot(x, ncode)}
\arguments{
  \item{x}{The genetic data file, specifically genepop format.}
  \item{ncode}{The coding type of the genetic data.}
}
\details{ This function returns the pairwise genetic differentiation (deta D) plot.}
\value{The result will return a genetic differentiation (deta D) correlation plot.}
\references{
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{qinxinghu@gmail.com}
\seealso{HierD; CorDdplot}
\examples{
download.file("https://sites.google.com/site/advancesinevolutionandgenetics
/Home/unit-3-1/Island.gen?attredirects=0&d=1",
"C:/Users/Public/Documents/Island.gen",cacheOK=FALSE)
DetaDcorplot("C:/Users/Public/Documents/Island.gen",ncode=3)
}
\keyword{Correlation plot}
\keyword{Plot genetic differentiation}
