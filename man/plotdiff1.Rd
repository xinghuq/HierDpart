\name{plotdiff1}
\alias{plotdiff1}
\title{Plot genetic differentiation (deta D) across loci}
\description{Plot genetic differentiation across (deta D, Oscar E, et al, 2018) loci.}
\usage{plotdiff1(x, ncode)}
\arguments{
  \item{x}{Genetic data, should be genepop format}
  \item{ncode}{The coding type of your genetic data.}
  }
\details{Function generates two genetic differentiation (deta D) plots, a box plot and a scatter plot}
\value{}
\references{
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
}
\author{qinxinghu@gmail.com}
\note{}
\seealso{DetaDcorplot}
\examples{
##plot genetic differentiation across loci
plotdiff1("Island.gen",ncode=3)}
\keyword{ "plotdiff1"}
