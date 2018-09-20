\name{qDplot}
\alias{qDplot}
\title{Plot the genetic diversity profiles (q=0,1,2).}
\description{This function calculates the diversity profiles and generates a line plot. You can get either one of the diversity index plot, richness (q=0), exponential Shannon entropy (q=1), and heterozygosity (q=2), or all of them.}
\usage{qDplot(x, q, ncode)}
\arguments{
  \item{x}{ The input genetic files, genepop format}
  \item{q}{ Hill numbers,  q can be o (richness), 1 (exponential Shannon entropy), 2 (heterozygosity), or "all".If q= "all", this will returns a line plot including three diversities.}
  \item{ncode}{ The coding type of your data.}
}
\details{Function returns a diversity plot.}
\value{}
\references{
1. Marcon, E., & Hérault, B. (2015). entropart: An R package to measure and partition diversity. Journal of Statistical Software, 67(8).

2. Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{qinxinghu@gmail.com}
\note{}
\seealso{qD}
\examples{# plot allelic richness
dDplot("Island.gen",q=0,ncode=3)
#plot all three diversity profiles
dDplot("Island.gen",q="all",ncode=3)}
\keyword{"qDplot"}
\keyword{"qD"}
