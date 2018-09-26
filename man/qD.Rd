\name{qD}
\alias{qD}
\title{Function to calculate genetic diversity profiles (q=0,1,2)}
\description{This function calculates genetic diversity profiles (q=0,1,2), including richness (q=0), exponential Shannon entropy (q=1), and heterozygosity (q=2).}
\usage{qD(x, q, ncode)}
\arguments{
  \item{x}{Input data, genepop format.}
  \item{q}{Hill numbers, q=0,1,2.}
  \item{ncode}{The coding type of the genetic file.}
}
\details{This function returns different genetic diversities, q is Hill numbers, specifically, richness (q=0), exponential Shannon entropy (q=1), and heterozygosity (q=2).}
\value{Result returns to the diversity profiles (q) as you defined.}
\references{
Marcon, E., & Herault, B. (2015). entropart: An R package to measure and partition diversity. Journal of Statistical Software, 67(8).

Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
}
\author{qinxinghu@gmail.com}
\examples{
##plot allelic richness
download.file("https://sites.google.com/site/advancesinevolutionandgenetics
/Home/unit-3-1/Island.gen?attredirects=0&d=1",
"C:/Users/Public/Documents/Island.gen",cacheOK=FALSE)
qD("C:/Users/Public/Documents/Island.gen",q=0,ncode=3)
}
\keyword{Diversity profiles}

