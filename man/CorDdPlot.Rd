\name{CorDdPlot}
\alias{CorDdPlot}
\title{ Plot genetic differentiation (deta D) and geographic distance.
}
\description{ This function plots genetic differentiation (deta D) and geographic distance.
}
\usage{
CorDdPlot(x, d, ncode)
}
\arguments{
  \item{x}{ The genetic file, here is genepop format file.
}
  \item{d}{ Geographic distance.
}
  \item{ncode}{ The type of your coding allele
}
}
\details{
}
\value{
\item{plot}{ Return to a linear plot for the relation between genetic differentiation (deta D) and geographic distance.}
\item{lm}{ Return to details of linear regression between genetic differentiation (deta D) and geographic distance.}
}
\references{
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
}
\author{qinxinghu@gmail.com
}
\note{
}
\seealso{
}
\examples{
CorDdPlot("Island.gen",d=FALSE,ncode=3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, d, ncode)
{
    data = COR_detaDd(x, d, ncode)
    require(ggplot2)
    a = plot(data$Dgeo, data$PairwiseDetaD, xlab = "Geographic Distance",
        ylab = "Genetic differentiation (detaD)")
    abline(lm(data$PairwiseDetaD ~ data$Dgeo))
    lm = lm(data$PairwiseDetaD ~ data$Dgeo)
    return(summary(lm))
  }
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
