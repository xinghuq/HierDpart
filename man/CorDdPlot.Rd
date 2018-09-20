\name{CorDdPlot}
\alias{CorDdPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot genetic differentiation (deta D) and geographic distance.
%%  ~~function to do ... ~~
}
\description{ This function plots genetic differentiation (deta D) and geographic distance.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CorDdPlot(x, d, ncode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The genetic file, here is genepop format file.
%%     ~~Describe \code{x} here~~
}
  \item{d}{ Geographic distance.
%%     ~~Describe \code{d} here~~
}
  \item{ncode}{ The type of your coding allele
%%     ~~Describe \code{ncode} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{plot}{ Return to a linear plot for the relation between genetic differentiation (deta D) and geographic distance.}
\item{lm}{ Return to details of linear regression between genetic differentiation (deta D) and geographic distance.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
