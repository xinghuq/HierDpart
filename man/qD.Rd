\name{qD}
\alias{qD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate genetic diversity profiles (q=0,1,2)
%%  ~~function to do ... ~~
}
\description{ This function calculates genetic diversity profiles (q=0,1,2),including richness (q=0), exponential Shannon entropy (q=1), and heterozygosity (q=2).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qD(x, q, ncode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Input data, genepop format
%%     ~~Describe \code{x} here~~
}
  \item{q}{ Hill numbers, q=0,1,2.
%%     ~~Describe \code{q} here~~
}
  \item{ncode}{ The coding type of the genetic file
%%     ~~Describe \code{ncode} here~~
}
}
\details{ This function returns different genetic diversities, q is Hill numbers, specifically, richness (q=0), exponential Shannon entropy (q=1), and heterozygosity (q=2).
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Result returns the diversity profiles (q) as your want.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Marcon, E., & Hérault, B. (2015). entropart: An R package to measure and partition diversity. Journal of Statistical Software, 67(8).
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
##plot allelic richness
qD("Island.gen",q=0,ncode=3)
#plot three diversity profiles
qD("Island.gen",q="all",ncode=3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, q, ncode)
{
    require(diveRsity)
    gp = ncode
    file = readGenepop(x, gp, bootstrap = FALSE)
    outfile = file$allele_freq
    npops = file$npops
    nloci = file$nloci
    D = as.data.frame(matrix(data = 0, ncol = npops, nrow = nloci))
    require(entropart)
    for (i in 1:nloci) {
        for (j in 1:npops) {
            D[i, j] = Diversity(outfile[[i]][, j], q)
        }
    }
    rownames(D) = file$loci_names
    colnames(D) = file$pop_names
    return(D)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ qD }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
