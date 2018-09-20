\name{qD}
\alias{qD}
\title{ Function to calculate genetic diversity profiles (q=0,1,2)
}
\description{ This function calculates genetic diversity profiles (q=0,1,2),including richness (q=0), exponential Shannon entropy (q=1), and heterozygosity (q=2).
}
\usage{
qD(x, q, ncode)
}
\arguments{
  \item{x}{ Input data, genepop format
}
  \item{q}{ Hill numbers, q=0,1,2.
}
\item{ncode}{ The coding type of the genetic file
}
}
\details{ This function returns different genetic diversities, q is Hill numbers, specifically, richness (q=0), exponential Shannon entropy (q=1), and heterozygosity (q=2).

}
\value{ Result returns the diversity profiles (q) as your want.
}
\references{
Marcon, E., & Hérault, B. (2015). entropart: An R package to measure and partition diversity. Journal of Statistical Software, 67(8).
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
}
\author{qinxinghu@gmail.com
}
\note{
}


\seealso{

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
\keyword{ "qD" }

