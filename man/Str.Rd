\name{Str}
\alias{Str}
\title{ function to generate the structure for HierDpart format
}
\description{ This function allows you to easily and automatically generate the hierarchical structure for your data. The structure files are HierDpart format working for function IDIP which is easy to construct by defining your aggregates and total number of subaggregate.
}
\usage{
Str(nreg, r, n)
}
\arguments{
  \item{nreg}{ The number of aggregates in your hierarchy.
}
  \item{r}{ A vector specifying number of populations or subpopulations in each aggregate.

}
  \item{n}{ An integer indicating the number of subaggregate.

}
}
\details{ This function is designed to fit the input format for function IDIP (Chao et al, 2017). It is more easy to understand and set the hierarchial structure using this fucntion.

}
\value{ A matrix that displays your hierarchical structure

}
\references{
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
Chao, A., & Chiu, C. H. (2017).User’s Guide for Online program iDIP (Information-based Diversity Partitioning).

}
\author{ qinxinghu@gmail.com

}
\note{

}


\seealso{

}
\examples{ # assume we are working on a metacommunity contains 16 communities, there are 4 aggregates in this metacommunity, and each aggregates contains 7, 4,3,2 subaggregates or communities.
 str1=Str(nreg=4, r=c(7,4,3,2),n=16)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (nreg, r, n)
{
    Str = data.frame(matrix(data = 0, ncol = n, nrow = 3))
    Str[1, ] = c(rep("ecosystem", times = n))
    if (length(r) != nreg)
        stop("number of regions should be equal to the number of  r  defined")
    Str[2, ] = c(rep(paste0("region", 1:nreg), r))
    Str[3, ] = c(paste("pop", 1:n))
    Str = as.matrix(Str)
    return(Str)
  }
}

\keyword{ "Str" }
\keyword{ "IDIP" }
