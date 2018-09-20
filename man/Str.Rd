\name{Str}
\alias{Str}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ function to generate the structure for HierDpart format
%%  ~~function to do ... ~~
}
\description{ This function allows you to easily and automatically generate the hierarchical structure for your data. The structure files are HierDpart format working for function IDIP which is easy to construct by defining your aggregates and total number of subaggregate.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Str(nreg, r, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nreg}{ The number of aggregates in your hierarchy.
%%     ~~Describe \code{nreg} here~~
}
  \item{r}{ A vector specifying number of populations or subpopulations in each aggregate.
%%     ~~Describe \code{r} here~~
}
  \item{n}{ An integer indicating the number of subaggregate.
%%     ~~Describe \code{n} here~~
}
}
\details{ This function is designed to fit the input format for function IDIP (Chao et al, 2017). It is more easy to understand and set the hierarchial structure using this fucntion.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ A matrix that displays your hierarchical structure
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
Chao, A., & Chiu, C. H. (2017).User’s Guide for Online program iDIP (Information-based Diversity Partitioning).
%% ~put references to the literature/web site here ~
}
\author{ qinxinghu@gmail.com
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
