\name{COR_detaDd}
\alias{COR_detaDd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate pairwise detaD and correlation between genetic differentiation (deta D) and geographic distance.
%%  ~~function to do ... ~~
}
\description{ This function calculates new pairwise genetic differentication, deta D (Oscar E, et al, 2018), as well as the correlation between genetic differentiation (deta D) and geographic distance.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
COR_detaDd(f, d, ncode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ A genetic data files, here the format is genepop format.
%%     ~~Describe \code{f} here~~
}
  \item{d}{ Pairwise geographic distance for calculating correlations. If no actual geographic distance, d=FALSE. The geographic distance should be a matrix that corresponds to the pairwise genetic differentiation matrix.
%%     ~~Describe \code{d} here~~
}
  \item{ncode}{ The code type of the genotype
%%     ~~Describe \code{ncode} here~~
}
}
\details{ This function returns lists of pairwise detaD matrix, and a matrix of theorical geographical distance if d is FALSE, as well as their correlation coefficient.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item {PairwiseDetaD} { The pairwise genetic differentiation-detaD matrix}
\item {Dgeo} { The pairwise geographical distance matrix}
\item {CordetaDd} {Person correlation coefficient of genetic differentiation (deta D) and geographic distance}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
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
\examples{
 d1=COR_detaDd("Island.gen",d=FALSE,ncode=3)
 print(d1)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (f, d, ncode)
{
    require(diveRsity)
    gp = ncode
    fr = readGenepop(f, gp, bootstrap = FALSE)
    af = fr$allele_freq
    DetaD = function(abun, struc) {
        n = sum(abun)
        N = ncol(abun)
        ga = rowSums(abun)
        gp = ga[ga > 0]/n
        G = sum(-gp * log(gp))
        H = nrow(struc)
        A = numeric(H - 1)
        W = numeric(H - 1)
        Diff = numeric(H - 1)
        wi = colSums(abun)/n
        W[H - 1] = -sum(wi[wi > 0] * log(wi[wi > 0]))
        pi = sapply(1:N, function(k) abun[, k]/sum(abun[, k]))
        Ai = sapply(1:N, function(k) -sum(pi[, k][pi[, k] > 0] *
            log(pi[, k][pi[, k] > 0])))
        A[H - 1] = sum(wi * Ai)
        if (H > 2) {
            for (i in 2:(H - 1)) {
                I = unique(struc[i, ])
                NN = length(I)
                ai = matrix(0, ncol = NN, nrow = nrow(abun))
                c
                for (j in 1:NN) {
                  II = which(struc[i, ] == I[j])
                  if (length(II) == 1) {
                    ai[, j] = abun[, II]
                  }
                  else {
                    ai[, j] = rowSums(abun[, II])
                  }
                }
                pi = sapply(1:NN, function(k) ai[, k]/sum(ai[,
                  k]))
                wi = colSums(ai)/sum(ai)
                W[i - 1] = -sum(wi * log(wi))
                Ai = sapply(1:NN, function(k) -sum(pi[, k][pi[,
                  k] > 0] * log(pi[, k][pi[, k] > 0])))
                A[i - 1] = sum(wi * Ai)
            }
        }
        Diff[1] = (G - A[1])/W[1]
        if (H > 2) {
            for (i in 2:(H - 1)) {
                Diff[i] = (A[i - 1] - A[i])/(W[i] - W[i - 1])
            }
        }
        Diff = Diff
        out = matrix(c(Diff), ncol = 1)
        return(out)
    }
    v1 = c("ecosystem", "region1", "pop1")
    v2 = c("ecosystem", "region1", "pop2")
    str = data.frame(v1, v2)
    str = as.matrix(str)
    npops = fr$npops
    nloci = fr$nloci
    Dmat = list()
    for (l in 1:nloci) {
        Dmat[[l]] = matrix(data = 0, nrow = npops, ncol = npops)
        for (i in 1:npops) {
            for (j in 1:npops) {
                Dmat[[l]][i, j] = DetaD((af[[l]][, c(i, j)]),
                  str)[2]
            }
        }
    }
    pairwiseDav = Reduce("+", Dmat)/length(Dmat)
    detaDmat = as.dist(pairwiseDav, diag = FALSE, upper = FALSE)
    if (d == TRUE) {
        if (is.matrix(d) == TRUE) {
            Dgeo = as.dist(d, diag = FALSE, upper = FALSE)
            COR_detaDd = cor(detaDmat, d, method = "pearson")
            return(COR_detaDd)
        }
        else {
            print("d must be a matrix")
        }
    }
    else {
        d == FALSE
        M = matrix(data = 0, nrow = npops, ncol = npops)
        colnames(M) = fr$pop_names
        rownames(M) = fr$pop_names
        for (i in 1:npops) {
            for (j in 1:npops) {
                M[i, j] = abs(i - j)
            }
        }
        Dgeo = as.dist(M, diag = FALSE, upper = FALSE)
        COR_detaDd = cor(detaDmat, Dgeo, method = "pearson")
        return(list(PairwiseDetaD = detaDmat, Dgeo = Dgeo, CordetaDd = COR_detaDd))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
