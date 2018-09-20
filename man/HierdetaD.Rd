\name{HierdetaD}
\alias{HierdetaD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate hierarchial genetic diversity and differentiation (deta D)
%%  ~~function to do ... ~~
}
\description{ This function calculates and decomposes hierarchial genetic diversity and differentiation (deta D) (Oscar,E et al, 2018) into α, β, γ diversity and deta D2 and deta D1 (genetic differentiation).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HierdetaD(x, nreg, r, ncode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The genetype data, genepop format
%%     ~~Describe \code{x} here~~
}
  \item{nreg}{ Number of aggregate in your hierarchical structure.
%%     ~~Describe \code{nreg} here~~
}
  \item{r}{ The number of subaggregate in your hierarchical structure.
%%     ~~Describe \code{r} here~~
}
  \item{ncode}{ Coding type of genetic data.
%%     ~~Describe \code{ncode} here~~
}
}
\details{ This function calculates the hierarchial genetic differentiation (deta D) for metapopulations structured in any hierarchy. In this function, you can input your genepop file (x) and structure (nreg, r) of your data in the arguments, and then returns a list of your hierarchial genetic diversity and differentiation (deta D ). The result contains D_gamma, D_alpha.2, D_alpha.1, D_beta.2, D_beta.1, Differentiation.2, Differentiation.1, which corresponds to different hierarchial levels.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Dst}{Result returns to a list of hierarchial diversities (q=1) and genetic differentiations, including D_gamma (total diversity), D_alpha.2 (average total diversity among aggregate (eg. region) within the whole ecosystem), D_alpha.1 (average diversity among populations), D_beta.2 (number of aggregate equivalent), D_beta.1 (number of population equivalent), Differentiation.2 (genetic differentiation between aggragates), Differentiation.1 (genetic differentiation between populations)}

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
HeirD=HierdetaD("Island.gen",nreg=4,r=c(7,4,2,3),ncode=3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, nreg, r, ncode)
{
    require(diveRsity)
    gp = ncode
    file = readGenepop(x, gp, bootstrap = FALSE)
    outfile = file$allele_freq
    npops = file$npops
    if (sum(r) != npops)
        stop("number of pops should be identical with your defined in region")
    Str = function(nreg, r, npops) {
        str = data.frame(matrix(data = 0, ncol = npops, nrow = 3))
        str[1, ] = c(rep("ecosystem", times = npops))
        if (length(r) != nreg)
            stop("number of regions should be equal to the number of  r  defined")
        str[2, ] = c(rep(paste0("region", 1:nreg), r))
        str[3, ] = c(paste("pop", 1:npops))
        str = as.matrix(str)
        return(str)
    }
    str = Str(nreg, r, npops)
    Hier_detaD = function(abun, struc) {
        n = sum(abun)
        N = ncol(abun)
        ga = rowSums(abun)
        gp = ga[ga > 0]/n
        G = sum(-gp * log(gp))
        H = nrow(struc)
        A = numeric(H - 1)
        W = numeric(H - 1)
        B = numeric(H - 1)
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
        B[1] = exp(G)/exp(A[1])
        if (H > 2) {
            for (i in 2:(H - 1)) {
                Diff[i] = (A[i - 1] - A[i])/(W[i] - W[i - 1])
                B[i] = exp(A[i - 1])/exp(A[i])
            }
        }
        Gamma = exp(G)
        Alpha = exp(A)
        Diff = Diff
        out = matrix(c(Gamma, Alpha, B, Diff), ncol = 1)
        rownames(out) <- c(paste0("D_gamma"), paste0("D_alpha.",
            (H - 1):1), paste0("D_beta.", (H - 1):1), paste0("Differentiation.",
            (H - 1):1))
        return(out)
    }
    Dst = data.frame(matrix(data = 0, ncol = 7, nrow = file$nloci))
    for (i in seq_along(outfile)) {
        Dst[i, ] = t(Hier_detaD(outfile[[i]], str))
    }
    rownames(Dst) = c(paste("Locus", 1:file$nloci))
    colnames(Dst) = c("D_gamma", "D_alpha.2", "D_alpha.1", "D_beta.2",
        "D_beta.1", "Differentiation.2", "Differentiation.1")
    return(Dst)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
