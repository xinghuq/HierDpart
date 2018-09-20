\name{plotdiff1}
\alias{plotdiff1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ PLot genetic differentiation (deta D) across loci
%%  ~~function to do ... ~~
}
\description{ Plot genetic differentiation across (deta D, Oscar E, et al, 2018) loci.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotdiff1(x, ncode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The genetic data, should be genepop format
%%     ~~Describe \code{x} here~~
}
  \item{ncode}{ The coding type of your genetic data.
%%     ~~Describe \code{ncode} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
##plot genetic differentiation across loci
plotdiff1("Island.gen",ncode=3)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, ncode)
{
    require(diveRsity)
    gp = ncode
    fr = readGenepop(x, gp, bootstrap = FALSE)
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
    pwdetaD = lapply(Dmat, function(x) {
        x[x != 0]
    })
    pwdetaD <- do.call(cbind, lapply(pwdetaD, as.vector))
    pwdetaD = unique(pwdetaD)
    Loci = c(1:fr$nloci)
    require(reshape2)
    pwd <- melt(pwdetaD)
    colnames(pwd) = c("pwpop", "Loci", "detaD")
    ylab = c("Pairwise Genetic Differentiation Across Loci (deta_D)")
    require(ggplot2)
    p1 <- ggplot(pwd, aes(x = Loci, y = detaD, colour = factor(Loci))) +
        ylab(ylab) + ylim(-0.25, 1) + geom_jitter(size = 4, alpha = 0.8) +
        scale_x_continuous("Loci", labels = as.character(pwd$Loci),
            breaks = pwd$Loci) + geom_hline(yintercept = mean(pwd$detaD),
        linetype = "dashed", color = "red", size = 1) + geom_smooth(method = loess,
        se = FALSE, fullrange = TRUE)
    p2 = boxplot(pwd$detaD ~ pwd$Loci, xlab = "Loci", ylab = "Diff1",
        main = "Pairwise Genetic Differentiation Across Loci (deta_D)")
    return(list(p1 = p1, p2 = p2))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plotdiff1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ DetaDcorplot }% __ONLY ONE__ keyword per line
