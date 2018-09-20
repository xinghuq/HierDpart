\name{IDIP}
\alias{IDIP}
\title{ Function to decompose diversity (genetic or species diversity) into
}
\description{ This function is from Information-based Diversity Partitioning (Chao et al, 2017). This function allows you to do analyses of decomposing information-based diversity under a specified multi-level hierarchical structure.
}
\usage{
IDIP(abun, struc)
}
\arguments{
  \item{abun}{ The count or frequency data, raw or relative species/allele abundances.
}
  \item{struc}{ The hierarchial structure.
}
}
\details{
}
\value{ Program iDIP can be applied to any arbitrary number of hierarchical levels. The output consists of a basic data summary and decomposition results, with the latter including (1) gamma (or total) diversity; alpha and beta diversity at each level; (2) proportion of total beta information (Shannon information) found at each level; (3) mean differentiation (dissimilarity) among aggregates at each level (Chao et al, 2017).
}
\references{
Chao, A., & Chiu, C. H. User’s Guide for Online program iDIP (Information-based Diversity Partitioning).
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
}
\author{

}
\note{
}


\seealso{

}
\examples{
str=Str(nreg=4,r=c(7,4,2,3),n=16)
abu=matrix(data=runif(16*3,min=0,max=1),nrow = 20,ncol = 16)
IDIP(abu,str)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (abun, struc)
{
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
    Prop = numeric(H - 1)
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
            pi = sapply(1:NN, function(k) ai[, k]/sum(ai[, k]))
            wi = colSums(ai)/sum(ai)
            W[i - 1] = -sum(wi * log(wi))
            Ai = sapply(1:NN, function(k) -sum(pi[, k][pi[, k] >
                0] * log(pi[, k][pi[, k] > 0])))
            A[i - 1] = sum(wi * Ai)
        }
    }
    total = G - A[H - 1]
    Diff[1] = (G - A[1])/W[1]
    Prop[1] = (G - A[1])/total
    B[1] = exp(G)/exp(A[1])
    if (H > 2) {
        for (i in 2:(H - 1)) {
            Diff[i] = (A[i - 1] - A[i])/(W[i] - W[i - 1])
            Prop[i] = (A[i - 1] - A[i])/total
            B[i] = exp(A[i - 1])/exp(A[i])
        }
    }
    Gamma = exp(G)
    Alpha = exp(A)
    Diff = Diff
    Prop = Prop
    out = matrix(c(Gamma, Alpha, B, Prop, Diff), ncol = 1)
    rownames(out) <- c(paste0("D_gamma"), paste0("D_alpha.",
        (H - 1):1), paste0("D_beta.", (H - 1):1), paste0("Proportion.",
        (H - 1):1), paste0("Differentiation.", (H - 1):1))
    return(out)
  }
}

\keyword{ "IDIP" }

