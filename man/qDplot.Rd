\name{qDplot}
\alias{qDplot}
\title{ Plot the genetic diversity profiles (q=0,1,2).
}
\description{ This function calculates the diversity profiles and generates a line plot. You can get either of the diversities, richness (q=0), exponential Shannon entropy (q=1), and heterozygosity (q=2), or all of them.
}
\usage{
qDplot(x, q, ncode)
}
\arguments{
  \item{x}{ The input genetic files, genepop format
}
  \item{q}{ Hill numbers, here you can plot q= "0", "1", "2", "all".If q= "all", function will gives you a line plot including three diversities.

}
  \item{ncode}{ The coding type of your data.

}
}
\details{

}
\value{ A diversity plot.

}
\references{
Marcon, E., & Hérault, B. (2015). entropart: An R package to measure and partition diversity. Journal of Statistical Software, 67(8).
Gaggiotti, O. E., Chao, A., Peres‐Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.

}
\author{

}
\note{

}



\seealso{

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, q, ncode)
{
    D00 = qD(x, q = 0, ncode)
    D11 = qD(x, q = 1, ncode)
    D22 = qD(x, q = 2, ncode)
    D0 = colMeans(D00)
    D1 = colMeans(D11)
    D2 = colMeans(D22)
    require(ggplot2)
    if (q == 0) {
        data = melt(D0)
        data = data.frame(popname = rownames(data), D0)
        p = ggplot(data, aes(x = popname, y = D0, group = 1)) +
            geom_line(size = 1) + ylim(0, (max(D0) + 5)) + ylab(paste("Diversity q=",
            q)) + xlab("Populations") + geom_point(size = 1) +
            theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }
    else if (q == 1) {
        data = melt(D1)
        data = data.frame(popname = rownames(data), D1)
        p = ggplot(data, aes(x = popname, y = D1, group = 1)) +
            geom_line(size = 1) + ylim(0, (max(D1) + 5)) + ylab(paste("Diversity q=",
            q)) + xlab("Populations") + geom_point(size = 1) +
            theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }
    else if (q == 2) {
        data = melt(D2)
        data = data.frame(popname = rownames(data), D2)
        p = ggplot(data, aes(x = popname, y = D2, group = 1)) +
            geom_line(size = 1) + ylim(0, (max(D2) + 5)) + ylab(paste("Diversity q=",
            q)) + xlab("Populations") + geom_point(size = 1) +
            theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }
    else if (q == "all") {
        qtotal = cbind(D0, D1, D2)
        rownames(qtotal) = gsub(",", "", rownames(qtotal))
        require(reshape2)
        data = melt(qtotal)
        colnames(data) = c("Populations", "qD", "Diversity")
        p = ggplot(data, aes(x = Populations, y = Diversity,
            group = qD, colour = qD)) + ylab("Diversity") + xlab("Populations") +
            geom_line(size = 1) + ylim(0, (max(data$Diversity) +
            5)) + geom_point(aes(shape = qD), size = 2) + theme(axis.text.x = element_text(angle = 30,
            hjust = 1))
    }
    else {
        stop("none of the q is workable")
    }
    return(p)
  }
}

\keyword{"qDplot" }
\keyword{"qD"}
