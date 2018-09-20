\name{COR_Fstd}
\alias{COR_Fstd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate pairwise Fst and correlation between Fst and geographic distance.
%%  ~~function to do ... ~~
}
\description{ This function calculates pairwise genetic differentication (Fst) , as well as the correlation between genetic differentiation (Fst) and geographic distance.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
COR_Fstd(x, d, ncode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{  A genetic data files, here the format is genepop format.
%%     ~~Describe \code{x} here~~
}
  \item{d}{ Pairwise geographic distance for calculating correlations. If no actual geographic distance, d=FALSE. The geographic distance should be a matrix that corresponds to the pairwise genetic differentiation matrix.
%%     ~~Describe \code{d} here~~
}
  \item{ncode}{ The code type of the genotype
%%     ~~Describe \code{ncode} here~~
}
}
\details{ This function returns lists of pairwise Fst matrix, and a matrix of theorical geographical distance if d is FALSE, as well as their correlation coefficient.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item {pwFst} { The pairwise genetic differentiation-Fst matrix}
\item {Dgeo} { The pairwise geographical distance matrix}
\item {COR_Fstd} {Person correlation coefficient of genetic differentiation (Fst) and geographic distance}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical F‐statistics. Molecular Ecology Notes, 5(1), 184-186.
Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), 1403-1405.
Weir, B. S.; Cockerham, C. Clark (1984). "Estimating F-Statistics for the Analysis of Population Structure". Evolution. 38 (6): 1358.
Evolution and the Genetics of Populations Volume 2: the Theory of Gene Frequencies, pg 294–295, S. Wright, Univ. of Chicago Press, Chicago, 1969.
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
 Island<- "https://www.github.com/xinghuq/HierDpart/data/Island.gen"
 d2=COR_Fstd("Island.gen",d=FALSE,ncode=3)
 print(d2)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, d, ncode)
{
    read.genepop <- function(file, ncode, quiet = FALSE) {
        require(adegenet)
        if (toupper(.readExt(file)) != "GEN")
            stop("File extension .gen expected")
        if (!quiet)
            cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")
        prevcall <- match.call()
        txt <- scan(file, sep = "\n", what = "character", quiet = TRUE)
        if (!quiet)
            cat("\nFile description: ", txt[1], "\n")
        txt <- txt[-1]
        txt <- gsub("\t", " ", txt)
        NA.char <- paste(rep("0", ncode), collapse = "")
        locinfo.idx <- 1:(min(grep("POP", toupper(txt))) - 1)
        locinfo <- txt[locinfo.idx]
        locinfo <- paste(locinfo, collapse = ",")
        loc.names <- unlist(strsplit(locinfo, "([,]|[\n])+"))
        loc.names <- trimws(loc.names)
        nloc <- length(loc.names)
        txt <- txt[-locinfo.idx]
        pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",
            toupper(txt))
        npop <- length(pop.idx)
        nocomma <- which(!(1:length(txt)) \%in\% grep(",", txt))
        splited <- nocomma[which(!nocomma \%in\% pop.idx)]
        if (length(splited) > 0) {
            for (i in sort(splited, decreasing = TRUE)) {
                txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
            }
            txt <- txt[-splited]
        }
        pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",
            toupper(txt))
        txt[length(txt) + 1] <- "POP"
        nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$",
            toupper(txt))) - 1
        pop <- factor(rep(1:npop, nind.bypop))
        txt <- txt[-c(pop.idx, length(txt))]
        temp <- sapply(1:length(txt), function(i) strsplit(txt[i],
            ","))
        ind.names <- vapply(temp, function(e) e[1], character(1))
        ind.names <- trimws(ind.names)
        vec.genot <- vapply(temp, function(e) e[2], character(1))
        vec.genot <- trimws(vec.genot)
        X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")),
            ncol = nloc, byrow = TRUE)
        if (any(duplicated(ind.names))) {
            rownames(X) <- .genlab("", nrow(X))
        }
        else {
            rownames(X) <- ind.names
        }
        colnames(X) <- loc.names
        pop.names.idx <- cumsum(table(pop))
        pop.names <- ind.names[pop.names.idx]
        levels(pop) <- pop.names
        if (!all(unique(nchar(X)) == (ncode * 2)))
            stop(paste("some alleles are not encoded with", ncode,
                "characters\nCheck 'ncode' argument"))
        res <- df2genind(X = X, pop = as.character(pop), ploidy = 2,
            ncode = ncode, NA.char = NA.char)
        res@call <- prevcall
        if (!quiet)
            cat("\n...done.\n\n")
        return(res)
    }
    g = read.genepop(x, ncode, quiet = FALSE)
    pairwise.fst <- function(x, pop = NULL, res.type = c("dist",
        "matrix")) {
        if (!is.genind(x))
            stop("x is not a valid genind object")
        if (!is.null(pop)) {
            pop(x) <- pop
        }
        temp <- pop(x)
        if (is.null(temp))
            stop("no grouping factor (pop) provided")
        if (length(levels(temp)) < 2) {
            warning("There is only one pop - returning NULL")
            return(NULL)
        }
        res.type <- match.arg(res.type)
        f1 <- function(pop1, pop2) {
            n1 <- nrow(pop1@tab)
            n2 <- nrow(pop2@tab)
            temp <- repool(pop1, pop2)
            b <- weighted.mean(Hs(temp), c(n1, n2))
            pop(temp) <- NULL
            a <- Hs(temp)
            return((a - b)/a)
        }
        lx <- seppop(x, treatOther = FALSE)
        temp <- pop(x)
        levPop <- levels(temp)
        allPairs <- combn(1:length(levPop), 2)
        if (!is.matrix(allPairs)) {
            allPairs <- matrix(allPairs, nrow = 2)
        }
        vecRes <- numeric()
        for (i in 1:ncol(allPairs)) {
            vecRes[i] <- f1(lx[[allPairs[1, i]]], lx[[allPairs[2,
                i]]])
        }
        squelres <- dist(1:length(levPop))
        res <- vecRes
        attributes(res) <- attributes(squelres)
        if (res.type == "matrix") {
            res <- as.matrix(res)
            lab <- popNames(x)
            colnames(res) <- rownames(res) <- lab
        }
        return(res)
    }
    PFst = pairwise.fst(g)
    npops = length(levels(g$pop))
    if (d == TRUE) {
        if (is.matrix(d) == TRUE) {
            Dgeo = as.dist(d, diag = FALSE, upper = FALSE)
            COR_Fstd = cor(PFst, Dgeo, method = "pearson")
            return(COR_Fstd)
        }
        else {
            print("d must be a matrix")
        }
    }
    else {
        d == FALSE
        M = matrix(data = 0, nrow = npops, ncol = npops)
        colnames(M) = levels(g$pop)
        rownames(M) = levels(g$pop)
        for (i in 1:npops) {
            for (j in 1:npops) {
                M[i, j] = abs(i - j)
            }
        }
        Dgeo = as.dist(M, diag = FALSE, upper = FALSE)
        if (class(PFst) != "matrix" & class(PFst) != "dist")
            stop("PFst has to be a matrix")
        if (class(Dgeo) != "matrix" & class(Dgeo) != "dist")
            stop("Dgeo has to be a matrix")
        if (sum(is.na(PFst)) != 0 | sum(is.na(Dgeo)) != 0)
            stop("Missing data in the dataset")
        if (length(PFst) != length(Dgeo))
            stop("Numbers of rows in PFst and Dgeo are not equal")
        COR_Fstd = cor(PFst, Dgeo, method = "pearson")
        return(list(pwFst = PFst, Dgeo = Dgeo, COR_Fstd = COR_Fstd))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
