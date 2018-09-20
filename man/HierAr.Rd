\name{HierAr}
\alias{HierAr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to calculate hierarchial allelic richness
%%  ~~function to do ... ~~
}
\description{ This function calculate hierarchial allelic richness for populations structured in hierarchy.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HierAr(x, nreg, r, ncode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The data file, genepop format
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
\details{ This function calculates the hierarchial allelic richness for metapopulations structured in any hierarchy. In this function, you can input your genepop file (x) and structure (nreg, r) of your data in the arguments, and then returns a matrix of your hierarchial Ar.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Ar_pop}{Allelic richness for each population across loci}
\item{Ar_reg}{Allelic richness for each aggregate ( for example, region)}
\item{Ar_overall}{Hierarchical allelic richness in different levels}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
El Mousadik A. and Petit R.J. (1996) High level of genetic differentiation for allelic richness among populations of the argan tree argania spinosa skeels endemic to Morocco. Theoretical and Applied Genetics, 92:832-839
Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical F‐statistics. Molecular Ecology Notes, 5(1), 184-186.
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
\examples{ # we want to calculate the hierarchial allelic richness of 16 populations structured in a 4 regions. Each region contains 7,4,2,3 populations.
HAr=HierAr=("Island.gen",nreg=4,r=c(7,4,2,3),ncode=3)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, nreg, r, ncode)
{
    read.genepop <- function(file, ncode, quiet = FALSE) {
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
    genfiles = read.genepop(x, ncode, quiet = TRUE)
    require(hierfstat)
    hfiles <- genind2hierfstat(genfiles)
    sampsize = summary(genfiles$pop)
    require(dplyr)
    npops = length(levels(genfiles$pop))
    if (length(r) != nreg)
        stop("Number of regions should be equal to the number defined in the level")
    if (sum(r) != npops)
        stop("Number of pops should be equal to the number defined in level")
    popr = list()
    rsample = list()
    for (i in 1:nreg) {
        popr[[i]] = list()
        popr[[i]] = as.factor(rep(paste("pop", i), times = sum(sampsize[(sum(head(r,
            i - 1)) + 1):(sum(head(r, i)))])))
        rsample[[i]] = sum(sampsize[(sum(head(r, i - 1)) + 1):(sum(head(r,
            i)))])
    }
    popeco = as.factor(rep("ecosystem", times = length(genfiles$pop)))
    rsample = as.data.frame(rsample)
    rsample = as.numeric(unlist(rsample))
    region = list()
    for (i in seq_along(r)) {
        region[[i]] = list()
        region[[i]] = hfiles[(sum(head(rsample, i - 1)) + 1):(sum(head(rsample,
            i))), ]
        region[[i]]$pop = factor(region[[i]]$pop)
    }
    arregion = list()
    hierar = list()
    hierarav = list()
    arecosystem = hfiles
    arecosystem$pop = factor(popeco)
    for (i in seq_along(r)) {
        arregion[[i]] = region[[i]]
        arregion[[i]]$pop = factor(popr[[i]])
        hierar[[i]] = list()
        hierarav[[i]] = list()
        hierar[[i]] = allelic.richness(arregion[[i]], min.n = NULL,
            diploid = TRUE)$Ar
        hierarav[[i]] = colMeans(hierar[[i]])
    }
    hierAr_R = do.call(cbind, lapply(hierarav, data.frame))
    hierAr_Rav = rowMeans(hierAr_R)
    hierarpop = allelic.richness(hfiles, min.n = NULL, diploid = TRUE)$Ar
    hierarpopav = colMeans(hierarpop)
    Arpopav = mean(hierarpopav)
    hierareco = allelic.richness(arecosystem, min.n = NULL, diploid = TRUE)$Ar
    hierarecoav = colMeans(hierareco)
    hierAr = cbind(hierarecoav, hierAr_Rav, Arpopav)
    colnames(hierAr) = c("Art", "Arr", "Arp")
    colnames(hierarpop) = c(paste("Arpop", 1:npops))
    colnames(hierAr_R) = c(paste("Ar_region", 1:nreg))
    return(list(Ar_pop = hierarpop, Ar_reg = hierAr_R, Ar_ovell = hierAr))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
