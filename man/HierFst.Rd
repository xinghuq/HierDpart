\name{HierFst}
\alias{HierFst}

\title{ Function to calculate hierarchial genetic differentiation (Fst, Weir B.S., 1996;Yang R.C.,1998).

}
\description{The function calculate hierarchial genetic differentiation (Fst, Weir B.S., 1996; Yang R.C., 1998). The hierarchy you set in the input arguments (nreg, r) will automatically generate your hierarchical structure. Based on the parameters you set, the result returns to a F-statistics matrix displaying low level (rows) of Fst within a higher level hierarchy (columns).

}
\usage{
HierFst(x, nreg, r, ncode)
}

\arguments{
  \item{x}{ The name of the data file, or the string giving the path to the file. The file should be a genepop object, with with the appropriate extension.

}
  \item{nreg}{ The number of aggregates in your hierarchy.

}
  \item{r}{ A vector specifying number of populations or subpopulations in each aggregate.

}
  \item{ncode}{ an integer indicating the number of characters used to code an allele.

}
}
\details{ This function was modified and improved from the varcomp.glob function from the hierfstat package (Jerome Goudet & Thibaut Jombart 2015). It was optimized to easily set the hierarchial structure while varcom.glob requires the structure data should be already correctly formated in accordance with hierfstat data. In this function, you can input your genepop file (x) and structure (nreg, r) of your data in the arguments, and then returns a matrix of your hierarchial Fst.

}
\value{ The result returns a matrix of hierarchical F-statistics type-coefficients at different levels.

}
\references{ Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.
Yang, R.C. (1998). Estimating hierarchical F-statistics. Evolution 52(4):950-956
Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and Fstatistics.
Molecular Ecology Notes. 5:184-186

}
\author{ qinxinghu@gmail.com

}
\note{ Note that you should arrange your population in the order corresponding to your structure set in nreg.

}



\seealso{

}
\examples{ # load data
             data(Island.gen)
             hFst=HierFst("Island.gen"",nreg=4,r=c(7,4,2,3),ncode=3)
             print(hFst)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, nreg, r, ncode)
{
require(adegenet)   
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
    require(dplyr)
    npops = length(levels(genfiles$pop))
    sampsize = summary(genfiles$pop)
    if (length(r) != nreg)
        stop("Number of regions should be equal to the number defined in the level")
    if (sum(r) != npops)
        stop("Number of pops should be equal to the number defined in level")
    popr = list()
    rsample = list()
    for (i in 1:nreg) {
        popr[[i]] = list()
        popr[[i]] = as.factor(rep(paste("", i), times = sum(sampsize[(sum(head(r,
            i - 1)) + 1):(sum(head(r, i)))])))
        rsample[[i]] = sum(sampsize[(sum(head(r, i - 1)) + 1):(sum(head(r,
            i)))])
    }
    region = as.numeric(unlist(popr))
    require(tibble)
    hfiles1 = add_column(hfiles, region, .before = 1)
    hfiles1[, 2] = as.numeric(hfiles1[, 2])
    HieFst = varcomp.glob(data.frame(hfiles1[, 1:2]), hfiles1[,
        -c(1, 2)])$F
    return(HieFst)
  }
}

\keyword{ HierFst}
\keyword{ Fst}
