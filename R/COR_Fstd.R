<<<<<<< HEAD

###### updated on 21-09-2018, the function for correlation of pirwise Fst and distance
=======
###### updated on 09-08-2018, the function for correlation of pirwise Fst and distance
>>>>>>> b0c8cb24be41d276ddfae5df6d537a8eb64f7ed6

COR_Fstd = function (x, d, ncode) {
  read.genepop1 <- function(file, ncode, quiet = FALSE) {
    adegenet::.readExt
    adegenet::.genlab
    adegenet::df2genind
    adegenet::is.genind
     adegenet::pop
     adegenet::repool
    adegenet::Hs
    adegenet::seppop
    adegenet::popNames
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
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
    npop <- length(pop.idx)
    nocomma <- which(!(1:length(txt)) %in% grep(",", txt))
    splited <- nocomma[which(!nocomma %in% pop.idx)]
    if (length(splited) > 0) {
      for (i in sort(splited, decreasing = TRUE)) {
        txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
      }
      txt <- txt[-splited]
    }
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
    txt[length(txt) + 1] <- "POP"
    nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))) - 1
    pop <- factor(rep(1:npop, nind.bypop))
    txt <- txt[-c(pop.idx, length(txt))]
    temp <- sapply(1:length(txt), function(i) strsplit(txt[i], ","))
    ind.names <- vapply(temp, function(e) e[1], character(1))
    ind.names <- trimws(ind.names)
    vec.genot <- vapply(temp, function(e) e[2], character(1))
    vec.genot <- trimws(vec.genot)
    X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")), ncol = nloc, byrow = TRUE)
    if (any(duplicated(ind.names))) {
      rownames(X) <- .genlab("", nrow(X))
    } else {
      rownames(X) <- ind.names
    }
    colnames(X) <- loc.names
    pop.names.idx <- cumsum(table(pop))
    pop.names <- ind.names[pop.names.idx]
    levels(pop) <- pop.names
    if (!all(unique(nchar(X)) == (ncode * 2)))
      stop(paste("some alleles are not encoded with", ncode, "characters\nCheck 'ncode' argument"))
    res <- df2genind(X = X, pop = as.character(pop), ploidy = 2, ncode = ncode, NA.char = NA.char)
    res@call <- prevcall
    if (!quiet)
      cat("\n...done.\n\n")
    return(res)
  }
  g = read.genepop1(x, ncode, quiet = FALSE)
     # genind file
  hierfstat::pairwise.fst
  PFst = pairwise.fst(g)
    npops = length(levels(g$pop))
    ## if M is matrix, then
    if (d == TRUE) {
        if (is.matrix(d) == TRUE) {
            Dgeo = as.dist(d, diag = FALSE, upper = FALSE)
            COR_Fstd = cor(PFst, Dgeo, method = "pearson")
            return(COR_Fstd)
        } else {
            print("d must be a matrix")
        }
    } else {
        d == FALSE
        M = matrix(data = 0, nrow = npops, ncol = npops)  ## nrow=npops, ncol=npops, which is 16 here
        colnames(M) = levels(g$pop)
        rownames(M) = levels(g$pop)
        for (i in 1:npops) {
            for (j in 1:npops) {
                M[i, j] = abs(i - j)
            }
        }

        Dgeo = as.dist(M, diag = FALSE, upper = FALSE)
        # if we want test IBD by mantel test ibd <- mantel.randtest(Dgen2,Dgeo1) plot(ibd) plot(Dgeo1, Dgen2)
        # abline(lm(Dgen2~Dgeo1), col='red',lty=2) if we want use the euclidean distance rather the real matrix distance
        # Dgeo1=dist(M, method = 'euclidean', diag = FALSE, upper = FALSE, p = 2)### this is euclidean distance
        # COR_Fsteud= cor(Dgen2,Dgeo1,method = 'pearson') #cor(c(matrix1), c(matrix2)).

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
