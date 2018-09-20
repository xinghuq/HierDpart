#### Modified HierFst

HierFst = function(x, nreg, r, ncode) {
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

    genfiles = read.genepop(x, ncode, quiet = TRUE)  # covert the genepop #files to genind files, we can also use read.genpop from adegent package
    require(hierfstat)
    hfiles <- genind2hierfstat(genfiles)  # convert into hieformat
    require(dplyr)
    npops = length(levels(genfiles$pop))
   # nloci = length(levels(genfiles$loc.fac))
    sampsize = summary(genfiles$pop)  ## sample size
    if (length(r) != nreg)
        stop("Number of regions should be equal to the number defined in the level")  ## number of pops per region
    if (sum(r) != npops)
        stop("Number of pops should be equal to the number defined in level")

    popr = list()
    rsample = list()
    for (i in 1:nreg) {
        popr[[i]] = list()
        popr[[i]] = as.factor(rep(paste("", i), times = sum(sampsize[(sum(head(r, i - 1)) + 1):(sum(head(r, i)))])))  ### be aware that times depend on the sample size and str on your data
        rsample[[i]] = sum(sampsize[(sum(head(r, i - 1)) + 1):(sum(head(r, i)))])
    }
    region = as.numeric(unlist(popr))
    require(tibble)
    hfiles1 = add_column(hfiles, region, .before = 1)
    hfiles1[, 2] = as.numeric(hfiles1[, 2])
    HieFst = varcomp.glob(data.frame(hfiles1[, 1:2]), hfiles1[, -c(1, 2)])$F

    return(HieFst)
}
