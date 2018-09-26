##### Here are scprits for getting hierchical heterozygosity ###################

HierHe = function(x, nreg, r, ncode) {
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

    genfiles = read.genepop1(x, ncode, quiet = TRUE)  # covert the genepop #files to genind files, we can also use read.genpop from adegent package
    hierfstat::genind2hierfstat
    hfiles <- genind2hierfstat(genfiles)  # convert into hieformat

    ## Here we add our hierchical information (regions-pops) to the data
    requireNamespace("dplyr")
    npops = length(levels(genfiles$pop))
    nloci = length(levels(genfiles$loc.fac))
    sampsize = summary(genfiles$pop)  ## sample size
    if (length(r) != nreg)
        stop("Number of regions should be equal to the number defined in the level")  ## number of pops per region
    if (sum(r) != npops)
        stop("Number of pops should be equal to the number defined in level")
    rsample = list()
    for (i in 1:nreg) {
        rsample[[i]] = sum(sampsize[(sum(head(r, i - 1)) + 1):(sum(head(r, i)))])
    }
    rsample = as.data.frame(rsample)
    rsample = as.numeric(unlist(rsample))
    region = list()
    hierHr = list()
    hierHrperloc = list()
    for (i in seq_along(r)) {
        region[[i]] = list()
        region[[i]] = hfiles[(sum(head(rsample, i - 1)) + 1):(sum(head(rsample, i))), ]
        region[[i]]$pop = factor(region[[i]]$pop)
        hierHr[[i]] = list()
        hierHrperloc[[i]] = list()
        hierHr[[i]] = t(basic.stats(region[[i]], diploid = TRUE, digits = 4)$overall)
        hierHrperloc[[i]] = basic.stats(region[[i]], diploid = TRUE, digits = 4)$perloc
    }

    Hrperloc = matrix(data = 0, ncol = nreg, nrow = nloci)
    for (j in seq(hierHrperloc)) {
        Hrperloc[, j] = hierHrperloc[[i]]$Ht
    }

    HierHr_T = t(basic.stats(hfiles, diploid = TRUE, digits = 4)$overall)
    HierHr_Tperloc = basic.stats(hfiles, diploid = TRUE, digits = 4)$perloc
    Hallperloc = cbind(HierHr_Tperloc$Ht, Hrperloc, HierHr_Tperloc$Hs)
    colnames(Hallperloc) = c("Ht", paste("Hr", 1:nreg), "Hp")
    rownames(Hallperloc) = c(paste("Locus", 1:nloci))
    hierHr = do.call(rbind, lapply(hierHr, data.frame))
    Hr = colMeans(hierHr)
    Hr = as.array(Hr)
    HierHr_T = as.data.frame(HierHr_T)
    HieHr = cbind(HierHr_T$Ht, Hr[3], HierHr_T$Hs)
    colnames(HieHr) = c("Ht", "Hr", "Hp")
    rownames(HieHr) = NULL
    return(list(HierHe_perloc = Hallperloc, HieHr_overall = HieHr))
}
