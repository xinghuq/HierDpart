
###### the function for correlation of pirwise DeltaD and distance


COR_DeltaDd = function(f, d, ncode) {
    diveRsity::readGenepop
    gp = ncode
    fr = readGenepop(f, gp, bootstrap = FALSE)
    af = fr$allele_freq
    ade4::mantel.randtest
    DeltaD = function(abun, struc) {
      ## Chao et al, 2017
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
        Ai = sapply(1:N, function(k) -sum(pi[, k][pi[, k] > 0] * log(pi[, k][pi[, k] > 0])))
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
                  } else {
                    ai[, j] = rowSums(abun[, II])
                  }
                }
                pi = sapply(1:NN, function(k) ai[, k]/sum(ai[, k]))
                wi = colSums(ai)/sum(ai)
                W[i - 1] = -sum(wi * log(wi))
                Ai = sapply(1:NN, function(k) -sum(pi[, k][pi[, k] > 0] * log(pi[, k][pi[, k] > 0])))
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
                Dmat[[l]][i, j] = DeltaD((af[[l]][, c(i, j)]), str)[2]  ### select two pops from allelefrequency
            }
        }
    }
    pairwiseDav = Reduce("+", Dmat)/length(Dmat)
    colnames(pairwiseDav) = fr$pop_names
    rownames(pairwiseDav) = fr$pop_names
    # library(popbio)
    DeltaDmat = as.dist(pairwiseDav, diag = FALSE, upper = FALSE)

    if (class(d) == "matrix" | class(d) == "dist") {
      if (length(DeltaDmat) != length(d))
        stop("Numbers of rows in DeltaD matrix and Dgeo are not equal")
      if (sum(is.na(DeltaDmat)) != 0 | sum(is.na(d)) != 0)
        stop("Missing data in the dataset")
            Dgeo = as.dist(d, diag = FALSE, upper = FALSE)
            #COR_DeltaDd = cor(DeltaDmat, Dgeo, method = "pearson")
            COR_DeltaDd=mantel.randtest(DeltaDmat, Dgeo, nrepet = 999)
        }

    if (is.null(d)==TRUE) {
        M = matrix(data = 0, nrow = npops, ncol = npops)
        colnames(M) = fr$pop_names
        rownames(M) = fr$pop_names
        for (i in 1:npops) {
            for (j in 1:npops) {
                M[i, j] = abs(i - j)
            }
        }
        Dgeo = as.dist(M, diag = FALSE, upper = FALSE)
        #ade4::mantel.randtest
        #COR_Fstd = cor.test(PFst, Dgeo, type = "pearson")
        COR_DeltaDd = mantel.randtest(DeltaDmat, Dgeo, nrepet = 999)
    }
    else {if (class(d) != "matrix" & class(d) != "dist")
      stop("d (Dgeo) has to be a matrix")}
    return(list(PairwiseDeltaD = DeltaDmat, Dgeo = Dgeo, CorDeltaDd = COR_DeltaDd))
}

