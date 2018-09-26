### Hierarchical Diversity decompostion for q =1

## hierarchical diversity over loci

HierD = function(x, nreg, r, ncode) {
    # r is level, pops per region
    diveRsity::readGenepop
    gp = ncode
    file = readGenepop(x, gp, bootstrap = FALSE)
    outfile = file$allele_freq
    npops = file$npops
    if (sum(r) != npops)
        stop("number of pops should be identical with your defined in region")
    Str = function(nreg, r, npops) {
        str = data.frame(matrix(data = 0, ncol = npops, nrow = 3))  #
        str[1, ] = c(rep("ecosystem", times = npops))
        if (length(r) != nreg)
            stop("number of regions should be equal to the number of  r  defined")
        str[2, ] = c(rep(paste0("region", 1:nreg), r))
        str[3, ] = c(paste("pop", 1:npops))
        str = as.matrix(str)
        return(str)
    }
    str = Str(nreg, r, npops)
    Hier_detaD = function(abun, struc) {
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
        B[1] = exp(G)/exp(A[1])
        if (H > 2) {
            for (i in 2:(H - 1)) {
                Diff[i] = (A[i - 1] - A[i])/(W[i] - W[i - 1])
                B[i] = exp(A[i - 1])/exp(A[i])
            }
        }
        Gamma = exp(G)
        Alpha = exp(A)
        Diff = Diff
        out = matrix(c(Gamma, Alpha, B, Diff), ncol = 1)
        rownames(out) <- c(paste0("D_gamma"), paste0("D_alpha.", (H - 1):1), paste0("D_beta.", (H - 1):1), paste0("Differentiation.",
            (H - 1):1))
        return(out)
    }

    Dst = data.frame(matrix(data = 0, ncol = 7, nrow = file$nloci))

    for (i in seq_along(outfile)) {
        Dst[i, ] = t(Hier_detaD(outfile[[i]], str))
    }
    rownames(Dst) = c(paste("Locus", 1:file$nloci))
    colnames(Dst) = c("D_gamma", "D_alpha.2", "D_alpha.1", "D_beta.2", "D_beta.1", "Differentiation.2", "Differentiation.1")
    return(Dst)
}
