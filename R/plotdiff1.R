#### Pairwise detaD across loci

plotdiff1 = function(x, ncode) {
    diveRsity::readGenepop
    gp = ncode
    fr = readGenepop(x, gp, bootstrap = FALSE)
    af = fr$allele_freq
    DetaD = function(abun, struc) {
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
                ## for every loci has one differentiation value
                Dmat[[l]][i, j] = DetaD((af[[l]][, c(i, j)]), str)[2]  ### select two pops from allelefrequency
            }
        }
    }
    pwdetaD = lapply(Dmat, function(x) {
        x[x != 0]
    })
    pwdetaD <- do.call(cbind, lapply(pwdetaD, as.vector))
    pwdetaD = unique(pwdetaD)
    Loci = c(1:fr$nloci)
    reshape2::melt
    pwd <- melt(pwdetaD)
    colnames(pwd) = c("pwpop", "Loci", "detaD")
    ylab = c("Pairwise Genetic Differentiation Across Loci (deta_D)")
    requireNamespace("ggplot2")
    p1 <- ggplot(pwd, aes(x = Loci, y = pwd$detaD, colour = factor(Loci))) + ylab(ylab) + ylim(-0.25, 1) + geom_jitter(size = 4,
        alpha = 0.8) + scale_x_continuous("Loci", labels = as.character(pwd$Loci), breaks = pwd$Loci) + geom_hline(yintercept = mean(pwd$detaD),
        linetype = "dashed", color = "red", size = 1) + geom_smooth(method = loess, se = FALSE, fullrange = TRUE)
    p2 = boxplot(pwd$detaD ~ pwd$Loci, xlab = "Loci", ylab = "Diff1", main = "Pairwise Genetic Differentiation Across Loci (deta_D)")
    return(list(p1 = p1, p2 = p2))
}

