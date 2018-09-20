## The function to get allele diversity (q=1,2,3)

qD = function(x, q, ncode) {
    # q=1,2,3
    require(diveRsity)
    gp = ncode
    file = readGenepop(x, gp, bootstrap = FALSE)
    outfile = file$allele_freq
    npops = file$npops
    nloci = file$nloci
    D = as.data.frame(matrix(data = 0, ncol = npops, nrow = nloci))
    require(entropart)
    for (i in 1:nloci) {
        for (j in 1:npops) {
            D[i, j] = Diversity(outfile[[i]][, j], q)  # n is the number of files, and i is loci, j is pops
        }
    }
    rownames(D) = file$loci_names
    colnames(D) = file$pop_names
    return(D)
}
