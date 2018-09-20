#### Deifine hierarchical structure for genetic data

Str = function(nreg, r, n) {
    Str = data.frame(matrix(data = 0, ncol = n, nrow = 3))  #
    Str[1, ] = c(rep("ecosystem", times = n))
    # should be equal to nreg
    if (length(r) != nreg)
        stop("number of regions should be equal to the number of  r  defined")
    Str[2, ] = c(rep(paste0("region", 1:nreg), r))
    Str[3, ] = c(paste("pop", 1:n))
    Str = as.matrix(Str)
    return(Str)
}
