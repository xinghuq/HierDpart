
CorDdPlot = function(x, d, ncode) {
    data = COR_detaDd(x, d, ncode)
    require(ggplot2)
     a= plot(data$Dgeo, data$PairwiseDetaD, xlab = "Geographic Distance", ylab = "Genetic differentiation (detaD)")
    abline(lm(data$PairwiseDetaD ~ data$Dgeo))
    lm = lm(data$PairwiseDetaD ~ data$Dgeo)
    return(summary(lm))
}
