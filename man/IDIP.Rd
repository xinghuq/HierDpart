\name{IDIP}
\alias{IDIP}
\title{Function to decompose diversity (genetic or species diversity) and differentiation into different hierarchical levels}
\description{This function comes from Information-based Diversity Partitioning (Chao et al, 2017). It allows you to decompose diversity under a specified multi-level hierarchical structure.
}
\usage{IDIP(abun, struc)}
\arguments{
  \item{abun}{The count or frequency data, raw or relative species/allele abundances.}
  \item{struc}{The hierarchial structure.}
}
\details{This function is definitely an useful tool to do analysis of any information based diversity decomposition. Whether the data is genetic allele/species count/abundance or other frequency data in ecology, chemistry or economy.}
\value{ IDIP can be applied to any arbitrary number of hierarchical levels. The output consists of a basic data summary and decomposition results, with the latter including (1) gamma (or total) diversity; alpha and beta diversity at each level; (2) proportion of total beta information (Shannon information) found at each level; (3) mean differentiation (dissimilarity) among aggregates at each level (Chao et al, 2017).}
\references{
Chao, A., & Chiu, C. H. User's Guide for Online program iDIP (Information-based Diversity Partitioning).

Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
}
\author{chao@stat.nthu.edu.tw}
\examples{
str=Str(nreg=4,r=c(7,4,2,3),n=16)
abu=matrix(data=runif(16*3,min=0,max=1),nrow = 20,ncol = 16)
IDIP(abu,str)
}
\keyword{Hierarchical structure decomposition}

