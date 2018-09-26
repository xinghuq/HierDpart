\name{Str}
\alias{Str}
\title{Function to generate the structure for HierDpart format}
\description{This function allows you to easily and automatically generate the hierarchical structure for your data. The structure files are HierDpart format working for function IDIP which is easy to construct by defining your aggregates and total number of subaggregate.}
\usage{Str(nreg, r, n)}
\arguments{
  \item{nreg}{The number of aggregates in your hierarchy.}
  \item{r}{A vector specifying number of populations or subpopulations in each aggregate.}
  \item{n}{An integer indicating the number of subaggregate.}
}
\details{This function is designed to fit the input format for function IDIP (Chao et al, 2017). It is more easy to understand and set the hierarchial structure using this fucntion.}
\value{A matrix that displays your hierarchical structure.}
\references{
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
Chao, A., & Chiu, C. H. (2017).User's Guide for Online program iDIP (Information-based Diversity Partitioning).}
\author{ qinxinghu@gmail.com}
\seealso{IDIP}
\examples{ #assume we are working on a metacommunity contains 16 communities,
   #there are 4 aggregates in this metacommunity, and each aggregates
   #contains 7, 4,3,2 subaggregates or communities.The structure will be:
 str1=Str(nreg=4, r=c(7,4,3,2),n=16)
 }
\keyword{Structure}
