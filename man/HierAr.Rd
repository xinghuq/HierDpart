\name{HierAr}
\alias{HierAr}
\title{Function to calculate hierarchical allelic richness}
\description{This function calculates hierarchical allelic richness for populations structured in any hierarchy.}
\usage{HierAr(x, nreg, r, ncode)}
\arguments{
  \item{x}{The data file, genepop format.}
  \item{nreg}{Number of aggregate in your hierarchical structure.}
  \item{r}{The number of subaggregate in your hierarchical structure.}
  \item{ncode}{ Coding type of genetic data.}
}
\details{This function calculates the hierarchical allelic richness for metapopulations structured in any hierarchy. In this function, you can input your genepop file (x) and structure (nreg, r) of your data in the arguments, and then returns a matrix of your hierarchical Ar. This is very efficient than other functions.}
\value{
\item{Ar_pop}{Allelic richness for each population per locus}
\item{Ar_reg}{Allelic richness for each aggregate (eg. region)}
\item{Ar_overall}{Hierarchical allelic richness in different levels}
}
\references{
El Mousadik A. and Petit R.J. (1996) High level of genetic differentiation for allelic richness among populations of the argan tree argania spinosa skeels endemic to Morocco. Theoretical and Applied Genetics, 92:832-839

Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5(1), 184-186.
}
\author{qinxinghu@gmail.com}
\seealso{HierHe; HierFst; HierD}
\examples{
#we want to calculate the hierarchical allelic richness of
#16 populations structured in a 4 regions.
#Each region contains 7,4,2,3 populations.
download.file("https://sites.google.com/site/advancesinevolutionandgenetics
/Home/unit-3-1/Island.gen?attredirects=0&d=1",
"C:/Users/Public/Documents/Island.gen",cacheOK=FALSE)
HAr=HierAr("C:/Users/Public/Documents/Island.gen", nreg=4, r=c(7,4,2,3), ncode=3)
}
\keyword{Hierarchical genetic decomposition}
