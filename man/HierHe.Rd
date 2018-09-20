\name{HierHe}
\alias{HierHe}
\title{Function to calculate hierarchical heterozygosity.}
\description{This function calculates and decomposes hierarchical heterozygosity into different levels. We can obtain allelic diversity (He) in any hierarchy.}
\usage{HierHe(x, nreg, r, ncode)}
\arguments{
  \item{x}{Genetic data, the format should be genepop file with proper extention.}
  \item{nreg}{Number of aggregate in the hierarchy.}
  \item{r}{The number of subaggregate in the hierarchy.}
  \item{ncode}{ Coding type of genetic data, whether 3, or 2.}
}
\details{In line with other genetic metrics, this function obatains hierarchical expected heterozygosity (He) for metapopulations. The result gives us hierarchical allelic diversity at different levels as you set in the input argument.}
\value{
  \item{HierHe_perloc}{List of hierarchical allelic diversity (heterozygosity) in different levels per locus.}
  \item{HieHr_overall}{Overall (mean) hierarchical allelic diversity (heterozygosity) in different levels.}
}
\references{
Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5(1), 184-186.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), 1403-1405.

Nei, M. (1977). F-statistics and analysis of gene diversity in subdivided populations. Annals of human genetics, 41(2), 225-233.
}
\author{qinxinghu@gmail.com}
\note{}
\seealso{HierAr,HierdetaD,HierFst}
\examples{
HeirHe=HierdetaD("Island.gen",nreg=4,r=c(7,4,2,3),ncode=3)}
\keyword{ "HierHe" }

