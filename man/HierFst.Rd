\name{HierFst}
\alias{HierFst}
\title{Function to calculate hierarchical genetic differentiation (Fst) (Weir B.S., 1996; Yang R.C.,1998).}
\description{This function calculates hierarchical genetic differentiation (Fst) based on Weir B.S., 1996 and Yang R.C., 1998). The hierarchy that set in the input arguments (nreg, r) will automatically generate your hierarchical structure. Based on the parameters you set, the result returns to a F-statistics matrix displaying one level (rows) of Fst within a higher level hierarchy (columns).}
\usage{HierFst(x, nreg, r, ncode)}
\arguments{
  \item{x}{The name of the data file, or the string giving the path to the file. The file should be a genepop object, with with the appropriate extension.}
  \item{nreg}{The number of aggregate in the hierarchy.}
  \item{r}{A vector specifying number of subaggregate (eg.subpopulations) within each aggregate.}
  \item{ncode}{An integer indicating the number of characters used to code an allele.}
}
\details{This function was modified and improved from the varcomp.glob function from the hierfstat package (Jerome Goudet & Thibaut Jombart 2015). It was optimized to easily set the hierarchical structure while varcom.glob requires the structure data should be already correctly formated in accordance with hierfstat data. In this function, you can input your genetic file (x) and set the structure (nreg, r) of your data in the arguments, and then it returns a matrix of hierarchical Fst.}
\value{The result returns a matrix of hierarchical F-statistics type-coefficients at different levels.}
\references{
Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

Yang, R.C. (1998). Estimating hierarchical F-statistics. Evolution 52 (4):950-956

Goudet J. (2005). Hierfstat, a package for R to compute and test variance components and Fs-tatistics.
Molecular Ecology Notes. 5:184-186}
\author{ qinxinghu@gmail.com}
\note{Note that you should arrange your population in the order corresponding to your structure set in nreg.}
\seealso{HierAr}
\examples{#
infile <- system.file('extdata', 'Island.gen',package='HierDpart')
locinfile <- 'Island.gen' ## file in user's directory not in R's extdata directory
check <- file.copy(infile,locinfile,overwrite=TRUE)
 hFst=HierFst(infile, nreg=4, r=c(7,4,2,3), ncode=3)
 print(hFst)}
\keyword{Hierarchical genetic decomposition}
