
# HierDpart

HierDpart is made for calculating hierarchical diversity and differentiation based on a unifying framework (Oscar, et al, 2018), whcih can investigate biological variation across metrixs and scales, from genes to ecosystem.

HierDpart is free R package: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 HierDpart works for calculating and decomposting real diversity (effective numbers qD, q=1,2,3 etc) and new differentiation (detaD) for any number of hierarchical levels based on a recently hieacrchical framework( Oscar,E et al, 2018).
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License along with HierDpart.  If not, see <http://www.gnu.org/licenses/>.

To install the development version of HierDpart 

You will need the package devtools to be able to install the devel version of HierDpart. 
>>>>>>> 35fe7ed9faf7bb7d8e3ed754a7f95ac57d5e787d
Install devtools:

install.packages("devtools")

To install HierDpart:

library(devtools)

install_github("xinghuq/HierDpart")

library("HierDpart")


# References
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.

Chao, A. and Chiu, C.-H. (2017) iDIP (Information-based Diversity Partitioning) Online: Software for partitioning Shannon diversity and phylogenetic diversity under multi-level hierarchical structures. Program and User's Guide published at http://chao.stat.nthu.edu.tw/wordpress/software_download/.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), 1403-1405.

Keenan, K., McGinnity, P., Cross, T. F., Crozier, W. W., & ProdÃ¶hl, P. A. (2013). diveRsity: an R package for the estimation and exploration of population genetics parameters and their associated errors. Methods in Ecology and Evolution, 4(8), 782-788.

Marcon, E., & HÃ©rault, B. (2015). entropart: An R package to measure and partition diversity. Journal of Statistical Software, 67(8).

Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical Fâ€statistics. Molecular Ecology Notes, 5(1), 184-186.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), 1403-1405.

