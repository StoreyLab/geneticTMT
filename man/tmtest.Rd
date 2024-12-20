% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmtest.R
\name{tmtest}
\alias{tmtest}
\title{This function calculates the p-value for the TMT given trio genotypes and 
children's phenotype. The phenotype can be either quantitative or dichotomous.}
\usage{
tmtest(
  Y,
  trio_geno,
  is_phased = FALSE,
  test = c("TMT", "Permutation"),
  n_per = 100
)
}
\arguments{
\item{Y}{The n-length vector of children' quantitative trait}

\item{trio_geno}{A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
for is_phased. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
`TRUE` for is_phased. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
The first row is assumed maternal, the second row is assumed paternal, and
the third (or together with the fourth) row is assumed child's genotypes.}

\item{is_phased}{Label for unphased/phased genotype. If `FALSE` (default), the
child's genotype is unphased. Else, the child's genotype is phased.}

\item{test}{A character string of the test type.}

\item{n_per}{A numeric number for the permutation times if \code{test=}'Permutation'}
}
\value{
A numeric result for the p-value.
}
\description{
This function calculates the p-value for the TMT given trio genotypes and 
children's phenotype. The phenotype can be either quantitative or dichotomous.
}
\note{
Undesired dimensions of the genotype matrices and any test type other 
    than 'TMT' or 'Permutation' for \code{test} will result in errors.
}
\examples{
trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
y_quant <- c(11.8,1.7,3.6,12.2)
tmt_p <- tmtest(y_quant, trio_example, test='TMT')

tmt_p_perm <- tmtest(y_quant, trio_example, test='Permutation', n_per=1000)

}
