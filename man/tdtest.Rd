% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdtest.R
\name{tdtest}
\alias{tdtest}
\title{Transmission-Disequilibrium Test}
\usage{
tdtest(Y, trio_geno, is_phased = FALSE)
}
\arguments{
\item{Y}{The n-length vector of children' quantitative trait}

\item{trio_geno}{A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
for is_phased. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
`TRUE` for is_phased. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
The first row is assumed maternal, the second row is assumed paternal, and
the third (or together with the fourth) row is assumed child's genotype.}

\item{is_phased}{Label for unphased/phased genotype. If `FALSE` (default), the
child's genotype is unphased. Else, the child's genotype is phased.}
}
\value{
A list of p value and McNemar's test statistic.
}
\description{
This function calculates the \eqn{\chi^2} test statistic and p value for TDT.
}
\note{
Undesired dimensions of the genotype matrices will result in errors.
}
\examples{
trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE) 
y_binary <- c(0,1,1,1)
tdt_list <- tdtest(y_binary, trio_example)

}
