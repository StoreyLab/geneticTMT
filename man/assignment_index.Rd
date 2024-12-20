% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assignment_index.R
\name{assignment_index}
\alias{assignment_index}
\title{Assignment Index}
\usage{
assignment_index(trio_geno, is_phased = FALSE)
}
\arguments{
\item{trio_geno}{A \eqn{3 \times n}{3 x n} matrix of trio genotype if `FALSE`
for \code{is_phased}. Or a \eqn{4 \times n}{4 x n} matrix of trio genotype if
`TRUE` for \code{is_phased}. Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
The first row is assumed maternal, the second row is assumed paternal, and
the third (or together with the fourth) row is assumed child's genotypes.}

\item{is_phased}{Label for unphased/phased genotype. If `FALSE` (default), the
child's genotype is unphased. Else, the child's genotype is phased.}
}
\value{
A named list of two vectors containing indicators for assigning to 
    treatment group 0 or treatment group 1.
}
\description{
This function constructs assignment indices based on trio genotypes.
}
\note{
Undesired dimensions of the genotype matrices will result in errors.
}
\examples{
trio_example <- matrix(c(1,1,2,0, 1,2,1,0, 1,1,1,0), nrow=3, byrow=TRUE)
list_assign <- assignment_index(trio_example, is_phased=FALSE)

}
