#' Example Trio Data
#' 
#'
#' @format A list of matrices:
#' \describe{
#'  \item{geno_m}{A matrix of the maternal genotypes.}
#'  \item{geno_p}{A matrix of the paternal genotypes.}
#'  \item{geno_c}{A matrix of the child's genotypes}
#'  \item{allele_m}{A matrix of the maternal transmitted alleles.}
#'  \item{allele_p}{A matrix of the paternal transmitted alleles.}
#'  \item{y}{A single row matrix of the child's phenotypes.}
#' }
#' @source Simulated data
"triolist"


#' Example Nuclear Family Data.
#' 
#'
#' @format A list of matrices:
#' \describe{
#'   \item{pheno_offspring}{Offspring Phenotype}
#'   \item{geno_offspring}{Offspring Genotypes}
#'   \item{geno_maternal}{Maternal Genotypes}
#'   \item{geno_paternal}{Paternal Genotypes}
#'   \item{fam_id_maternal}{Maternal Family ID}
#'   \item{fam_id_paternal}{Paternal Family ID}
#'   \item{fam_id_offspring}{Offspring Family ID}
#' }
#'
#' @source Simulated data
"familydata"

NULL
