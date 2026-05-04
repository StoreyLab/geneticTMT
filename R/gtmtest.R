# The Generalized Transmission-Mean Test (gTMT)
#'
#' This function calculates the p-value for the gTMT given nuclear family
#' genotypes and offspring phenotype, allowing multiple offspring per family.
#' The phenotype can be either quantitative or dichotomous.
#'
#' @param Y The J-length vector of offspring phenotypes
#' @param G A K-length vector of offspring genotype at a genetic loci.
#'     Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#' @param Zm A J-length vector of maternal genotype at a genetic loci.
#'     Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#' @param Zp A J-length vector of paternal genotype at a genetic loci.
#'     Must contain only values \eqn{(0, 1, 2)}{0, 1, 2}.
#' @param fam_id_o A K-length vector of family id for offspring.
#' @param fam_id_m A J-length vector of maternal family id.
#' @param fam_id_p A J-length vector of paternal family id.
#' @param test A character string of the test type.
#' @param n_per A numeric number for the permutation times if \code{test=}'Permutation'
#'
#' @importFrom stats pnorm
#'
#' @return A numeric result for the p-value.
#'
#' @note Undesired dimensions of the genotype matrices and any test type other
#'     than 'gTMT' or 'Permutation' for \code{test} will result in errors.
#'
#' @examples
#' y_example <- c(9.3,0.9,20.1,7.8,17.5,-0.7,2.4,23.6,-1.1,6.0,0.2,7.7,8.1) # offspring phenotype
#' g_example <- c(1,0,2,1,2,0,0,2,0,1,0,1,1) # offspring genotype
#' zm_example <- c(1,2,1,0,1,1,0,1) # maternal genotype
#' zp_example <- c(1,1,1,0,1,0,1,2) # paternal genotype
#' fam_id_o <- c(1,1,1,2,2,3,4,5,5,6,7,7,8) # offspring family id
#' fam_id_m <- c(1,2,3,4,5,6,7,8) # maternal family id
#' fam_id_p <- c(1,2,3,4,5,6,7,8) # paternal family id
#'
#' gtmt_p <- gtmtest(y_example,
#'                   g_example, zm_example, zp_example,
#'                   fam_id_o, fam_id_m, fam_id_p, test='gTMT')
#'
#' gtmt_p_perm <- gtmtest(y_example,
#'                        g_example, zm_example, zp_example,
#'                        fam_id_o, fam_id_m, fam_id_p,
#'                        test='Permutation', n_per=1000)
#'
#'
#' @export
gtmtest <- function(Y, G, Zm, Zp, fam_id_o, fam_id_m, fam_id_p,
                    test=c('gTMT','Permutation'), n_per=100){
  test <- match.arg(test)
  if(test == 'gTMT'){
    d_var <- calc_dgtmt(Y, G, Zm, Zp, fam_id_o, fam_id_m, fam_id_p, estimate_variance=TRUE)
    dgtmt <- as.numeric(d_var['dgtmt'])
    varhat <- as.numeric(d_var['varhat'])
    pvalue <- ifelse(is.na(varhat), NaN, 2 * stats::pnorm(-abs(dgtmt/varhat)))
  } else if (test == 'Permutation'){
    dgtmt <- calc_dgtmt(Y, G, Zm, Zp, fam_id_o, fam_id_m, fam_id_p, estimate_variance=FALSE)
    if( is.na(dgtmt) ){
      pvalue <- NaN
    } else {
      perm_loop <- function(i){
        set.seed(i)
        dgtmt_per <- calc_dgtmt(sample(Y, size=length(Y), replace=FALSE),
                                G, Zm, Zp, fam_id_o, fam_id_m, fam_id_p)
        return( abs(dgtmt_per) >= abs(dgtmt) )
      }
      p_perm <- lapply(c(1:n_per), FUN=perm_loop)
      p_perm <- unlist(p_perm)
      pvalue <- sum(p_perm) / length(p_perm)
    }
  }
  return(pvalue)
}
