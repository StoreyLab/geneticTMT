# The Generalized Transmission-Mean Test Statistic and Sampling Variance Estimate
#'
#' This function computes the unbiased gTMT estimand given nuclear family
#' genotypes and offspring phenotype. The phenotype can be either quantitative
#' or dichotomous.
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
#' @param estimate_variance Label for estimate sampling variance of the gTMT
#'     statistic or not. If `FALSE` (default), only calculate the gTMT statistic
#'     child's genotype is unphased. Else, the child's genotype is phased.
#'
#' @return A numeric value of the gTMT statistics if \code{estimate_variance} is
#'     `TRUE`, or a two-length vector containing the numeric value of the gTMT
#'     statistic and corresponding sampling variance estimate if
#'     \code{estimate_variance} is FALSE.
#'
#' @note Undesired dimensions of the genotype matrices will result in errors.
#'
#' @examples
#' y_example <- c(9.3,0.9,20.1,7.8,17.5,-0.7,2.4,23.6,-1.1,6.0,0.2,7.7,8.1) # offspring phenotype
#' g_example <- c(1,0,2,1,2,0,0,2,0,1,0,1,1) # offspring genotype
#' zm_example <- c(1,2,1,0,1,1,0,1) # maternal genotype
#' zp_example <- c(1,1,1,0,1,0,1,2) # paternal genotype
#' fam_id_o <- c(1,1,1,2,2,3,4,5,5,6,7,7,8) # offspring family id
#' fam_id_m <- c(1,2,3,4,5,6,7,8) # maternal family id
#' fam_id_p <- c(1,2,3,4,5,6,7,8) # paternal family id
#' dgtmt <- calc_dgtmt(y_example,
#'                     g_example, zm_example, zp_example,
#'                     fam_id_o, fam_id_m, fam_id_p)
#' d_var <- calc_dgtmt(y_example,
#'                     g_example, zm_example, zp_example,
#'                     fam_id_o, fam_id_m, fam_id_p,
#'                     estimate_variance=TRUE)
#'
#' @export
calc_dgtmt <- function(Y, G, Zm, Zp, fam_id_o, fam_id_m, fam_id_p,
                       estimate_variance=FALSE){
  # construct assignment vectors
  fam_geno <- matrix(c(Zm[match(fam_id_o, fam_id_m)], Zp[match(fam_id_o, fam_id_p)],G),
                     nrow=3, byrow=TRUE)
  list_assign <- assignment_index(fam_geno)
  W0 <- list_assign$W0
  W1 <- list_assign$W1
  # estimate mu to center trait value
  mu_hat <- calc_mu_hat(Y, W0, W1)
  # sizes of treatment and control
  N <- sum(W1) + sum(W0)
  # calculate the tmt estimand
  dgtmt <- ifelse(N>0, 2/N*(sum((Y-mu_hat)*W1) - sum((Y-mu_hat)*W0)), NaN)
  if(estimate_variance){
    T0 <- (W0==1) & (W1==0)
    T1 <- (W0==0) & (W1==1)
    T00 <- (W0==2) & (W1==0)
    T11 <- (W0==0) & (W1==2)
    N0 <- sum(T0)
    N1 <- sum(T1)
    N00 <- sum(T00)
    N11 <- sum(T11)
    y <- as.numeric(Y)
    mu_hat_0 <- ifelse(N0>0, sum(y[T0])/N0, NaN)
    mu_hat_1 <- ifelse(N1>0, sum(y[T1])/N1, NaN)
    mu_hat_00 <- ifelse(N00>0, 2*sum(y[T00])/N00, NaN)
    mu_hat_11 <- ifelse(N11>0, 2*sum(y[T11])/N11, NaN)
    sigma2_hat_0 <- ifelse(N0>1, sum((y[T0]-mu_hat_0)**2)/(N0-1), 0)
    sigma2_hat_1 <- ifelse(N1>1, sum((y[T1]-mu_hat_1)**2)/(N1-1), 0)
    sigma2_hat_00 <- ifelse(N00>1, sum((2*y[T00]-mu_hat_00)**2) / (N00-1), 0)
    sigma2_hat_11 <- ifelse(N11>1, sum((2*y[T11]-mu_hat_11)**2) / (N11-1), 0)
    if((N0<1) & (N1<1) & (N00<1) & (N11<1)){
      varhat <- NaN
      warning('no heterozygous parental genotype')
    } else {
      varhat <- 4/(N**2)*(N0*sigma2_hat_0 + N1*sigma2_hat_1 + N00*sigma2_hat_00 + N11*sigma2_hat_11)
    }
    stat_out <- c(dgtmt, varhat)
    names(stat_out) <- c('dgtmt','varhat')
    return(stat_out)
  } else {
    return(dgtmt)
  }
}
