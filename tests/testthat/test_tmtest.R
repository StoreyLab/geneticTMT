# load example data
load('./triolist.RData')
load('./pre.RData')

trio_phased <- matrix(c(triolist$geno_m[1,], 
                        triolist$geno_p[1,], 
                        triolist$allele_m[1,], 
                        triolist$allele_p[1,]),
                      nrow=4, byrow=TRUE)

trio_unphased <- matrix(c(triolist$geno_m[1,], 
                          triolist$geno_p[1,], 
                          triolist$geno_c[1,]),
                        nrow=3, byrow=TRUE)

list_assign <- assignment_index(trio_unphased)

test_that("expect precomputed values: dtmt", {
  # expect error if genotypes matrices dimensions do not match the phase/unphase status
  expect_error( tmtest(triolist$y, trio_phased, is_phased=FALSE, test='TMT') )
  expect_error( tmtest(triolist$y, trio_unphased, is_phased=TRUE, test='TMT') )
  # expect error if test type out of categories
  expect_error( tmtest(triolist$y, trio_unphased, test='TDT') )
  # check population mean estimate
  expect_equal( calc_mu_hat(triolist$y, list_assign$W0, list_assign$W1), mu_hat)
  # check tmt statistics
  expect_equal( calc_dtmt(triolist$y, trio_unphased),  dtmt)
  expect_equal( calc_dtmt_nc(triolist$y, trio_unphased),  dtmtnc)
  expect_equal( calc_dtmt_pc(triolist$y, trio_unphased, mu=126),  dtmtpc)
  # check sampling variance estimate
  expect_equal( calc_varhat_dtmt(triolist$y, trio_unphased), varhat_dtmt )
  # check p values
  expect_equal( tmtest(triolist$y, trio_unphased, test='TMT'), ptmt)
  expect_equal( tmtest(triolist$y, trio_unphased, test='TMT'), ptmt_perm)
  # unphased and phased genotypes should yield the same result
  expect_equal( assignment_index(trio_unphased, is_phased=FALSE)$W0,
                assignment_index(trio_phased, is_phased=TRUE)$W0)
  expect_equal( assignment_index(trio_unphased, is_phased=FALSE)$W1,
                assignment_index(trio_phased, is_phased=TRUE)$W1)
  expect_equal( calc_dtmt(triolist$y, trio_unphased, is_phased=FALSE),  
                calc_dtmt(triolist$y, trio_phased, is_phased=TRUE))
  expect_equal( calc_varhat_dtmt(triolist$y, trio_unphased, is_phased=FALSE),  
                calc_varhat_dtmt(triolist$y, trio_phased, is_phased=TRUE))
  expect_equal( tmtest(triolist$y, trio_unphased, is_phased=FALSE, test='TMT'), 
                tmtest(triolist$y, trio_phased, is_phased=TRUE, test='TMT'))
})