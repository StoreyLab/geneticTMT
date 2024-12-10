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

# convert quantitative to dichotomous
y_binary <- y_quantitative <- triolist$y
y_binary <- ifelse(y_binary>median(y_binary), 1, 0)

test_that("expect precomputed values: dtdt", {
  # expect error if genotypes matrices dimensions do not match the phase/unphase status
  expect_error( tdtest(y_binary, trio_phased, is_phased=FALSE) )
  expect_error( tdtest(y_binary, trio_unphased, is_phased=TRUE) )
  # expect error if phenotypes not binary
  expect_error( tdtest(y_quantitative, trio_unphased) )
  # check tdt statistics
  expect_equal( calc_dtdt(y_binary, trio_unphased),  dtdt)  
  # check p values
  expect_equal( tdtest(y_binary, trio_unphased)$p, ptdt)
  # unphased and phased genotypes should yield the same result
  expect_equal( calc_dtdt(y_binary, trio_unphased, is_phased=FALSE),  
                calc_dtdt(y_binary, trio_phased, is_phased=TRUE))
  expect_equal( tdtest(y_binary, trio_unphased, is_phased=FALSE), 
                tdtest(y_binary, trio_phased, is_phased=TRUE))
})