# geneticTMT

![R-CMD-check.yaml](https://github.com/StoreyLab/geneticTMT/actions/workflows/R-CMD-check.yaml/badge.svg)

The `geneticTMT` R package implements the Transmission Mean Test (TMT) and the Transmission Disequilibrium Test (TDT) to infer causal genotype-phenotype relationships for population-sampled nuclear families.

## Installation

```R
install.packages("devtools")
devtools::install_github("StoreyLab/geneticTMT")
```

## Examples

### Input trio data

`geneticTMT` handles trio genotypes and phenotypes such as the example below:

```R
library(geneticTMT)
# load example trio genotypes and phenotypes
data(triolist)
# trio genotypes can be either phased or unphased
trio_unphased <- matrix(c(triolist$geno_m[1,],
                          triolist$geno_p[1,],
                          triolist$geno_c[1,]),
                        nrow=3, byrow=TRUE)
# phased data include parental transmitted alleles                      
trio_phased <- matrix(c(triolist$geno_m[1,],
                        triolist$geno_p[1,],
                        triolist$allele_m[1,],
                        triolist$allele_p[1,]),
                      nrow=4, byrow=TRUE)

```

### TMT

Calculate the TMT statistic as an unbiased estimate of the average causal effect of the target locus on the child's trait.

```R
dtmt <- calc_dtmt(triolist$y, trio_unphased)
```

Calculate p-value via the main test function

```R
ptmt <- tmtest(triolist$y, trio_unphased, test='TMT')
```


### TDT

TDT requires dichotomous phenotypes.

```R
# convert quantitative to dichotomous
y_binary <- triolist$y
y_binary <- ifelse(y_binary>median(y_binary), 1, 0)
```

Calculate the TDT statistic as an unbiased estimate of the average causal effect of the target locus on the child's trait.

```R
dtdt <- calc_dtdt(y_binary, trio_unphased)
```

Calculate p-value via the main test function

```R
ptdt <- tdtest(y_binary, trio_unphased)$p
```


## Citations
