# geneticTMT

![R-CMD-check.yaml](https://github.com/StoreyLab/geneticTMT/actions/workflows/R-CMD-check.yaml/badge.svg)

The `geneticTMT` R package implements the Transmission Mean Test (TMT) for population-sampled parent–child trios and the generalized Transmission Mean Test (gTMT) for population-sampled nuclear families to identify causal genotype–phenotype relationships.

## Installation

```R
install.packages("pak")
pak::pak("StoreyLab/geneticTMT")
```

## Examples

### Input trio data

`geneticTMT` handles trio genotypes and phenotypes such as the example below:

```R
library(geneticTMT)
# load example trio genotypes and phenotypes
data(triolist)
# load example nuclear family genotypes and offspring phenotypes
data(familydata)
```

### TMT

Prepare trio data at a target locus.

```R
# trio genotypes can be either phased or unphased
trio_unphased <- matrix(c(triolist$geno_m[1,],
                          triolist$geno_p[1,],
                          triolist$geno_c[1,]),
                        nrow=3, byrow=TRUE)
```

Calculate the TMT statistic as an unbiased estimate of the average causal effect of the target locus on the child's trait.

```R
dtmt <- calc_dtmt(triolist$y, trio_unphased)
```

Calculate p-value via the main test function

```R
ptmt <- tmtest(triolist$y, trio_unphased, test='TMT')
```

### gTMT

Calculate the gTMT statistic as an unbiased estimate of the average causal effect of the target locus on offspring phenotype.

```R
# calculate the gTMT statistic across 20 example loci
d_per_locus <- function(i){
  d_gtmt <- calc_dgtmt(Y=familydata$pheno_offspring, G=familydata$geno_offspring[i,], 
                       Zm=familydata$geno_maternal[i,], Zp=familydata$geno_paternal[i,], 
                       fam_id_o=familydata$fam_id_offspring, fam_id_m=familydata$fam_id_maternal, 
                       fam_id_p=familydata$fam_id_paternal, estimate_variance=FALSE)
  return(d_gtmt)
}
gtmt_stat <- unlist(lapply((1:20), d_per_locus))
```

Calculate the p-values.

```R
# test all 20 example loci
p_per_locus <- function(i){
  p_gtmt <- gtmtest(Y=familydata$pheno_offspring, G=familydata$geno_offspring[i,], 
                    Zm=familydata$geno_maternal[i,], Zp=familydata$geno_paternal[i,], 
                    fam_id_o=familydata$fam_id_offspring, fam_id_m=familydata$fam_id_maternal, 
                    fam_id_p=familydata$fam_id_paternal, test='gTMT')
  return(p_gtmt)
}
pvalue_gtmt <- unlist(lapply(1:20, p_per_locus))
```


## Citations

Yushi Tang, Irineo Cabreros, John D Storey.  2026.  "Identifying causal genotype–phenotype relationships for population-sampled parent–child trios." Genetic Epidemiology 50(1): e70027. [doi:10.1002/gepi.70027](https://doi.org/10.1002/gepi.70027). bioRxiv [doi:10.1101/2024.12.10.627752](https://doi.org/10.1101/2024.12.10.627752) 2024-12-11.

Yushi Tang, John D Storey.  2025.  "A generalized test of genotype–phenotype causality in population-sampled nuclear families." bioRxiv [doi:10.64898/2025.12.29.696865](https://doi.org/10.64898/2025.12.29.696865) 2025-12-29.
