
# LNPar

The goal of LNPar is to estimate and test for a Lognormal-Pareto Mixture

## Installation

You can install the development version of LNPar from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("marco-bee/lnpar")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(LNPar)
ySim <- rLnormParMix(100,.5,0,1,4,1.5)
mixFit <- LPfit(ySim,90,0)
```
