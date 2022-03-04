
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saeHB.hnb

<!-- badges: start -->
<!-- badges: end -->

This package provides a function and datasets for area level of Small
Area Estimation under Hurdle Negative Binomial Model using Hierarchical
Bayesian (HB) Method.

## Author

Raka Ikmana, Azka Ubaidillah

## Maintaner

Raka Ikmana <221810548@stis.ac.id>

## Function

-   `HurdleNB()` This function gives small area estimator under Hurdle
    Negative Binomial Model and is implemented to variable of
    interest (y) that assumed to be a HNB Distribution. The value of
    variable of interest must be a non-negative data count.

## Installation

You can install the development version of `saeHB.hnb` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rakaikmana/saeHB.hnb")
```

## Example

This is a basic example of using `HurdlenNB()` function to make an
estimate based on synthetic data in this package

``` r
library(saeHB.hnb)
## For data without any non-sampled area
data(dataHNB)       # Load dataset

## For data with non-sampled area use dataHNBNs
## Fitting model
result <- HurdleNB(y ~ x1 + x2, data=dataHNB)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 50
#>    Unobserved stochastic nodes: 159
#>    Total graph size: 2128
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 50
#>    Unobserved stochastic nodes: 159
#>    Total graph size: 2128
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 50
#>    Unobserved stochastic nodes: 159
#>    Total graph size: 2128
#> 
#> Initializing model
```

Small Area mean Estimates

``` r
result$Est
```

Estimated model coefficient

``` r
result$coefficient
```

Estimated random effect variances

``` r
result$refVar
#>             [,1]
#> a.var.u 0.703366
#> a.var.v 2.417690
```

## References

-   Andika, A., Abdullah, S., & Nurrohmah, S. (2019). “Hurdle Negative
    Binomial Regression Model”. Proceeding of ICSA 2019, p: 57-68.
    <doi:10.29244/icsa.2019.pp57-68>.
-   Hilbe, J. M. (2011). Negative Binomial Regression 2nd Edition. New
    York : Cambridge University Press. <doi:10.1017/CBO9780511973420>.
-   Ntzoufras, I. (2009). Bayesian Modelling Using WinBUGS. New Jersey :
    John Wiley & Sons, Inc.
-   Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
    Jersey : John Wiley & Sons, Inc. <doi:10.1002/9781118735855>.
