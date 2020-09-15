
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PSweight

This `PSweight` Package is to perform propensity score weighting
analysis for causal inference. Two main modules are included to assist
the design and analysis of observational studies. In the design module,
the `SumStat` function is used to generate distributional plots of the
estimated propensity scores and balance diagnostics after propensity
score weighting. The `summary` and `plot` functions are available to
tabulate and plot weighted balance statistics for visual comparisons. In
the analysis module, the `PSweight` function, the average potential
outcomes for each treatment group is estimated using weighting, and the
`summary` function generates point estimates, standard errors and
confidence intervals for the desired causal contrasts of interest. The
current version of `PSweight` package includes the following types of
weights: the overlap weights (ATO), the inverse probability of treatment
weights (ATE), the average treatment effect among the treated weights
(ATT), the matching weights (ATM) and the entropy weights (ATEN), and
allows for binary and multiple (categorical) treatments. In addition to
the simple weighting estimator, the package also implements the
augmented weighting estimator that combines weighting and outcome
regression. For binary outcomes, both the additive and ratio estimands
(causal relative risk and odds ratio) are considered, and variance is
estimated by either the sandwich method or nonparametric bootstrap. To
allow for additional flexibility in specifying the propensity score and
outcome models, the package can also work with user-supplied propensity
score estimates and outcome predictions through `ps.estimate` and
`out.estimate`, and provide a sandwich standard error that ignores the
variability in estimating these nuisances.

## Installation

You can install the released version of PSweight from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PSweight")
```

## Example

This is a basic example on design:

``` r
library(PSweight)
example("SumStat")
#> 
#> SumStt> data("psdata")
#> 
#> SumStt> # the propensity model
#> SumStt> ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
#> 
#> SumStt> # using SumStat to estimate propensity scores
#> SumStt> msstat <- SumStat(ps.formula, trtgrp="2", data=psdata,
#> SumStt+    weight=c("IPW","overlap","treated","entropy","matching"))
#> 
#> SumStt> summary(msstat)
#> unweighted result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.213  0.180 -0.074 0.270
#> cov2 -0.186  0.146 -0.080 0.332
#> cov3  0.009 -0.054 -0.019 0.062
#> cov4  0.150 -0.448  0.464 0.530
#> cov5  0.998  0.738  1.351 0.452
#> cov6  0.501  0.507  0.495 0.024
#> 
#> IPW result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.022  0.002  0.000 0.016
#> cov2 -0.023 -0.030 -0.027 0.007
#> cov3 -0.007 -0.030 -0.031 0.024
#> cov4  0.003  0.027  0.001 0.015
#> cov5  1.003  1.146  1.013 0.097
#> cov6  0.495  0.486  0.489 0.018
#> 
#> overlap result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.041 -0.043 -0.036 0.005
#> cov2 -0.042 -0.066 -0.055 0.025
#> cov3 -0.011 -0.013 -0.029 0.017
#> cov4  0.095  0.096  0.093 0.002
#> cov5  0.961  0.987  0.976 0.020
#> cov6  0.491  0.489  0.487 0.007
#> 
#> treated result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1  0.124  0.180  0.209 0.059
#> cov2  0.117  0.146  0.126 0.029
#> cov3  0.013 -0.054 -0.071 0.083
#> cov4 -0.435 -0.448 -0.454 0.011
#> cov5  0.747  0.738  0.764 0.026
#> cov6  0.498  0.507  0.486 0.043
#> 
#> entropy result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.039 -0.029 -0.027 0.009
#> cov2 -0.040 -0.056 -0.048 0.016
#> cov3 -0.009 -0.020 -0.028 0.018
#> cov4  0.060  0.070  0.058 0.007
#> cov5  0.986  1.072  1.001 0.061
#> cov6  0.493  0.487  0.488 0.012
#> 
#> matching result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.024 -0.043 -0.021 0.016
#> cov2 -0.027 -0.057 -0.045 0.031
#> cov3 -0.018 -0.008 -0.036 0.028
#> cov4  0.165  0.167  0.161 0.004
#> cov5  0.932  0.923  0.944 0.017
#> cov6  0.484  0.494  0.489 0.019
#> 
#> 
#> SumStt> # importing user-supplied propensity scores "e.h"
#> SumStt> fit <- nnet::multinom(formula=ps.formula, data=psdata, maxit=500, trace=FALSE)
#> 
#> SumStt> e.h <- fit$fitted.values
#> 
#> SumStt> varname <- c("cov1","cov2","cov3","cov4","cov5","cov6")
#> 
#> SumStt> msstat0 <- SumStat(zname="trt", xname=varname, data=psdata, ps.estimate=e.h,
#> SumStt+    trtgrp="2",  weight=c("IPW","overlap","treated","entropy","matching"))
#> 
#> SumStt> summary(msstat0)
#> unweighted result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.213  0.180 -0.074 0.270
#> cov2 -0.186  0.146 -0.080 0.332
#> cov3  0.009 -0.054 -0.019 0.062
#> cov4  0.150 -0.448  0.464 0.530
#> cov5  0.998  0.738  1.351 0.452
#> cov6  0.501  0.507  0.495 0.024
#> 
#> IPW result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.022  0.002  0.000 0.016
#> cov2 -0.023 -0.030 -0.027 0.007
#> cov3 -0.007 -0.030 -0.031 0.024
#> cov4  0.003  0.027  0.001 0.015
#> cov5  1.003  1.146  1.013 0.097
#> cov6  0.495  0.486  0.489 0.018
#> 
#> overlap result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.041 -0.043 -0.036 0.005
#> cov2 -0.042 -0.066 -0.055 0.025
#> cov3 -0.011 -0.013 -0.029 0.017
#> cov4  0.095  0.096  0.093 0.002
#> cov5  0.961  0.987  0.976 0.020
#> cov6  0.491  0.489  0.487 0.007
#> 
#> treated result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1  0.124  0.180  0.209 0.059
#> cov2  0.117  0.146  0.126 0.029
#> cov3  0.013 -0.054 -0.071 0.083
#> cov4 -0.435 -0.448 -0.454 0.011
#> cov5  0.747  0.738  0.764 0.026
#> cov6  0.498  0.507  0.486 0.043
#> 
#> entropy result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.039 -0.029 -0.027 0.009
#> cov2 -0.040 -0.056 -0.048 0.016
#> cov3 -0.009 -0.020 -0.028 0.018
#> cov4  0.060  0.070  0.058 0.007
#> cov5  0.986  1.072  1.001 0.061
#> cov6  0.493  0.487  0.488 0.012
#> 
#> matching result
#>      Mean 1 Mean 2 Mean 3   SMD
#> cov1 -0.024 -0.043 -0.021 0.016
#> cov2 -0.027 -0.057 -0.045 0.031
#> cov3 -0.018 -0.008 -0.036 0.028
#> cov4  0.165  0.167  0.161 0.004
#> cov5  0.932  0.923  0.944 0.017
#> cov6  0.484  0.494  0.489 0.019
```

This is a basic example on analysis:

``` r
example("PSweight")
#> 
#> PSwght> data("psdata")
#> 
#> PSwght> # the propensity and outcome models
#> PSwght> ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
#> 
#> PSwght> out.formula<-Y~cov1+cov2+cov3+cov4+cov5+cov6
#> 
#> PSwght> # without augmentation
#> PSwght> ato1<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,weight = 'overlap')
#> 
#> PSwght> summary(ato1)
#> Original group value:  1, 2, 3 
#> 
#> contrast 
#>             1  2 3
#> Contrast 1 -1  1 0
#> Contrast 2 -1  0 1
#> Contrast 3  0 -1 1
#> 
#>             Estimate Std.Error   Lower.CL   Upper.CL
#> Contrast 1 -1.241607 0.1673446 -1.5695962 -0.9136176
#> Contrast 2  1.124818 0.1709932  0.7896776  1.4599588
#> Contrast 3  2.366425 0.2585402  1.8596956  2.8731546
#> 
#> PSwght> # augmented weighting estimator
#> PSwght> ato2<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,
#> PSwght+                augmentation = TRUE,out.formula = out.formula,family = 'gaussian',weight = 'overlap')
#> 
#> PSwght> summary(ato2)
#> Original group value:  1, 2, 3 
#> 
#> contrast 
#>             1  2 3
#> Contrast 1 -1  1 0
#> Contrast 2 -1  0 1
#> Contrast 3  0 -1 1
#> 
#>             Estimate Std.Error   Lower.CL   Upper.CL
#> Contrast 1 -1.238193 0.1217147 -1.4767496 -0.9996368
#> Contrast 2  1.165513 0.1348439  0.9012242  1.4298027
#> Contrast 3  2.403707 0.1829589  2.0451138  2.7622995
```
