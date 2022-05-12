
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

## Link to Paper

You can access the [paper](https://arxiv.org/pdf/2010.08893.pdf) with
examples on real-world dateset.

## Installation

You can install the released version of PSweight from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PSweight")
```

## Update

The version 1.1.7 includes module for cluster design. Please check out the help page for PSweight_cl and SumStat_cl.


## Downloads

![Downloads
Status](https://cranlogs.r-pkg.org/badges/grand-total/PSweight)

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
#> SumStt> #summary(msstat)
#> SumStt> 
#> SumStt> # importing user-supplied propensity scores "e.h"
#> SumStt> # fit <- nnet::multinom(formula=ps.formula, data=psdata, maxit=500, trace=FALSE)
#> SumStt> # e.h <- fit$fitted.values
#> SumStt> # varname <- c("cov1","cov2","cov3","cov4","cov5","cov6")
#> SumStt> # msstat0 <- SumStat(zname="trt", xname=varname, data=psdata, ps.estimate=e.h,
#> SumStt> #  trtgrp="2",  weight=c("IPW","overlap","treated","entropy","matching"))
#> SumStt> # summary(msstat0)
#> SumStt> 
#> SumStt> 
#> SumStt> 
#> SumStt>
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
#> 
#> Closed-form inference: 
#> 
#> Original group value:  1, 2, 3 
#> 
#> Contrast: 
#>             1  2 3
#> Contrast 1 -1  1 0
#> Contrast 2 -1  0 1
#> Contrast 3  0 -1 1
#> 
#>            Estimate Std.Error      lwr      upr  Pr(>|z|)    
#> Contrast 1 -1.24161   0.16734 -1.56960 -0.91362 1.177e-13 ***
#> Contrast 2  1.12482   0.17099  0.78968  1.45996 4.764e-11 ***
#> Contrast 3  2.36643   0.25854  1.85970  2.87315 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> PSwght> # augmented weighting estimator, takes longer time to calculate sandwich variance
#> PSwght> # ato2<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,
#> PSwght> #              augmentation = TRUE,out.formula = out.formula,family = 'gaussian',weight = 'overlap')
#> PSwght> # summary(ato2)
#> PSwght> 
#> PSwght> 
#> PSwght> 
#> PSwght>
```
