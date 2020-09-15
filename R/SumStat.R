#' Calculate summary statistics for propensity score weighting
#'
#' \code{SumStat} is used to generate distributional plots of the estimated propensity scores and balance
#' diagnostics after propensity score weighting.
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the propensity score model to be fitted. Additional details of model specification are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param ps.estimate an optional matrix or data frame containing estimated (generalized) propensity scores for each observation. Typically, this is an N by J matrix, where N is the number of observations and J is the total number of treatment levels. Preferably, the column names of this matrix should match the names of treatment level, if column names are missing or there is a mismatch, the column names would be assigned according to the alphabatic order of treatment levels. A vector of propensity score estimates is also allowed in \code{ps.estimate}, in which case a binary treatment is implied and the input is regarded as the propensity to receive the last category of treatment by alphabatic order, unless otherwise stated by \code{trtgrp}.
#' @param trtgrp an optional character defining the "treated" population for estimating the average treatment effect among the treated (ATT). Only necessary if \code{weight = "treated"}. This option can also be used to specify the treatment (in a two-treatment setting) when a vector argument is supplied for \code{ps.estimate}. Default value is the last group in the alphebatic order.
#' @param Z an optional vector specifying the values of treatment, only necessary when the covariate matrix \code{covM} is provided instead of \code{data}.
#' @param covM an optional covariate matrix or data frame including covariates, their interactions and higher-order terms. When the covariate matrix \code{covM} is provided, the balance statistics are generated according to each column of this matrix.
#' @param zname an optional character specifying the name of the treatment variable in \code{data}.
#' @param xname an optional vector of characters including the names of covariates in \code{data}.
#' @param data an optional data frame containing the variables in the propensity score model. If not found in data, the variables are taken from \code{environment(formula)}.
#' @param weight a character or vector of characters including the types of weights to be used. \code{"IPW"} specifies the inverse probability weights for estimating the average treatment effect among the combined population (ATE). \code{"treated"} specifies the weights for estimating the average treatment effect among the treated (ATT). \code{"overlap"} specifies the (generalized) overlap weights for estimating the average treatment effect among the overlap population (ATO), or population at clinical equipoise. \code{"matching"} specifies the matching weights for estimating the average treatment effect among the matched population (ATM). \code{"entropy"} specifies the entropy weights for the average treatment effect of entropy weighted population (ATEN). Default is \code{"overlap"}.
#' @param delta trimming threshold for estimated (generalized) propensity scores. Should be no larger than 1 / number of treatment groups. Default is 0, corresponding to no trimming.
#'
#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable (identical to the variable name used to specify \code{zname}) and \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment}. \code{ps.formula} specifies logistic or multinomial logistic
#' models for estimating the propensity scores, when \code{ps.estimate} is \code{NULL}.
#'
#' When comparing two treatments, \code{ps.estimate} can either be a vector or a two-column matrix of estimated
#' propensity scores. If a vector is supplied, it is assumed to be the propensity scores to receive the treatment, and
#' the treatment group corresponds to the last group in the alphebatic order, unless otherwise specified by \code{trtgrp}.
#' When comparing multiple (J>=3) treatments, \code{ps.estimate} needs to be specified as an N by J matrix,
#' where N indicates the number of observations, and J indicates the total number of treatments.
#' This matrix specifies the estimated generalized propensity scores to receive each of the J treatments.
#' In general, \code{ps.estimate} should have column names that indicate the level of the treatment variable,
#' which should match the levels given in \code{Z}.
#' If column names are empty or there is a mismatch, the column names will be created following
#' the alphebatic order of treatmentlevels. The rightmost coulmn of \code{ps.estimate} is then assumed
#' to be the treatment group when estimating ATT (\code{"treated"}). \code{trtgrp} can also be used to specify the treatment
#' group for estimating ATT.
#'
#' To generate balance statistics, one can directly specify \code{Z} and \code{covM} to indicate the treatment levels and
#' covariate matrix. Alternatively, one can supply \code{data}, \code{zname}, and \code{xname} to indicate the
#' same information. When both are specified, the function will prioritize inputs from \code{Z} and \code{covM}.
#' When \code{ps.estimate} is not \code{NULL}, argument \code{zname}.
#'
#' Current version of \code{PSweight} allows for five types of propensity score weights used to estimate ATE (\code{"IPW"}), ATT {(\code{"treated"})}, and
#' ATO{(\code{"overlap"})}, ATM {\code{"matching"}} and ATEN \code{"entropy"}. These weights are members of a larger class of balancing weights defined in Li, Morgan, and Zaslavsky (2018).
#' When there is a practical violation of the positivity assumption, \code{delta} defines the symmetric
#' propensity score trimming rule following Crump et al. (2009). With multiple treatments, \code{delta} defines the
#' multinomial trimming rule introduced in Yoshida et al. (2019). The overlap weights can also be considered as
#' a data-driven continuous trimming strategy without specifying trimming rules, see Li, Thomas and Li (2019).
#' Additional details on balancing weights and generalized overlap weights for multiple treatment groups are provided in
#' Li and Li (2019). For details about matching weights and entropy weights, please refer to Li and Greene (2013) and Zhou, Matsouaka and Thomas (2020).
#'
#' @return SumStat returns a \code{SumStat} object including a list of the following value:
#' treatment group, propensity scores, propensity score weights, effective sample sizes,
#' and balance statistics. A summary of \code{SumStat} can be obtained with \code{\link{summary.SumStat}}.
#'
#' \describe{
#' \item{\code{ trtgrp}}{a character indicating the treatment group.}
#'
#' \item{\code{ propensity}}{a data frame of estimated propensity scores.}
#'
#' \item{\code{ ps.weights}}{a data frame of propensity score weights.}
#'
#' \item{\code{ ess}}{a table of effective sample sizes. This serves as a conservative measure to
#' characterize the variance inflation or precision loss due to weighting, see Li and Li (2019).}
#'
#' \item{\code{ unweighted.sumstat}}{A list of tables including covariate means and variances
#' by treatment group and standardized mean differences.}
#'
#' \item{\code{ ATE.sumstat}}{If \code{"IPW"} is included in \code{weight}, this is a list of summary statistics using inverse probability weighting.}
#'
#' \item{\code{ ATT.sumstat}}{If \code{"treated"} is included in \code{weight}, this is a list of summary statistics using the ATT weights.}
#'
#' \item{\code{ ATO.sumstat}}{If \code{"overlap"} is included in \code{weight}, this is a list of summary statistics using the overlap weights.}
#'
#' \item{\code{ ATM.sumstat}}{If \code{"matching"} is included in \code{weight}, this is a list of summary statistics using the matching weights.}
#'
#' \item{\code{ ATEN.sumstat}}{If \code{"entropy"} is included in \code{weight}, this is a list of summary statistics using the entropy weights.}
#'
#' \item{\code{ trim}}{If \code{delta > 0}, this is a table summarizing the number of observations before and after trimming.}
#'
#' }
#'
#' @references
#' Crump, R. K., Hotz, V. J., Imbens, G. W., Mitnik, O. A. (2009).
#' Dealing with limited overlap in estimation of average treatment effects. Biometrika, 96(1), 187-199.
#'
#' Li, L., Greene, T. (2013).
#' A weighting analogue to pair matching in propensity score analysis. The International Journal of Biostatistics, 9(2), 215-234.
#'
#' Li, F., Morgan, K. L., Zaslavsky, A. M. (2018).
#' Balancing covariates via propensity score weighting.
#' Journal of the American Statistical Association, 113(521), 390-400.
#'
#'
#' Li, F., Thomas, L. E., Li, F. (2019).
#' Addressing extreme propensity scores via the overlap weights. American Journal of Epidemiology, 188(1), 250-257.
#'
#' Yoshida, K., Solomon, D.H., Haneuse, S., Kim, S.C., Patorno, E., Tedeschi, S.K., Lyu, H.,
#' Franklin, J.M., Stürmer, T., Hernández-Díaz, S. and Glynn, R.J. (2019).
#' Multinomial extension of propensity score trimming methods: A simulation study.
#' American Journal of Epidemiology, 188(3), 609-616.
#'
#'
#' Li, F., Li, F. (2019). Propensity score weighting for causal inference with multiple treatments.
#' The Annals of Applied Statistics, 13(4), 2389-2415.
#'
#' Zhou, Y., Matsouaka, R. A., Thomas, L. (2020).
#' Propensity score weighting under limited overlap and model misspecification. Statistical Methods in Medical Research (Online)
#'
#'
#' @export
#'
#' @examples
#'
#' data("psdata")
#' # the propensity model
#' ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
#'
#' # using SumStat to estimate propensity scores
#' msstat <- SumStat(ps.formula, trtgrp="2", data=psdata,
#'    weight=c("IPW","overlap","treated","entropy","matching"))
#' summary(msstat)
#'
#' # importing user-supplied propensity scores "e.h"
#' fit <- nnet::multinom(formula=ps.formula, data=psdata, maxit=500, trace=FALSE)
#' e.h <- fit$fitted.values
#' varname <- c("cov1","cov2","cov3","cov4","cov5","cov6")
#' msstat0 <- SumStat(zname="trt", xname=varname, data=psdata, ps.estimate=e.h,
#'    trtgrp="2",  weight=c("IPW","overlap","treated","entropy","matching"))
#' summary(msstat0)
#'
#' @import nnet
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
#'
#'
SumStat<- function(ps.formula,ps.estimate=NULL,trtgrp=NULL,Z=NULL,covM=NULL,zname=NULL,xname=NULL,data=NULL,weight="overlap",delta=0){
  if (is.null(ps.estimate)){
    SumStat_f(ps.formula,ps.estimate=NULL,trtgrp=trtgrp,data=data,weight=weight,delta=delta)
  }else{
    SumStat_p(ps.estimate=ps.estimate,Z=Z,covM=covM,zname=zname,xname=xname,
              trtgrp=trtgrp,data=data,weight=weight,delta=delta)
  }

}

