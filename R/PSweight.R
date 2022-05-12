#' Estimate average causal effects by propensity score weighting
#'
#' The function \code{PSweight} is used to estimate the average potential outcomes corresponding to
#' each treatment group among the target population. The function currently implements
#' the following types of weights: the inverse probability of treatment weights (IPW: target population is the combined population),
#' average treatment effect among the treated weights (treated: target population is the population receiving a specified treatment),
#' overlap weights (overlap: target population is the overlap population at clinical equipoise), matching weights (matching: target population
#' is population obtained under 1:1 matching), entropy weights (entropy: target population is the population weighted by the entropy function).
#' Augmented propensity score weighting estimators are also allowed, with propensity scores and outcome model estimates either estimated
#' within the function, or supplied by external routines.
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the propensity score model to be fitted. Additional details of model specification
#' are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param ps.estimate an optional matrix or data frame containing estimated (generalized) propensity scores for
#' each observation. Typically, this is an N by J matrix, where N is the number of observations and J is the
#' total number of treatment levels. Preferably, the column name of this matrix should match the name of treatment level,
#' if column name is missing or there is a mismatch, the column names would be assigned according to alphabatic order
#' of the treatment levels. A vector of propensity score estimates is also allowed in \code{ps.estimate}, in which
#' case a binary treatment is implied and the input is regarded as the propensity to receive the last category of
#' treatment by alphabatic order, unless otherwise stated by \code{trtgrp}.
#' @param trtgrp an optional character defining the "treated" population for estimating the average treatment
#' effect among the treated (ATT). Only necessary if \code{weight = "treated"}. This option can also be used to specify
#' the treatment (in a two-treatment setting) when a vector argument is supplied for \code{ps.estimate}.
#' Default value is the last group in the alphebatic order.
#' @param zname an optional character specifying the name of the treatment variable in \code{data}.
#' @param yname an optional character specifying the name of the outcome variable in \code{data}.
#' @param data an optional data frame containing the variables in the propensity score model
#' and outcome model (if augmented estimator is used). If not found in data, the variables are
#' taken from \code{environment(formula)}.
#' @param weight a character or vector of characters including the types of weights to be used.
#' \code{"IPW"} specifies the inverse probability of treatment weights for estimating the average treatment
#' effect among the combined population. \code{"treated"} specifies the weights for estimating the
#' average treatment effect among the treated. \code{"overlap"} specifies the (generalized) overlap weights
#' for estimating the average treatment effect among the overlap population, or population at
#' clinical equipoise. \code{"matching"} specifies the matching weights for estimating the average treatment effect
#' among the matched population (ATM). \code{"entropy"} specifies the entropy weights for the average treatment effect
#' of entropy weighted population (ATEN). Default is \code{"overlap"}.
#' @param delta trimming threshold for estimated (generalized) propensity scores.
#' Should be no larger than 1 / number of treatment groups. Default is 0, corresponding to no trimming.
#' @param augmentation logical. Indicate whether augmented weighting estimators should be used.
#' Default is \code{FALSE}.
#' @param bootstrap logical. Indaicate whether bootstrap is used to estimate the standard error
#' of the point estimates. Default is \code{FALSE}.
#' @param R an optional integer indicating number of bootstrap replicates. Default is \code{R = 50}.
#' @param out.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the outcome model to be fitted. Additional details of model specification
#' are given under "Details". This argument is optional if \code{out.estimate} is not \code{NULL}.
#' @param out.estimate an optional matrix or data frame containing estimated potential outcomes
#' for each observation. Typically, this is an N by J matrix, where N is the number of observations
#' and J is the total number of treatment levels. Preferably, the column name of this matrix should
#' match the name of treatment level, if column name is missing or there is a mismatch,
#' the column names would be assigned according to alphabatic order of the treatment levels, with a
#' similar mechanism as in \code{ps.estimate}.
#' @param family a description of the error distribution and link function to be used in the outcome model.
#' Only required if \code{out.formula} is provided. Supported distributional families include
#' \code{"gaussian" (link = identity)}, \code{"binomial" (link = logit)} and \code{"poisson" (link = log)}.
#' See \code{\link{family}} in \code{\link{glm}} for more details. Default is \code{"gaussian"}.
#' @param ps.method a character to specify the method for estimating propensity scores. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param ps.control a list to specify additional options when \code{method} is set to \code{"gbm"} or \code{"SuperLearner"}.
#' @param out.method a character to specify the method for estimating the outcome regression model. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param out.control a list to specify additional options when \code{out.method} is set to \code{"gbm"} or \code{"SuperLearner"}.
#'
#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable (identical to the variable name used to specify \code{zname}) and \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment}. Similarly, a typical form for \code{out.formula} is
#' \code{outcome ~ terms} where \code{outcome} is the outcome variable (identical to the variable name
#' used to specify \code{yname}) and \code{terms} is a series of terms which specifies a linear
#' predictor for \code{outcome}. Both \code{ps.formula} and \code{out.formula} by default specify generalized
#' linear models when \code{ps.estimate} and/or \code{out.estimate} is \code{NULL}. The argument \code{ps.method} and \code{out.method} allow user to choose
#' model other than glm to fit the propensity score and outcome regression models for augmentation. Additional argument in the \code{gbm()} function can be supplied through the \code{ps.control} and \code{out.control} argument. Please refer to the user manual of the \code{gbm} package for all the
#' allowed arguments.  \code{"SuperLearner"} is also allowed in the \code{ps.method} and \code{out.method} arguments. Currently, the SuperLearner method only supports binary treatment with the default method set to \code{"SL.glm"}. The estimation approach is default to \code{"method.NNLS"} for both propensity and outcome regression models.
#' Prediction algorithm and other tuning parameters can also be passed through\code{ps.control} and \code{out.control}. Please refer to the user manual of the \code{SuperLearner} package for all the allowed specifications.
#'
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
#' the alphebatic order of values in \code{Z}, and the rightmost coulmn of \code{ps.estimate} is assumed
#' to be the treatment group, when estimating ATT. \code{trtgrp} can also be used to specify the treatment
#' group for estimating ATT. The same mechanism applies to \code{out.estimate}, except that the input for \code{out.estimate}
#' must be an N by J matrix, where each row corresponds to the estimated potential outcomes (corresponding to each treatment)
#' for each observation.
#'
#' The argument \code{zname} and/or \code{yname} is required when \code{ps.estimate}
#' and/or \code{out.estimate} is not \code{NULL}.
#'
#' Current version of \code{PSweight} allows for five types of propensity score weights used to estimate ATE (IPW), ATT (treated) and
#' ATO (overlap), ATM (matching) and ATEN (entropy). These weights are members of larger class of balancing weights defined in Li, Morgan, and Zaslavsky (2018).
#' Specific definitions of these weights are provided in Li, Morgan, and Zaslavsky (2018), Li and Greene (2013), Zhou, Matsouaka and Thomas (2020).
#' When there is a practical violation of the positivity assumption, \code{delta} defines the symmetric
#' propensity score trimming rule following Crump et al. (2009). With multiple treatments, \code{delta} defines the
#' multinomial trimming rule introduced in Yoshida et al. (2019). The overlap weights can also be considered as
#' a data-driven continuous trimming strategy without specifying trimming rules, see Li, Thomas and Li (2019).
#' Additional details on balancing weights and generalized overlap weights for multiple treatment groups are provided in
#' Li and Li (2019).
#'
#' If \code{augmentation = TRUE}, an augmented weighting estimator will be implemented. For binary treatments, the augmented
#' weighting estimator is presented in Mao, Li and Greene (2018). For multiple treatments, the augmented weighting estimator is
#' mentioned in Li and Li (2019), and additional details will appear in our ongoing work (Zhou et al. 2020+). When
#' \code{weight = "IPW"}, the augmented estimator is also referred to as a doubly-robust (DR) estimator.
#'
#' When \code{bootstrap = TRUE}, the variance will be calculated by nonparametric bootstrap, with \code{R} bootstrap
#' replications. The default of \code{R} is 50. Otherwise, the variance will be calculated using the sandwich variance
#' formula obtained in the M-estimation framework.
#'
#' @return PSweight returns a \code{PSweight} object containing a list of the following values:
#' estimated propensity scores, average potential outcomes corresponding to each treatment,
#' variance-covariance matrix of the point estimates, the label for each treatment group,
#' and estimates in each bootstrap replicate if \code{bootstrap = TRUE}.
#' A summary of PSweight can be obtained with \code{\link{summary.PSweight}}.
#'
#' \describe{
#'
#' \item{\code{ trtgrp}}{a character indicating the treatment group.}
#'
#' \item{\code{ propensity}}{a data frame of estimated propensity scores.}
#'
#' \item{\code{ muhat}}{ average potential outcomes by treatment groups, with reference to specific target populations.}
#'
#' \item{\code{ covmu}}{ variance-covariance matrix of \code{muhat}.}
#'
#' \item{\code{ muboot}}{ an optional list of point estimates in each bootstrap replicate \code{bootstrap = TRUE}.}
#'
#' \item{\code{ group}}{ a table of treatment group labels corresponding to the output point estimates \code{muhat}.}
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
#' Mao, H., Li, L., Greene, T. (2019). Propensity score weighting analysis and treatment effect discovery.
#' Statistical Methods in Medical Research, 28(8), 2439-2454.
#'
#' Li, F., Thomas, L. E., Li, F. (2019).
#' Addressing extreme propensity scores via the overlap weights. American Journal of Epidemiology, 188(1), 250-257.
#'
#' Yoshida, K., Solomon, D.H., Haneuse, S., Kim, S.C., Patorno, E., Tedeschi, S.K., Lyu, H.,
#' Franklin, J.M., Stürmer, T., Hernández-Díaz, S. and Glynn, R.J. (2019).
#' Multinomial extension of propensity score trimming methods: A simulation study.
#' American Journal of Epidemiology, 188(3), 609-616.
#'
#' Li, F., Li, F. (2019). Propensity score weighting for causal inference with multiple treatments.
#' The Annals of Applied Statistics, 13(4), 2389-2415.
#'
#' Zhou, Y., Matsouaka, R. A., Thomas, L. (2020).
#' Propensity score weighting under limited overlap and model misspecification. Statistical Methods in Medical Research, 29(12), 3721-3756.
#'
#'
#'
#' @export
#'
#' @examples
#' data("psdata")
#' # the propensity and outcome models
#' ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
#' out.formula<-Y~cov1+cov2+cov3+cov4+cov5+cov6
#'
#' # without augmentation
#' ato1<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,weight = 'overlap')
#' summary(ato1)
#'
#' # augmented weighting estimator, takes longer time to calculate sandwich variance
#' # ato2<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata,
#' #              augmentation = TRUE,out.formula = out.formula,family = 'gaussian',weight = 'overlap')
#' # summary(ato2)
#'
#' @import nnet
#' @import MASS
#' @import numDeriv
#' @importFrom  stats binomial coef cov formula glm lm model.matrix model.extract model.frame plogis poisson predict qnorm quantile sd as.formula printCoefmat
#' @importFrom  utils capture.output combn tail
#' @importFrom  graphics hist legend
#'
PSweight<-function(ps.formula=NULL,ps.estimate=NULL,trtgrp=NULL,zname=NULL,yname,data,weight='overlap',delta=0,augmentation=FALSE,bootstrap=FALSE,R=50,out.formula=NULL,out.estimate=NULL,family='gaussian',ps.method='glm',ps.control=list(),out.method='glm',out.control=list()){

  #extract zname
  if(!is.null(ps.formula)){
    ps.formula<-as.formula(ps.formula)
    zname<-all.vars(ps.formula)[1]
  }

  data[zname]<-as.character(unlist(data[zname]))
  categoryz1<-unique(unlist(data[zname]))
  z1<-as.numeric(factor(unlist(data[zname])))
  oldlevel1<-categoryz1[order(unique(z1))]
  ncate<-length(categoryz1)

  #trim the data
  if(delta>0){
    trimobj<-do.call(PStrim,list(data=data,ps.formula = ps.formula, zname=zname, ps.estimate=ps.estimate,delta=delta,optimal=FALSE,out.estimate=out.estimate,method=ps.method,ps.control=ps.control))
    data<-trimobj$data
    ps.estimate<-trimobj$ps.estimate
    out.estimate<-trimobj$out.estimate
  }


  if(ncate==2){

    do.call(binest,list(ps.formula=ps.formula,ps.estimate=ps.estimate,zname=zname,yname=yname,data=data,trtgrp=trtgrp,augmentation=augmentation,bootstrap=bootstrap,R=R,out.formula=out.formula,out.estimate=out.estimate,family=family,weight=weight,ps.method=ps.method,ps.control=ps.control,out.method=out.method,out.control=out.control))

  }else{

    do.call(mulest,list(ps.formula=ps.formula,ps.estimate=ps.estimate,zname=zname,yname=yname,data=data,trtgrp=trtgrp,augmentation=augmentation,bootstrap=bootstrap,R=R,out.formula=out.formula,out.estimate=out.estimate,family=family,weight=weight,ps.method=ps.method,ps.control=ps.control,out.method=out.method,out.control=out.control))
  }
}


