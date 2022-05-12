#' Estimate average causal effects by propensity score weighting for a binary treatment with clustering.
#'
#' The function \code{PSweight_cl} is used to estimate the average potential outcomes corresponding to
#' each treatment group among the target population with two-level data. The function currently implements
#' the following types of weights: the inverse probability of treatment weights (IPW: target population is the combined population),
#' average treatment effect among the treated weights (treated: target population is the population receiving a specified treatment),
#' overlap weights (overlap: target population is the overlap population at clinical equipoise), matching weights (matching: target population
#' is population obtained under 1:1 matching), entropy weights (entropy: target population is the population weighted by the entropy function).
#' Augmented propensity score weighting estimators are also allowed, with propensity scores and outcome model estimated
#' within the function through mixed effect model.
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the propensity score model to be fitted. Additional details of model specification
#' are given under "Details".
#' @param trtgrp an optional character defining the "treated" population for estimating the average treatment
#' effect among the treated (ATT). Only necessary if \code{weight = "treated"}. This option can also be used to specify
#' the treatment (in a two-treatment setting) when a vector argument is supplied for \code{ps.estimate}.
#' Default value is the last group in the alphebatic order.
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
#' @param bs_level an optional character defining the cluster level (name of the variable) for each bootstrap resampling.
#' Default is \code{NULL}.
#' @param R an optional integer indicating number of bootstrap replicates. Default is \code{R = 50}.
#' @param out.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the outcome model to be fitted. Additional details of model specification
#' are given under "Details". Different from the out.formula in \code{PSweight}, this \code{\link{formula}} should include the treatment label with cooresponding cluster forms.
#' @param family a description of the error distribution and link function to be used in the outcome model.
#' Only required if \code{out.formula} is provided. Supported distributional families include
#' \code{"gaussian" (link = identity)}, \code{"binomial" (link = logit)} and \code{"poisson" (link = log)}.
#' See \code{\link{family}} in \code{\link{glm}} for more details. Default is \code{"gaussian"}.
#' @param nAGQ integer scalar - the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood.
#' Defaults to 1, corresponding to the Laplace approximation. Please refer to lme4 package for more details.
#'
#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms+1|clusters} where \code{treatment} is the treatment
#' variable  and \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment} and cluster level effects. Similarly, a typical form for \code{out.formula} is
#' \code{outcome ~ treatment+terms+1|cluster} where \code{outcome} is the outcome variable (identical to the variable name
#' used to specify \code{yname}); \code{terms} is a series of terms which specifies a linear
#' predictor for \code{outcome}; \code{clusters} is the random effects term for clusters. Both \code{ps.formula} and \code{out.formula} by default specify generalized
#' linear mixed effect models.
#'
#'
#' Current version of \code{PSweight_cl} allows for five types of propensity score weights used to estimate ATE (IPW), ATT (treated) and
#' ATO (overlap), ATM (matching) and ATEN (entropy). These weights are members of larger class of balancing weights defined in Li, Morgan, and Zaslavsky (2018).
#' Specific definitions of these weights are provided in Li, Morgan, and Zaslavsky (2018), Li and Greene (2013), Zhou, Matsouaka and Thomas (2020).
#' When there is a practical violation of the positivity assumption, \code{delta} defines the symmetric
#' propensity score trimming rule following Crump et al. (2009). The overlap weights can also be considered as
#' a data-driven continuous trimming strategy without specifying trimming rules, see Li, Thomas and Li (2019).
#' Additional details on balancing weights and generalized overlap weights for multiple treatment groups are provided in
#' Li and Li (2019).
#'
#' If \code{augmentation = TRUE}, an augmented weighting estimator will be implemented. For binary treatments, the augmented
#' weighting estimator is presented in Mao, Li and Greene (2018). When
#' \code{weight = "IPW"}, the augmented estimator is also referred to as a doubly-robust (DR) estimator.
#'
#' When \code{bootstrap = TRUE}, the variance will be calculated by nonparametric bootstrap, with \code{R} bootstrap
#' replications. \code{bs_level} needs to be specified as the variable name for the cluster in order to conduct cluster
#' level resampling and maintaining the cluster level coorelation. The default value \code{NULL} treat each observation independently.
#' The default of \code{R} is 50. Otherwise, the variance will be calculated using the sandwich variance
#' formula obtained in the M-estimation framework.
#'
#' @return PSweight_cl returns a \code{PSweight} object containing a list of the following values:
#' estimated propensity scores, average potential outcomes corresponding to each treatment,
#' variance-covariance matrix of the point estimates, the label for each treatment group,
#' and estimates in each bootstrap replicate if \code{bootstrap = TRUE}.
#' A summary of PSweight_cl can be obtained with \code{\link{summary.PSweight}}.
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
#' Li, F., Zaslavsky, A. M., Landrum, M. B. (2013).
#' Propensity score weighting with multilevel data. Statistics in Medicine, 32(19), 3373-3387.
#'
#' Fuentes, A., Lüdtke, O., Robitzsch, A. (2021).
#' Causal inference with multilevel data: A comparison of different propensit score weighting appropaches. Multivariate Behavioral Research, 1-24.
#'
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
#'
#' Li, F., Li, F. (2019). Propensity score weighting for causal inference with multiple treatments.
#' The Annals of Applied Statistics, 13(4), 2389-2415.
#'
#' Zhou, Y., Matsouaka, R. A., Thomas, L. (2020).
#' Propensity score weighting under limited overlap and model misspecification. Statistical Methods in Medical Research 29(12), 3721-3756.
#'
#'
#'
#' @export
#'
#' @examples
#' #data("psdata_cl")
#' #ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6+(1|clt)
#' #ato_cl<-PSweight(ps.formula = ps.formula,yname = 'Y',data = psdata_cl)
#' #summary(ato_cl)
#'
#'
#' @import lme4
#' @import nnet
#' @import MASS
#' @import numDeriv
#' @importFrom  stats binomial coef cov formula glm lm model.matrix model.extract model.frame plogis poisson predict qnorm quantile sd as.formula printCoefmat
#' @importFrom  utils capture.output combn tail
#' @importFrom  graphics hist legend
#'
PSweight_cl<-function(ps.formula=NULL,trtgrp=NULL,yname,data,weight='overlap',delta=0,augmentation=FALSE,bootstrap=FALSE,bs_level=NULL,R=50,out.formula=NULL,family='gaussian',nAGQ=1L){

  #extract zname
  ps.formula<-as.formula(ps.formula)
  zname<-all.vars(ps.formula)[1]


  #trim the data
  if(delta>0){
    fit.e<-lme4::glmer(formula = ps.formula, data=data, nAGQ=nAGQ,family="binomial")
    e.h <- as.numeric(predict(fit.e,newdata=data,type="response"))
    ps.estimate<-cbind(1-e.h,e.h)

    trimobj<-do.call(PStrim,list(data=data,ps.estimate=ps.estimate,zname=zname,delta=delta))
    data<-trimobj$data
  }


  data[zname]<-as.character(unlist(data[zname]))
  if (!is.null(bs_level)){
    data[,bs_level]<-as.character(data[,bs_level])
  }
  do.call(binest_cl,list(ps.formula=ps.formula,zname=zname,yname=yname,data=data,trtgrp=trtgrp,augmentation=augmentation,bootstrap=bootstrap,bs_level=bs_level,R=R,out.formula=out.formula,family=family,weight=weight,nAGQ=nAGQ))
}



