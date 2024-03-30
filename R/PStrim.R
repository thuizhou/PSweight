#' Trim the input data and propensity estimate
#'
#' Trim the original data and propensity estimate according to symmetric propensity score trimming rules.
#'
#' @param data an optional data frame containing the variables required by \code{ps.formula}.
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the propensity score model to be fitted. Additional details of model specification are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param ps.estimate an optional matrix or data frame containing estimated (generalized) propensity scores for each observation. Typically, this is an N by J matrix, where N is the number of observations and J is the total number of treatment levels. Preferably, the column name of this matrix should match the name of treatment level, if column name is missing or there is a mismatch, the column names would be assigned according to alphabatic order of the treatment levels. A vector of propensity score estimates is also allowed in \code{ps.estimate}, in which case a binary treatment is implied and the input is regarded as the propensity to receive the last category of treatment by alphabatic order, unless otherwise stated by \code{trtgrp}.
#' @param zname an optional character specifying the name of the treatment variable in \code{data}. Unless \code{ps.formula} is specified, \code{zname} is required.
#' @param delta trimming threshold for estimated (generalized) propensity scores. Should be no larger than 1 / number of treatment groups. Default is 0, corresponding to no trimming.
#' @param optimal an logical argument indicating if optimal trimming should be used. Default is \code{FALSE}.
#' @param out.estimate an optional matrix or data frame containing estimated potential outcomes
#' for each observation. Typically, this is an N by J matrix, where N is the number of observations
#' and J is the total number of treatment levels. Preferably, the column name of this matrix should
#' match the name of treatment level, if column name is missing or there is a mismatch,
#' the column names would be assigned according to alphabatic order of the treatment levels, with a
#' similar mechanism as in \code{ps.estimate}.
#' @param method a character to specify the method for estimating propensity scores. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param ps.control a list to specify additional options when \code{method} is set to \code{"gbm"} or \code{"SuperLearner"}.
#'
#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable (identical to the variable name used to specify \code{zname}) and \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment}. \code{ps.formula} specifies a
#' model for estimating the propensity scores, when \code{ps.estimate} is \code{NULL}.
#'  \code{"glm"} is the default method for propensity score estimation. Logistic regression will be used for binary outcomes,
#' and multinomial logistic regression will be used for outcomes with more than two categories. The alternative method option of \code{"gbm"} serves as an API to call the \code{gbm()} function from the
#' \code{gbm} package. Additional argument in the \code{gbm()} function can be supplied through the \code{ps.control=list()} argument in \code{SumStat()}. Please refer to the user manual of the \code{"gbm"} package for all the
#' allowed arguments. Currently, models for binary or multinomial treatment will be automatically chosen based on the number of treatment categories.
#' \code{"SuperLearner"} is also allowed in the \code{method} argument to call the \code{SuperLearner()} function in \code{SuperLearner} package.
#' Currently, the SuperLearner method only support binary treatment with the default method set to \code{"SL.glm"}. The estimation approach is default to \code{"method.NNLS"}.
#' Prediction algorithm and other tuning parameters can also be passed through \code{ps.control=list()}. Please refer to the user manual of the \code{SuperLearner} package for all the allowed specifications.
#'
#'
#' When comparing two treatments, \code{ps.estimate} can either be a vector or a two-column matrix of estimated
#' propensity scores. If a vector is supplied, it is assumed to be the propensity scores to receive the treatment, and
#' the treatment group corresponds to the last group in the alphebatic order, unless otherwise specified by \code{trtgrp}.
#' When comparing multiple (J>=3) treatments, \code{ps.estimate} needs to be specified as an N by J matrix,
#' where N indicates the number of observations, and J indicates the total number of treatments.
#' This matrix specifies the estimated generalized propensity scores to receive each of the J treatments.
#' The same mechanism applies to \code{out.estimate}, except that the input for \code{out.estimate}
#' must be an N by J matrix, where each row corresponds to the estimated potential outcomes (corresponding to each treatment)
#' for each observation.
#'
#' With binary treatments, \code{delta} defines the symmetric propensity score trimming rule following Crump et al. (2009).
#' With multiple treatments, \code{delta} defines the symmetric multinomial trimming rule introduced in Yoshida et al. (2019).
#' With binary treatments and when \code{optimal} equals \code{TRUE}, the trimming function implements the optimal
#' symmetric trimming rule in Crump et al. (2009). The optimal trimming threshold \code{delta} is then returned.
#' With multiple treatments and \code{optimal} equals \code{TRUE}, the trimming function implements the optimal trimming rule in Yang et al. (2016).
#' The optimal cutoff \code{lambda}, which defines the acceptable upper bound for the sum of inverse generalized propensity scores, is
#' returned. See Yang et al. (2016) and Li and Li (2019) for details.
#'
#' The argument \code{zname} is required when \code{ps.estimate} is not \code{NULL}.
#'
#' @return PStrim returns a list of the following values:
#' \describe{
#'
#' \item{\code{ data}}{a data frame of trimmed data. }
#'
#' \item{\code{ trim_sum}}{a table summrizing the number of cases by treatment groups before and after trimming. }
#'
#' \item{\code{ ps.estimate}}{ a data frame of propensity estimate after trimming. }
#'
#' \item{\code{ delta}}{ an optional output of trimming threshold for symmetric trimming. }
#'
#' \item{\code{ lambda}}{ an optional output trimming threshold for optimal trimming with multiple treatment groups. }
#'
#' \item{\code{ out.estimate}}{ a data frame of estimated potential outcomes after trimming. }
#' }
#'
#'
#' @references
#' Crump, R. K., Hotz, V. J., Imbens, G. W., Mitnik, O. A. (2009).
#' Dealing with limited overlap in estimation of average treatment effects. Biometrika, 96(1), 187-199.
#'
#' Yoshida, K., Solomon, D.H., Haneuse, S., Kim, S.C., Patorno, E., Tedeschi, S.K., Lyu, H.,
#' Franklin, J.M., Stürmer, T., Hernández-Díaz, S. and Glynn, R.J. (2019).
#' Multinomial extension of propensity score trimming methods: A simulation study.
#' American Journal of Epidemiology, 188(3), 609-616.
#'
#' Yang, S., Imbens, G. W., Cui, Z., Faries, D. E., Kadziola, Z. (2016). Propensity score matching
#' and subclassification in observational studies with multi-level treatments. Biometrics, 72(4), 1055-1065.
#'
#' Li, F., Li, F. (2019). Propensity score weighting for causal inference with multiple treatments.
#' The Annals of Applied Statistics, 13(4), 2389-2415.
#'
#'
#' @export
#'
#' @examples
#' data("psdata")
#'
#' # the propensity model
#' ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6
#'
#' # trim the original data by setting the threshold of propensity as 0.05
#' PStrim(data=psdata, ps.formula=ps.formula, delta=0.05)
#' PStrim(data=psdata, ps.formula=ps.formula, optimal=TRUE)
#'
#' @import nnet
#' @import MASS
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd uniroot
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend



PStrim<-function(data,ps.formula=NULL,zname=NULL,ps.estimate=NULL,delta=0,optimal=FALSE,out.estimate=NULL,method='glm',ps.control=list()){

  #extract zname
  if(!is.null(ps.formula)){
    ps.formula<-as.formula(ps.formula)
    zname<-all.vars(ps.formula)[1]
  }


  datatmp<-data
  z<-as.factor(unlist(data[zname]))
  datatmp[zname]<-z
  ncate<-length(levels(z))
  n<-length(z)

  if(!optimal){
    if(is.null(ps.estimate)){
      #number of groups
      if(1/ncate<=delta){
        warning('invalid trimming, return original data')
      }else{
        propensity<-do.call(PSmethod,list(ps.formula = ps.formula, method=method, data=datatmp, ncate=ncate, ps.control=ps.control))$e.h
        psidx<-apply(as.matrix(propensity),1,function(x) min(x)>delta)
        if(length(unique(z[psidx]))==ncate){
          data<-data[psidx,]
          out.estimate<-out.estimate[psidx,]
        }else{
          warning('One or more groups removed after trimming, reset your delta, trimming not applied')
        }
      }
    }else{
      if(1/ncate<=delta){
        warning('invalid trimming, return original data')
      }else{
        if(length(ps.estimate)==n){
          ps.estimate<-cbind(1-ps.estimate,ps.estimate)
        }

        psidx<-apply(as.matrix(ps.estimate),1,function(x) min(x)>delta)

        if(length(unique(z[psidx]))==ncate){
          data<-data[psidx,]
          ps.estimate<-ps.estimate[psidx,]
          out.estimate<-out.estimate[psidx,]
        }
      }
    }
  }

  if(optimal){
    if (delta>0) warning('delta is ignored when optimal trimming is used')
    if (ncate==2){
      if(is.null(ps.estimate)){
        propensity<-do.call(PSmethod,list(ps.formula = ps.formula, method=method, data=datatmp,ncate=ncate,ps.control=ps.control))$e.h
      }else{
        if(length(ps.estimate)==n){
          ps.estimate<-cbind(1-ps.estimate,ps.estimate)
        }
        propensity<-ps.estimate
      }

      # find optimal trimming rule
      k <- 2*mean(1/(propensity[,1]*propensity[,2]))-max(1/(propensity[,1]*propensity[,2]))
      if(k >= 0){
        tolg <- 0
      } else {
        sum.wt <- as.numeric(1/(propensity[,1]*propensity[,2]))
        trim.fun <- function(x){
          sum.wt.trim <- sum.wt[sum.wt <= x]
          return(x - 2*mean(sum.wt.trim))
        }
        trim.fun <- Vectorize(trim.fun)
        ran<-range(sum.wt)
        lambda <- uniroot(trim.fun, lower = ran[1], upper= ran[2])$root
        tolg <- 1/2 - sqrt(1/4 - 1/lambda)
      }
    }else{
      if(is.null(ps.estimate)){
        propensity<-do.call(PSmethod,list(ps.formula = ps.formula, method=method, data=datatmp,ncate=ncate,ps.control=ps.control))$e.h
      }else{
        propensity<-ps.estimate
      }

      # find optimal trimming rule
      sum.wt <- as.numeric(rowSums(1/propensity))
      trim.fun <- function(x){
        sum.wt.trim <- sum.wt[sum.wt <= x]
        return(x - 2*mean(sum.wt.trim) / mean(sum.wt <= x))
      }
      trim.fun <- Vectorize(trim.fun)
      if(trim.fun(max(rowSums(1/propensity))) < 0){
        lambda <- max(rowSums(1/propensity))+1
      } else{
        ran<-range(rowSums(1/propensity))
        lambda <- uniroot(trim.fun, lower=ran[1], upper=ran[2])$root
      }
      keep <- (sum.wt <= lambda)
    }

    if (ncate==2){
      delta<-tolg
      psidx<-apply(as.matrix(propensity),1,function(x) min(x)>delta)
    }else{ psidx<-keep }

    if(length(unique(z[psidx]))==ncate){
      data<-data[psidx,]
    }else{
      warning('One or more groups removed after trimming, reset your delta, trimming not applied')
    }


    ps.estimate<-ps.estimate[psidx,]
    out.estimate<-out.estimate[psidx,]
  }

  trimmed<-table(datatmp[zname])-table(data[zname])
  remained<-table(data[zname])
  trim_sum<-rbind(trimmed,remained)

  if (ncate>2 & optimal){
    out<-list(data=data,trim_sum=trim_sum,ps.estimate=ps.estimate,delta=NULL,lambda=lambda, out.estimate=out.estimate)
  }
  else{
    out<-list(data=data,trim_sum=trim_sum,ps.estimate=ps.estimate,delta=delta,lambda=NULL, out.estimate=out.estimate)
  }

  class(out)<-'PStrim'
  out

}



