#' Trim the input data and propensity estimate
#'
#' Trim the original data and propensity estimate according to symmetric propensity score trimming rules.
#'
#' @param data an optional data frame containing the variables required by \code{ps.formula}.
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the propensity score model to be fitted. Additional details of model specification are given under ‘Details’. This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param ps.estimate an optional matrix or data frame containing estimated (generalized) propensity scores for each observation. Typically, this is an N by J matrix, where N is the number of observations and J is the total number of treatment levels. Preferably, the column name of this matrix should match the name of treatment level, if column name is missing or there is a mismatch, the column names would be assigned according to alphabatic order of the treatment levels. A vector of propensity score estimates is also allowed in \code{ps.estimate}, in which case a binary treatment is implied and the input is regarded as the propensity to receive the last category of treatment by alphabatic order, unless otherwise stated by \code{trtgrp}.
#' @param zname an optional character specifying the name of the treatment variable in \code{data}.
#' @param delta trimming threshold for estimated (generalized) propensity scores. Should be no larger than 1 / number of treatment groups. Default is 0, corresponding to no trimming.
#'
#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable (identical to the variable name used to specify \code{zname}) and \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment}. \code{ps.formula} specifies generalized linear
#' model for estimating the propensity scores, when \code{ps.estimate} is \code{NULL}. See \code{glm} for
#' more details on generalized linear models.
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
#' group for estimating ATT.
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
#' \item{\code{ ps.estimate}}{ an optional dataframe of propensity estimate after trimming if propensity estimate is imported. }
#' }
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
#'
#' @import nnet
#' @import MASS
#' @import numDeriv
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
#'



PStrim<-function(data,ps.formula=NULL,zname=NULL,ps.estimate=NULL,delta=0){

  #extract zname
  if(typeof(ps.formula)!="character"){
    ps.formula<-Reduce(paste0,deparse(ps.formula))
  }
  ps.formula<-gsub(" ","",ps.formula)

  if(is.null(zname)){
    zname<-unlist(strsplit(ps.formula,'~'))[[1]][1]
  }


  datatmp<-data
  z<-as.factor(unlist(data[zname]))
  datatmp[zname]<-z
  ncate<-length(unique(z))

  if(is.null(ps.estimate)){
    #number of groups
    if(1/ncate<=delta){
      warning('invalid trimming, return original data')
    }else{
      if(ncate==2){
        fittrim <- glm(ps.formula, family = binomial(link = "logit"),data=datatmp)
        propensity<-cbind(1-fittrim$fitted.values,fittrim$fitted.values)
      }else{
        fittrim<- nnet::multinom(formula = ps.formula, data=datatmp,maxit = 500, Hess = TRUE, trace = FALSE)
        propensity<-fittrim$fitted.values
      }
      psidx<-apply(as.matrix(propensity),1,function(x) min(x)>delta)
      if(length(unique(z[psidx]))==ncate){
        data<-data[psidx,]
      }else{
        warning('One or more groups removed after trimming, reset your delta, trimming not applied')
      }
    }
  }else{
    if(1/ncate<=delta){
      warning('invalid trimming, return original data')
    }else{
      if(is.vector(ps.estimate)){
        ps.estimate<-cbind(1-ps.estimate,ps.estimate)
      }else if(ncol(ps.estimate)==1){
        ps.estimate<-cbind(1-ps.estimate,ps.estimate)
      }
      psidx<-apply(as.matrix(ps.estimate),1,function(x) min(x)>delta)

      if(length(unique(z[psidx]))==ncate){
        data<-data[psidx,]
        ps.estimate<-ps.estimate[psidx,]
      }else{
        warning('One or more groups removed after trimming, reset your delta, trimming not applied')
      }
    }
  }

  trimmed<-table(datatmp[zname])-table(data[zname])
  remained<-table(data[zname])
  trim_sum<-rbind(trimmed,remained)
  return(list(data=data,trim_sum=trim_sum,ps.estimate=ps.estimate))
}



