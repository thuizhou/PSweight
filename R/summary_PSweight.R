#' Summarize a PSweight object
#'
#' \code{summary.PSweight} is used to summarize the results from \code{\link{PSweight}}.
#' The output contains the average causal effects defined by specific contrasts, as well as their
#' standard error estimates.
#'
#' @param object a PSweight object obtained from the \code{\link{PSweight}} function.
#' @param contrast a vector or matrix specifying the causal contrast of interest. The average causal effects will be
#' defined by such contrats. For multiple treatments, the contrast parameters are explained in Li and Li (2019)
#' for estimating general causal effects. Default is all pairwise contrasts between any two treatment groups.
#' @param type a character specifying the target estimand. The most commonly seen additive estimand is specified
#' by \code{type = "DIF"}, abbreviated for weighted difference-in-means. This is the usual pairwise average treatment
#' effects as those defined in Li, Morgan, and Zaslavsky (2018) and Li and Li (2019). For binary (or count outcomes), we also
#' allow two ratio estimands: causal relative risk (\code{type = "RR"}) and causal odds ratio (\code{type = "OR"}).
#' Estimates for these two ratio estimands will be reported on the log scale (log relative risk and log
#' odds ratio) to improve the approximate for asymptotic normality. With binary outcomes, \code{"DIF"} is the same
#' as the average causal risk difference. Default is "DIF" if left empty.
#' @param ... further arguments passed to or from other methods.
#'
#' @details For the \code{contrast} argument, one specifies the contrast of interest and thus defines the target estimand
#' for comparing treatments. For example, if there are three treatment levels: A, B, and C, the contrast A-C
#' (i.e., E[Y(A)] - E[Y(C)]) can be specified by \code{c(1,0,-1)}. The contrasts of A-C and B-C can be
#' jointly specified by \code{rbind(c(1,0,-1), c(0,1,-1))}.
#'
#' For estimating the causal relative risk (\code{type = "RR"}), the contrast is specified at the log scale. For example,
#' the contrast A-C (specified by \code{c(1,0,-1)}) implies the estimation of log\{E[Y(A)]\} - log\{E[Y(C)]\}. For estimating the causal odds
#' ratio, the contrast is specified at the log odds scale. For example, the contrast A-C (specified by \code{c(1,0,-1)})
#' implies the estimation of log\{E[Y(A)]/E[1-Y(A)]\} - log\{E[Y(C)]/E[1-Y(C)]\}.
#'
#' The variance of the contrasts will be estimated by the delta method (if sandwich variance is used, or
#' \code{bootstrap = FALSE}), or nonparametric bootstrap (if \code{bootstrap = TRUE}). Details will be given in
#' Zhou et al. (2020+).
#'
#' The argument \code{type} takes one of three options: \code{"DIF"}, \code{"RR"}, or \code{"RR"}, with \code{"DIF"} as
#' the default option. Typically, \code{"RR"} is relavent for binary or count outcomes, and \code{"OR"} is relavent
#' only for binary outcomes. \code{"DIF"} applies to all types of outcomes.
#'
#' @return A list of following values:
#'
#' \describe{
#' \item{\code{ trtgrp}}{ a character indicating the treatment group, or target population under ATT weights.}
#'
#' \item{\code{ estimates}}{ a matrix of point estimates, standard errors and 95% confidence intervals
#' for contrasts of interest.}
#'
#' \item{\code{ bootestimates}}{ a list of data frames containing estimated contrasts in each bootstrap replicate,
#' if bootstrap is used to estimate standard errors.}
#'
#' \item{\code{ contrast}}{ a table listing the specified contrasts of interest.}
#'
#' \item{\code{ group}}{ a table of treatment group labels corresponding to the output point estimates, provided in results
#' obtained from \code{\link{PSweight}}.}
#'
#' }
#'
#' @references
#'
#' Li, F., Morgan, K. L., Zaslavsky, A. M. (2018).
#' Balancing covariates via propensity score weighting.
#' Journal of the American Statistical Association, 113(521), 390-400.
#'
#' Li, F., Li, F. (2019). Propensity score weighting for causal inference with multiple treatments.
#' The Annals of Applied Statistics, 13(4), 2389-2415.
#'
#'
#' @export
#'
#' @examples
#'
#' ## For examples, run: example(PSweight).
#'
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
#'
summary.PSweight<-function(object,contrast=NULL,type='DIF',...){

  muhat<-object$muhat
  covmu<-object$covmu
  muboot<-object$muboot
  trtgrp<-object$trtgrp
  group<-object$group


  #groups
  ngrp<-length(group)


  #transform contrast into matrix
  if(is.vector(contrast)){
    contrast<-t(contrast)
  }

  #error message for wrong contrast
  if(!is.null(contrast)){
    if(dim(contrast)[2]!=ngrp){
      cat('Contract length not equal to treatment groups, please check','\n')
      cat('\n')
      contrast<-NULL
    }
    cat('\n')
  }




  #if all contrast
  if(is.null(contrast)){
    ncst<-ngrp*(ngrp-1)/2
    contrasttmp<-matrix(0,ngrp-1,ngrp)
    contrasttmp[,1]<--1
    contrasttmp[1:(ngrp-1),2:ngrp]<-diag(1,ngrp-1,ngrp-1)
    contrast<-contrasttmp
    if(ngrp>2){
      for(i in 1:(ngrp-2)){
        ntmp<-nrow(contrasttmp)-1
        contrasttmp<-cbind(0,contrasttmp)[1:ntmp,1:ngrp]
        contrast<-rbind(contrast,contrasttmp)
      }
    }
  }

  rownames(contrast)<-paste('Contrast',1:nrow(contrast))
  colnames(contrast)<-group


  #use bootstrap or not
  if(is.null(object$muboot)){
    if(type=='DIF'){
      est.h<-c(contrast%*%muhat)
      se.h<-sqrt(diag(contrast%*%covmu%*%t(contrast)))
      lcl<-est.h-qnorm(0.975)*se.h
      ucl<-est.h+qnorm(0.975)*se.h
    }else if(type=='RR'){
      tranmuhat<-log(muhat)
      tranest<-c(contrast%*%tranmuhat)
      trancovmu<-diag(1/muhat)%*%covmu%*%diag(1/muhat)
      transe<-sqrt(diag(contrast%*%trancovmu%*%t(contrast)))
      tranlcl<-tranest-qnorm(0.975)*transe
      tranucl<-tranest+qnorm(0.975)*transe
      est.h<-(tranest)
      se.h<-transe
      lcl<-(tranlcl)
      ucl<-(tranucl)
    }else if(type=='OR'){
      tranmuhat<-log(muhat/(1-muhat))
      tranest<-c(contrast%*%tranmuhat)
      trancovmu<-diag(1/(muhat*(1-muhat)))%*%covmu%*%diag(1/(muhat*(1-muhat)))
      transe<-sqrt(diag(contrast%*%trancovmu%*%t(contrast)))
      tranlcl<-tranest-qnorm(0.975)*transe
      tranucl<-tranest+qnorm(0.975)*transe
      est.h<-(tranest)
      se.h<-transe
      lcl<-(tranlcl)
      ucl<-(tranucl)
    }else{
      cat('type not found')
      cat('\n')
    }
  }else{
    if(type=='DIF'){
      samp<-muboot%*%t(contrast)
      est.h<-apply(samp,2,mean)
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025))
      ucl<-apply(samp,2,function(x) quantile(x,0.975))
    }else if(type=='RR'){
      samp<-log(muboot)%*%t(contrast)
      est.h<-apply(samp,2,mean)
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025))
      ucl<-apply(samp,2,function(x) quantile(x,0.975))
    }else if(type=='OR'){
      samp<-log(muboot/(1-muboot))%*%t(contrast)
      est.h<-apply(samp,2,mean)
      se.h<-sqrt(diag(cov(samp)))
      lcl<-apply(samp,2,function(x) quantile(x,0.025))
      ucl<-apply(samp,2,function(x) quantile(x,0.975))
    }else{
      cat('type not found')
      cat('\n')
    }
  }

  estimates<-cbind(est.h,se.h,lcl,ucl)
  if(is.vector(estimates)){
    t(estimates)
  }

  colnames(estimates)<-c("Estimate","Std.Error","Lower.CL","Upper.CL")
  rownames(estimates)<-rownames(contrast)
  if(is.null(object$muboot)){
    bootestimates<-NULL
  }else{
    if(is.vector(samp)){
      samp<-t(samp)
    }
    bootestimates<-samp
    colnames(bootestimates)<-rownames(contrast)
    rownames(bootestimates)<-NULL
  }
  out<-list(estimates=estimates,bootestimates=bootestimates,contrast=contrast,group=group,trtgrp=trtgrp)
  class(out)<-'PSweightsum'
  out
}
