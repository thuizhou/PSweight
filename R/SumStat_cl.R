#' Calculate summary statistics for propensity score weighting with clustering (for binary treatment only)
#'
#' \code{SumStat_cl} is used to generate distributional plots of the estimated propensity scores and balance
#' diagnostics after propensity score weighting with two-level data.
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the propensity score model to be fitted. Additional details of model specification are given under "Details".
#' @param trtgrp an optional character defining the "treated" population for estimating the average treatment effect among the treated (ATT). Only necessary if \code{weight = "treated"}. This option can also be used to specify the treatment (in a two-treatment setting). Default value is the last group in the alphebatic order.
#' @param data an data frame containing the variables in the propensity score model. If not found in data, the variables are taken from \code{environment(formula)}.
#' @param weight a character or vector of characters including the types of weights to be used. \code{"IPW"} specifies the inverse probability weights for estimating the average treatment effect among the combined population (ATE). \code{"treated"} specifies the weights for estimating the average treatment effect among the treated (ATT). \code{"overlap"} specifies the (generalized) overlap weights for estimating the average treatment effect among the overlap population (ATO), or population at clinical equipoise. \code{"matching"} specifies the matching weights for estimating the average treatment effect among the matched population (ATM). \code{"entropy"} specifies the entropy weights for the average treatment effect of entropy weighted population (ATEN). Default is \code{"overlap"}.
#' @param delta trimming threshold for estimated (generalized) propensity scores. Should be no larger than 1 / number of treatment groups. Default is 0, corresponding to no trimming.
#' @param nAGQ integer scalar - the number of points per axis for evaluating the adaptive Gauss-Hermite approximation to the log-likelihood.
#' Defaults to 1, corresponding to the Laplace approximation. Please refer to lme4 package for more details.
#' @details A typical form for \code{ps.formula} is \code{treatment ~ terms+1|clusters} where \code{treatment} is the treatment
#' variable, \code{terms} is a series of terms
#' which specifies a linear predictor for \code{treatment}, and \code{clusters} is the cluster indicator. The current version supports two-level models and the random-effects term is required to be the last piece in the formula. \code{ps.formula} specifies a mixed-effects logistic regression
#' model for estimating propensity scores. The treatment group corresponds to the last group in the alphebatic order, unless otherwise specified by \code{trtgrp}.
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
#' @return SumStat_cl returns a \code{SumStat} object including a list of the following value:
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
#' Li, F., Thomas, L. E., Li, F. (2019).
#' Addressing extreme propensity scores via the overlap weights. American Journal of Epidemiology, 188(1), 250-257.
#'
#' Li, F., Li, F. (2019). Propensity score weighting for causal inference with multiple treatments.
#' The Annals of Applied Statistics, 13(4), 2389-2415.
#'
#' Zhou, Y., Matsouaka, R. A., Thomas, L. (2020).
#' Propensity score weighting under limited overlap and model misspecification. Statistical Methods in Medical Research 29(12), 3721-3756.
#'
#' Li, F., Zaslavsky, A. M., & Landrum, M. B. (2013). Propensity score weighting with multilevel data. Statistics in Medicine, 32(19), 3373-3387.
#'
#' @export
#'
#' @examples
#'
#' data("psdata_cl")
#' # the propensity model
#' # ps.formula<-trt~cov1+cov2+cov3+cov4+cov5+cov6+(1|clt)
#'
#' # using SumStat to estimate propensity scores
#' # msstat <- SumStat_cl(ps.formula, trtgrp="1", data=psdata_cl,
#' #   weight=c("IPW","overlap","treated","entropy","matching"))
#' #summary(msstat)
#'
#' @import nnet SuperLearner gbm lme4
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
#'
#'
SumStat_cl<- function(ps.formula=NULL,trtgrp=NULL,data=NULL,weight="overlap",delta=0,nAGQ=1L){

    #extract z name
    zname<-all.vars(ps.formula)[1]
    xname<-all.vars(ps.formula)[-c(1,length(all.vars(ps.formula)))]

    #set ordered group
    facz<-as.factor(unlist(data[,zname]))
    data[,zname]<-facz

    #treatment label
    #creat a dictionary for the original and recoded values in Z
    dic<-levels(facz)
    ncate<-length(dic) #number of categories
    if (!is.null(trtgrp)) trt<-which(dic==trtgrp)

    if (ncate==1) stop("Treatment variable needs to have more than 1 category.","\n")


    #trim the data
    if(delta>0){
      fit.e<-lme4::glmer(formula = ps.formula, data=data, nAGQ=nAGQ,family="binomial")
      e.h <- as.numeric(predict(fit.e,newdata=data,type="response"))
      ps.estimate<-cbind(1-e.h,e.h)
      did<-which(ps.estimate[,1]<delta | ps.estimate[,1]> (1-delta))
      datatmp<-data
      data<-data[-did,]
      trimmed<-table(datatmp[zname])-table(data[zname])
      remained<-table(data[zname])
      trim_sum<-rbind(trimmed,remained)
    }

    #fit propensity score model
    fit.e <-lme4::glmer(formula = ps.formula, data=data, nAGQ=nAGQ,family="binomial")
    e.h <- as.numeric(predict(fit.e,newdata=data,type="response"))
    e.h <- cbind((1-e.h),e.h)
    #post-trimming processing
    z<-as.numeric(data[,zname])
    n<-length(z) #total obs
    data["zindex"]<-z
    covM<-data[,xname]

    #transform factor/categorical column, if any
    cat_index<-c()
    for(i in 1:ncol(covM)){
      if(inherits(covM[,i],"factor")) {
        if(length(table(covM[,i]))==1){
          stop("invariate covariate supplied ","\n")
        }else{
          covmtmp<-covM[,i]
          lev<-unique(as.character(as.factor(covmtmp)))
          cols <- paste0(names(covM)[i], lev)
          covM[cols] <-0
          ncols=length(cols)
          for(j in 1:ncols){
            covM[cols[j]]<-1.*(covmtmp==lev[j])
          }
          cat_index<-c(cat_index,i)
        }
        #stop("Non-numeric values exist in covariates, please transform","\n")
      }
    }

    if(!is.null(cat_index)){
      covM<-covM[,-cat_index]
    }
    covM<-as.matrix(covM)


    #dimension of data
    n<-dim(data)[1]
    ncate<-length(unique(data$trt))

    #correspondence btw Z and zindex
    dic<-rep(NA,ncate)
    for (i in 1:ncate) dic[i]<-as.character(unique(unlist(data[zname])[which(data$zindex==i)]))
    if (!is.null(trtgrp)) trt<-which(dic==trtgrp)


  #use weight_gen() to obtain weights according to user's specification
  weight_gen<-function(AT){
    if(AT =='overlap'){
      tilt.h<-(1/apply(1/e.h,1,sum))
      allwt<-(1/e.h)*tilt.h
      wt<-rep(0,n)
      wt1<-rep(0,n)
      for(i in 1:ncate){
        wt[z==i]<-allwt[z==i,i]
        wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])}
    }

    else if (AT == 'IPW'){
      tilt.h<-rep(1,n)
      allwt<-1/e.h
      wt<-rep(0,n)
      wt1<-rep(0,n)
      for(i in 1:ncate){
        wt[z==i]<-allwt[z==i,i]
        wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])}
    }

    else if (AT == 'matching'){
      tilt.h<-apply(e.h, 1, min)
      allwt<-tilt.h/e.h
      wt<-rep(0,n)
      wt1<-rep(0,n)
      for(i in 1:ncate){
        wt[z==i]<-allwt[z==i,i]
        wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])}
    }

    else if (AT == 'entropy'){
      e.hclip<- pmax(e.h,1e-6)
      tilt.h<-(-apply(e.hclip*log(e.hclip) ,1,sum))
      allwt<-tilt.h/e.hclip
      wt<-rep(0,n)
      wt1<-rep(0,n)
      for(i in 1:ncate){
        wt[z==i]<-allwt[z==i,i]
        wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])}
    }

    else if (AT == 'treated'){
      if (is.null(trtgrp)){
        tilt.h<-e.h[,ncate]
        allwt<-tilt.h/e.h
        wt<-rep(0,n)
        wt1<-rep(0,n)
        for(i in 1:ncate){
          wt[z==i]<-allwt[z==i,i]
          wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])
        }
      }
      else{
        if (!trtgrp%in%dic) {
          warning("trtgrp not found in the treatment variable, argument ignored")
          trt<-ncate
        }
        tilt.h<-e.h[,trt]
        allwt<-tilt.h/e.h
        wt<-rep(0,n)
        wt1<-rep(0,n)
        for(i in 1:ncate){
          wt[z==i]<-allwt[z==i,i]
          wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])
        }
      }
    }

    else if (AT=="none"){
      tilt.h<-w<-rep(1,n)
      allwt<-data.frame(matrix(rep(1,ncate*n),ncol=ncate,nrow=n))
      wt<-rep(0,n)
      wt1<-rep(0,n)
      for(i in 1:ncate){
        wt[z==i]<-allwt[z==i,i]
        wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])}
    }else{
      stop("weight option not recognized","\n")
    }

    return(data.frame(wt,wt1,tilt.h))
  }

  #use wstat() function to calculate summary statistics--asd
  wstat<-function(covM,z,w,h){
    mres<-vres<-vres0<-NULL
    ncate<-length(table(z))
    hstd<-h/sum(h)

    pres<-as.numeric(colSums(covM*hstd))
    for(k in 1:ncate){
      nk<-sum(z==k)
      covMk<-covM[z==k,,drop=F]
      wk<-w[z==k]
      mk<-as.numeric(colSums(covMk*wk))
      mres<-cbind(mres,mk)
      vks<-as.numeric(colSums((covMk-rep(1,nk)%*%t(mk))^2*wk))
      vk<-sqrt(vks/(1-sum(wk^2)))
      vres<-cbind(vres,vk)
      vk0<-as.numeric(apply(covMk,2,sd))
      vres0<-cbind(vres0,vk0)
    }
    colnames(mres)<-paste0("Mean ",dic[1:ncate])
    colnames(vres)<-paste0("Weighted SD ",dic[1:ncate])
    colnames(vres0)<-paste0("Unweighted SD ", dic[1:ncate])
    rownames(mres)<-rownames(vres)<-rownames(vres0)<-colnames(covM)

    sx2<-1/ncate*rowSums(vres^2)
    sx20<-1/ncate*rowSums(vres0^2)

    #PSD
    psd<-(mres-pres)/sqrt(sx2)
    psd0<-(mres-pres)/sqrt(sx20)
    rownames(psd)<-rownames(psd0)<-colnames(covM)
    colnames(psd)<-paste0("PSD weighted var ",dic[1:ncate])
    colnames(psd0)<-paste0("PSD unweighted var ",dic[1:ncate])

    #ASD
    G<-combn(1:ncate,2)
    asd<-matrix(NA,ncol=dim(G)[2],nrow=dim(covM)[2])
    asd0<-matrix(NA,ncol=dim(G)[2],nrow=dim(covM)[2])
    rownames(asd)<-rownames(asd0)<-colnames(covM)
    colnames(asd)<-paste0("ASD weighted var ", apply(G,2,paste,collapse="-"))
    colnames(asd0)<-paste0("ASD unweighted var ", apply(G,2,paste,collapse="-"))


    for (g in 1:dim(G)[2]){
      asd[,g]<-(mres[,G[,g],drop=F][,1]-mres[,G[,g],drop=F][,2])/sqrt(sx2)
      asd0[,g]<-(mres[,G[,g],drop=F][,1]-mres[,G[,g],drop=F][,2])/sqrt(sx20)
    }

    #return(cbind(mres=mres,vres=vres,vres0=vres0,ASD.weighted.var=asd,
    #            ASD.unweighted.var=asd0,PSD.weighted.var=psd,PSD.unweighted.var=psd0))
    return(cbind(mres,vres,vres0,asd,asd0,psd,psd0))

  }

  #Matrix for effective sample size
  eff.sample.size<-matrix(rep(NA,ncate*(length(weight)+1)),nrow=ncate,ncol=(length(weight)+1))
  colnames(eff.sample.size)<-c("unweighted",weight)
  rownames(eff.sample.size)<-dic
  eff.sample.size[,1]<-as.numeric(table(data$zindex))

  uw<-weight_gen(AT="none")
  unweighted<-wstat(covM,z=data$zindex,w=uw[,2],h=uw[,3])

  ps.weights<-data.frame(Z=data[zname],zindex=data$zindex)
  #add all weight types and output weights and effective sample sizes
  for (i in 1:length(weight)){
    ATweights<-weight_gen(weight[i])
    ps.weights[weight[i]]<-ATweights[,2]

    sum_stat_AT<- wstat(covM,data$zindex,ATweights[,2],ATweights[,3])
    assign(weight[i], sum_stat_AT, envir = as.environment(-1))

    #add effective sample sizes
    for (j in 1:ncate){
      eff.sample.size[j,i+1]<-sum(ATweights[,1]*(data$zindex==j))^2/sum((ATweights[,1]*(data$zindex==j))^2)
    }
  }

  #output
  if (is.null(trtgrp)) trtgrp=dic[ncate]
  output<-list(trtgrp=trtgrp, propensity=e.h, ps.weights=ps.weights, ess=eff.sample.size, unweighted.sumstat=unweighted)

  for (i in 1:(length(weight))){
    output[[paste0(weight[i],".sumstat")]]<-get(weight[i])
  }

  if (delta>0) output[["trim"]]<-trim_sum

  class(output)<-'SumStat'
  output
}


