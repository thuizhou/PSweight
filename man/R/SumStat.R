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
#' @param method a character to specify the method for estimating propensity scores. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param ps.control a list to specify additional options when \code{method} is set to \code{"gbm"} or \code{"SuperLearner"}.
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
#' \code{"glm"} is the default method for propensity score estimation. Logistic regression will be used for binary outcomes,
#' and multinomial logistic regression will be used for outcomes with more than two categories. The alternative method option of \code{"gbm"} serves as an API to call the \code{gbm()} function from the
#' \code{gbm} package. Additional argument in the \code{gbm()} function can be supplied through the \code{ps.control=list()} argument in \code{SumStat()}. Please refer to the user manual of the gbm package for all the
#' allowed arguments. Currently, models for binary or multinomial treatment will be automatically chosen based on the number of treatment categories.
#' \code{"SuperLearner"} is also allowed in the \code{method} argument to pass the propensity score estimation to the \code{SuperLearner()} function in SuperLearner package.
#' Currently, the SuperLearner method only supports binary treatment with the default method set to \code{"SL.glm"}. The estimation approach is default to \code{"method.NNLS"} in the \code{SumStat()} function.
#' Prediction algorithm and other tuning parameters can also be passed through \code{ps.control=list()} to \code{SumStat()}. Please refer to the user manual of the \code{SuperLearner} package for all the allowed specifications.
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
#' Greenwell B., Boehmke B.,Cunningham J, GBM Developers (2020) gbm: Generalized Boosted Regression Models. Cran: https://cran.r-project.org/web/packages/gbm/index.html
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
#' Polley E., LeDell E., Kennedy C., Lendle S., van der Laan M. (2019) SuperLearner: Super Learner Prediction. Cran: https://cran.r-project.org/web/packages/SuperLearner/index.html
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
#' Propensity score weighting under limited overlap and model misspecification. Statistical Methods in Medical Research (Online)
#'
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
#' #summary(msstat)
#'
#' # importing user-supplied propensity scores "e.h"
#' fit <- nnet::multinom(formula=ps.formula, data=psdata, maxit=500, trace=FALSE)
#' e.h <- fit$fitted.values
#' varname <- c("cov1","cov2","cov3","cov4","cov5","cov6")
#' msstat0 <- SumStat(zname="trt", xname=varname, data=psdata, ps.estimate=e.h,
#'    trtgrp="2",  weight=c("IPW","overlap","treated","entropy","matching"))
#' summary(msstat0)
#'
#' @import nnet SuperLearner gbm
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
#'
#'
SumStat<- function(ps.formula=NULL,ps.estimate=NULL,trtgrp=NULL,Z=NULL,covM=NULL,zname=NULL,xname=NULL,data=NULL,weight="overlap",delta=0,method='glm',ps.control=list()){
  #if not user-supplied weights
  if (is.null(ps.estimate)){
    #extract z name
    zname<-all.vars(ps.formula)[1]

    #set ordered group
    facz<-as.factor(unlist(data[,zname]))
    data[,zname]<-facz

    #treatment label
    #creat a dictionary for the original and recoded values in Z
    dic<-levels(facz)
    ncate<-length(dic) #number of categories
    if (!is.null(trtgrp)) trt<-which(dic==trtgrp)

    if (ncate==1) stop("Treatment variable needs to have more than 1 category.","\n")

    #number of categories and data frame size

    #trim the data
    if(delta>0){
      tmp<-do.call(PStrim,list(data=data,ps.formula = ps.formula, zname=zname,delta=delta,optimal=FALSE,method=method,ps.control=ps.control))
      data<-tmp[[1]]
      trim_sum<-tmp[[2]]
    }


    #fit propensity score model
    e.h<-do.call(PSmethod,c(list(ps.formula = ps.formula, method=method, data=data,ncate=ncate),ps.control))$e.h

    #post-trimming processing
    z<-as.numeric(data[,zname])
    n<-length(z) #total obs
    data["zindex"]<-z

  }else{ #user supplied weights
    #print message if more than nesscary info supplied
    if ((!is.null(Z))&&(!is.null(zname))==T) warning("both Z and zname supplied, zname ignored")
    if ((!is.null(xname))&&(!is.null(covM))==T) warning("both covM and xname supplied, xname ignored")
    if (!is.null(ps.formula)) warning("ps.estimate supplied and ps.formula ignored. Please supply covariate information through covM or both xname and data")
    if (length(ps.control>0)) warning ("ps.control argument ignored for user-supplied ps.estiamte")

    #extract input info
    if (!is.null(Z)){
      trts<-sort(trimws(unique(Z)))
      if (is.null(covM)){
        if (!is.null(xname)){
          if (is.null(data)) {
            stop ("need to specify data, or covM","\n")
          }
          else if (!all(xname%in%colnames(data)==T)){
            stop("one or more item in xname not found in covM","\n")
          }else{
            covM=data[,xname,drop=F]
          }
        }else {
          stop("need to specify covM, or xname and data","\n")
        }
      }
    }else if(!is.null(zname)){

      if (is.null(data)){
        stop ("need to specify zname and data, or trt Z","\n")
      }else if (!zname%in%colnames(data)){
        stop("zname not found in data","\n")
      }
      else{
        Z<-data[,zname]
      }
      if (is.null(covM)){
        if (!is.null(xname)){
          if (!all(xname%in%colnames(data)==T)){
            stop("one or more item in xname not found in covM","\n")
          }else{
            covM=data[,xname,drop=F]
          }
        }else {
          stop("need to specify covM, or xname and data","\n")
        }
      }
    }else { #no Z or zname
      stop("need to specify trt Z, or zname and data","\n")
    }

    #transform factor/categorical column, if any
    cat_index<-c()
    for(i in 1:ncol(covM)){
      if(class(covM[,i])=="factor") {
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

    #data validity checks
    trts<-sort(trimws(unique(Z)))

    if (length(trts)==1) {
      stop("treatment variable needs to have more than 1 category.","\n")
    }

    #one column matrix or one row matrix
    if(is.matrix(ps.estimate)){
      if(sum(dim(ps.estimate)==1)>0){
        ps.estimate<-c(ps.estimate)
      }
    }

    #check if vector
    if (is.vector(ps.estimate)) {
      ps.estimate=data.frame(1-ps.estimate,ps.estimate)
      colnames(ps.estimate)=NULL
      if(!is.null(trtgrp)){
        newname<-unique(c(trtgrp,trts))[2:1]
        colnames(ps.estimate)<-newname
      }
    }
    if (dim(ps.estimate)[1]!=dim(covM)[1]) {
      stop("number of obs do not match across Z, covM, data, and ps.estimate","\n")
    }
    if (dim(ps.estimate)[1]!=length(Z)) {
      stop("number of obs do not match across Z, covM, data, and ps.estimate","\n")
    }
    if (length(Z)!=dim(covM)[1]) {
      stop("number of obs do not match across Z, covM, data, and ps.estimate","\n")
    }

    if (dim(ps.estimate)[2]!=length(trts)) {
      stop("number of trt in Z/zname does not match number of trts/columns in ps.estimate","\n")
    }

    if (is.null(colnames(ps.estimate))) {
      colnames(ps.estimate)=trts
    }else if (sum(trimws(colnames(ps.estimate))%in%trts)<length(trts)){
      warning("colnames of ps.estimate do not match values in Z; colnames of ps.estimate reassigned alphebatically")
      colnames(ps.estimate)=trts
    }else {
      ps.estimate<-ps.estimate[,trts]
    }

    if(is.null(colnames(covM))){
      colnames(covM)<-paste('VAR',1:(dim(covM)[2]))
    }

    #construct data for potential trimming
    data<-data.frame(Z,covM)
    covM<-as.matrix(covM)
    zname<-"Z"
    xname<-colnames(covM)
    categoryz<-trts

    #trim the data
    if(delta>0){
      tmp<-PStrim(data=data,ps.formula=NULL,zname=zname,ps.estimate=ps.estimate,delta=delta)
      data<-tmp[[1]]
      covM<-as.matrix(data[,-1])
      trim_sum<-tmp[[2]]
      ps.estimate<-tmp[[3]]
    }

    #complete the rest info
    z<-as.numeric(factor(unlist(data[zname])))
    data["zindex"]<-z
    e.h<-ps.estimate

    #dimension of data
    n<-dim(data)[1]
    ncate<-length(categoryz)

    #correspondance btw Z and zindex
    dic<-rep(NA,ncate)
    for (i in 1:ncate) dic[i]<-as.character(unique(unlist(data[zname])[which(data$zindex==i)]))
    if (!is.null(trtgrp)) trt<-which(dic==trtgrp)

  }



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

  if(is.null(ps.estimate)){
    covM<-as.data.frame(model.matrix(formula(ps.formula),data))
    covM<-as.matrix(covM)
    if (ncol(covM)>1)
    {   if(unique(covM[,1])==1) {
      covM<-covM[,-1,drop=F]
    }
    }
  }

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


