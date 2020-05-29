#SumStat_f() for not using user-supplied weights

##INPUT
# ps.formula: the propensity model
# data: the input dataframe(can have factor column)
# weight: vector of weighting methods of interest (can be "ATE","ATO","ATT")
# trtgrp: specify the reference group(only valid in ATT)
# delta: trim the data according to propensity

##OUTPUT
# object SumStat

#############################################################################################
SumStat_f<- function(ps.formula,ps.estimate=NULL,trtgrp=NULL,data,weight=c("ATO"),delta=0){
  #extract z name
  if(typeof(ps.formula)!="character"){
    ps.formula=Reduce(paste0,deparse(ps.formula))
  }
  ps.formula=gsub(" ","",ps.formula)
  zname=unlist(strsplit(ps.formula,'~'))[[1]][1]

  data[,zname]<-as.factor(data[,zname])

  categoryz<-sort(unique(unlist(data[zname])))
  if (length(categoryz)==1) stop("Treatment variable needs to have more than 1 category.","\n")

  #trim the data
  if(delta>0){
    tmp<-PStrim(data=data,ps.formula=ps.formula,zname=zname,delta=delta)
    data<-tmp[[1]]
    trim_sum<-tmp[[2]]
  }

  #number of categories and dataframe size
  z<-as.numeric(factor(unlist(data[zname])))
  data["zindex"]<-z

  ncate<-length(categoryz)
  n<-dim(data)[1]
  tabz<-table(unlist(data[zname]),data$zindex)
  #creat a dictionary for the original and recoded values in Z
  dic<-rep(NA,ncate)
  for (i in 1:ncate) dic[i]<-as.character(unique(unlist(data[zname])[which(data$zindex==i)]))
  if (!is.null(trtgrp)) trt<-which(dic==trtgrp)

  #estimate ps
  if(ncate==2){
    fit <- glm(formula = ps.formula, data=data,family = binomial(link = "logit"))
    e.h <- fit$fitted.values
    e.h <- cbind(1-e.h,e.h)
  }else{
    fit <- nnet::multinom(formula = ps.formula, data=data,maxit = 500, Hess = TRUE, trace = FALSE)
    e.h <- fit$fitted.values
  }
  e.h1<-e.h
  colnames(e.h1)<-categoryz

  #use weight_gen() to obtain weights according to user's specification
  weight_gen<-function(AT){
    if(AT =='ATO'){
      tilt.h<-(1/apply(1/e.h,1,sum))
      allwt<-(1/e.h)*tilt.h
      wt<-rep(0,n)
      wt1<-rep(0,n)
      for(i in 1:ncate){
        wt[z==i]<-allwt[z==i,i]
        wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])}
    }

    else if (AT == 'ATE'){
      tilt.h<-rep(1,n)
      allwt<-1/e.h
      wt<-rep(0,n)
      wt1<-rep(0,n)
      for(i in 1:ncate){
        wt[z==i]<-allwt[z==i,i]
        wt1[z==i]<-allwt[z==i,i]/sum(allwt[z==i,i])}
    }

    else if (AT == 'ATT'){
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
       if (!trtgrp%in%categoryz) {
          warning("trtgrp not found in the treatment variable, argument ignored")
          trt=ncate
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
    mres=vres=vres0=NULL
    ncate<-length(table(z))
    hstd=h/sum(h)

    pres=as.numeric(colSums(covM*hstd))
    for(k in 1:ncate){
      nk=sum(z==k)
      covMk=covM[z==k,,drop=F]
      wk=w[z==k]
      mk=as.numeric(colSums(covMk*wk))
      mres=cbind(mres,mk)
      vks=as.numeric(colSums((covMk-rep(1,nk)%*%t(mk))^2*wk))
      vk=sqrt(vks/(1-sum(wk^2)))
      vres=cbind(vres,vk)
      vk0=as.numeric(apply(covMk,2,sd))
      vres0=cbind(vres0,vk0)
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
    colnames(psd)<-colnames(psd0)<-dic[1:ncate]

    #ASD
    G<-combn(1:ncate,2)
    asd<-matrix(NA,ncol=dim(G)[2],nrow=dim(covM)[2])
    asd0<-matrix(NA,ncol=dim(G)[2],nrow=dim(covM)[2])
    rownames(asd)<-rownames(asd0)<-colnames(covM)
    colnames(asd)<-colnames(asd0)<-apply(G,2,paste,collapse="-")
    for (g in 1:dim(G)[2]){
      asd[,g]<-(mres[,G[,g],drop=F][,1]-mres[,G[,g],drop=F][,2])/sqrt(sx2)
      asd0[,g]<-(mres[,G[,g],drop=F][,1]-mres[,G[,g],drop=F][,2])/sqrt(sx20)
    }

    return(list(mres=mres,vres=vres,vres0=vres0,ASD.weighted.var=asd,
                ASD.unweighted.var=asd0,PSD.weighted.var=psd,PSD.unweighted.var=psd0))
  }

  #Matrix for effective sample size
  eff.sample.size<-matrix(rep(NA,ncate*(length(weight)+1)),nrow=ncate,ncol=(length(weight)+1))
  colnames(eff.sample.size)<-c("unweighted",weight)
  rownames(eff.sample.size)<-dic
  eff.sample.size[,1]<-as.numeric(table(data$zindex))

  covM<-as.data.frame(model.matrix(formula(ps.formula),data))
  covM<-as.matrix(covM)
  if (ncol(covM)>1)
  {   if(unique(covM[,1])==1) {
        covM<-covM[,-1,drop=F]
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
  output<-list(trtgrp=trtgrp, propensity=e.h1, ps.weights=ps.weights, ess=eff.sample.size, unweighted.sumstat=unweighted)

  for (i in 1:(length(weight))){
    output[[paste0(weight[i],".sumstat")]]<-get(weight[i])
  }

  if (delta>0) output[["trim"]]<-trim_sum

  class(output)<-'SumStat'
  output
}



