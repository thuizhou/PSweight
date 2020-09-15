#SumStat_p() for using user-supplied weights



#############################################################################################
SumStat_p<- function(ps.estimate=NULL,trtgrp=NULL,zname=NULL,xname=NULL,Z=NULL,covM=NULL,data=NULL,weight=c("overlap"),delta=0){

  #print message if more than nesscary info supplied
  if ((!is.null(Z))&&(!is.null(zname))==T) warning("both Z and zname supplied, zname ignored")
  if ((!is.null(xname))&&(!is.null(covM))==T) warning("both covM and xname supplied, xname ignored")

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
  else if (AT == "treated"){
    if (is.null(trtgrp))
    {
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
      if (!trtgrp%in%categoryz)
      { warning("trtgrp not found in the treatment variale, argument ignored ")
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
    colnames(asd)<-colnames(asd0)<-apply(G,2,paste,collapse="")
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

  uw<-weight_gen(AT="none")
  unweighted<-wstat(covM,z=data$zindex,w=uw[,2],h=uw[,3])

  ps.weights<-data.frame(Z=data[zname],zindex=data$zindex)
  #add all weight types and output weights and effective sample sizes
  for (i in 1:length(weight)){
    ATweights<-weight_gen(weight[i])
    ps.weights[weight[i]]<-ATweights[,2]

    sum_stat_AT<- wstat(covM,data$zindex,ATweights[,2],ATweights[,3])
    assign(weight[i], sum_stat_AT, envir=as.environment(-1))

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
