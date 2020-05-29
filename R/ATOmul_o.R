# Multiple ATO with outcome value supplied



ATOmul_o<-function(ps.formula,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=200,out.formula=NULL,out.estimate=NULL,family=NULL,delta=0){

  #extract the name for treatment
  if(typeof(ps.formula)!="character"){
    ps.formula=Reduce(paste0,deparse(ps.formula))
  }
  ps.formula=gsub(" ","",ps.formula)
  zname=unlist(strsplit(ps.formula,'~'))[[1]][1]



  y=unlist(data[yname])



  categoryz<-unique(unlist(data[zname]))
  z<-as.numeric(factor(unlist(data[zname])))
  oldlevel<-categoryz[order(unique(z))]


  #dataframe used for propensity model excluding yname
  data1<-data[,!names(data)%in%yname]
  data1[zname]<-z

  # summary statistics
  tabz<-table(z)
  ncate<-length(tabz)
  n<-length(y)

  # estimate ps
  fit <- multinom(formula = ps.formula, data=data1,maxit = 500, Hess = TRUE, trace = FALSE)

  W<- model.matrix(formula(ps.formula),data1)                      # design matrix (including intercept)
  e.h <- fit$fitted.values
  beta.h<-as.numeric(t(coef(fit)))

  #use the harmonic mean to get the weight
  htilt.h<-(1/apply(1/e.h,1,sum))
  allwt<-(1/e.h)*htilt.h
  wt<-rep(0,n)
  for(i in 1:ncate){
    wt[z==i]<-allwt[z==i,i]
  }

  #the point estimate
  mu.h<-c()
  for(i in 1:ncate){
    ytmp<-y[z==i]
    wttmp<-wt[z==i]
    mu.h<-c(mu.h,sum(ytmp*wttmp)/sum(wttmp))
  }

  #augmentation and no bootstrap
  if(sum(colnames(out.estimate)%in%categoryz)<2){
    out.estimate<-as.matrix(out.estimate)
    m.h<-as.matrix(out.estimate)
    warning("wrong column name set for out.estimate, treatment set as: ",oldlevel)
  }else{
    out.estimate<-as.matrix(out.estimate)
    out.estimate<-out.estimate[,match(oldlevel,colnames(out.estimate))]
    m.h<-out.estimate
  }


  #calculate the augmentation term and updata mu
  augz.h<-c()
  augh.h<-c()
  for(i in 1:ncate){
    mtmp<-m.h[,i]
    wttmp<-wt[z==i]

    augztmp<-sum(mtmp[z==i]*wttmp)/sum(wttmp)
    aughtmp<-sum(htilt.h*mtmp)/sum(htilt.h)
    augz.h<-c(augz.h,augztmp)
    augh.h<-c(augh.h,aughtmp)
  }

  muhat<-mu.h+augh.h-augz.h
  names(muhat)<-oldlevel

  #Estimate the sandwich variance after augmentation
  p<-ncol(W)
  #define the m estimator
  phi<-function(theta){
    mu<-theta[1:ncate]
    augz<-theta[(ncate+1):(2*ncate)]
    augh<-theta[(2*ncate+1):(3*ncate)]
    beta<-theta[(3*ncate+1):(3*ncate+(ncate-1)*p)]

    e<-rep(1,n)
    for(i in 1:(ncate-1)){
      etmp=exp(W%*%beta[((i-1)*p+1):(i*p)])
      e<-cbind(e,etmp)
    }
    e<-e/(apply(e,1,sum))
    htilt<-(1/apply(1/e,1,sum))
    wtt<-(1/e)*htilt



    fmu<-c()
    fbeta<-c()
    faugz<-c()
    faugh<-c()

    for(i in 1:ncate){
      fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
      fmu<-rbind(fmu,fmutmp)
      faugztmp<-(m.h[,i]-augz[i])*wtt[,i]*(z==i)
      faughtmp<-htilt*(m.h[,i]-augh[i])
      faugz<-rbind(faugz,faugztmp)
      faugh<-rbind(faugh,faughtmp)
    }
    for(i in 2:ncate){
      fbetatmp<-((z==i)-e[,i])*W
      fbeta<-cbind(fbeta,fbetatmp)
    }


    f<-rbind(fmu,faugz,faugh,t(fbeta))
    return(f)
  }
  #define score function
  mphi<-function(theta){
    rowMeans(phi(theta))
  }

  #define the meat B, covariance operator
  Omega<-function(theta){
    phis<-phi(theta)
    return(tcrossprod(phis)/n)
  }

  #choose conservative or not
  conser<-1

  tryCatch( {
    theta.h<-c(mu.h,augz.h,augh.h,beta.h)
    Atheta<-jacobian(mphi,theta.h)
    invAtheta <- solve(Atheta)
    conser<-0
  },error = function(w) {
    warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
    })
  #if not pd than use conservative
  if(conser==1){
    phi<-function(theta){
      mu<-theta[1:ncate]
      augz<-theta[(ncate+1):(2*ncate)]
      augh<-theta[(2*ncate+1):(3*ncate)]
      htilt<-(1/apply(1/e.h,1,sum))
      wtt<-(1/e.h)*htilt
      fmu<-c()
      faugz<-c()
      faugh<-c()

      for(i in 1:ncate){
        fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
        fmu<-rbind(fmu,fmutmp)
        faugztmp<-(m.h[,i]-augz[i])*wtt[,i]*(z==i)
        faughtmp<-htilt*(m.h[,i]-augh[i])
        faugz<-rbind(faugz,faugztmp)
        faugh<-rbind(faugh,faughtmp)
      }

      f<-rbind(fmu,faugz,faugh)
      return(f)
    }
    theta.h<-c(mu.h,augz.h,augh.h)
    Atheta<-jacobian(mphi,theta.h)
    invAtheta <- solve(Atheta)
  }

  V<-invAtheta%*%Omega(theta.h)%*%t(invAtheta)/n
  a<-c()
  for(i in 1:ncate){
    atmp<-rep(0,length(theta.h))
    atmp[c(i,2*ncate+i)]<-1
    atmp[ncate+i]<--1
    a<-rbind(a,atmp)
  }
  covmu<-a%*%V%*%t(a)
  muboot<-NULL
  colnames(covmu)<-rownames(covmu)<-oldlevel


  colnames(e.h)<-oldlevel
  out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
  class(out)<-'PSweight'
  out
}






