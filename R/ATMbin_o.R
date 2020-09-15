# Binary ATM with outcome value supplied


ATMbin_o <- function(ps.formula,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=200,out.formula=NULL,out.estimate=NULL,family=NULL,delta=0){

  #extract z name
  if(typeof(ps.formula)!="character"){
    ps.formula=Reduce(paste0,deparse(ps.formula))
  }
  ps.formula=gsub(" ","",ps.formula)
  zname=unlist(strsplit(ps.formula,'~'))[[1]][1]


  #set ordered group
  categoryz<-unique(unlist(data[zname]))
  z<-as.numeric(factor(unlist(data[zname])))-1
  oldlevel<-categoryz[order(unique(z))]

  data[zname]<-z
  data1<-data[,colnames(data)!=yname]
  y=unlist(data[yname])


  # summary statistics
  n1 <- sum(z)
  n0 <- sum(1-z)
  n <- n0 + n1


  #estimated ps
  fit <- glm(ps.formula, family = binomial(link = "logit"),data=data1)
  e.h <- as.numeric(fit$fitted.values)
  beta.h<-as.numeric(coef(fit))
  W<- model.matrix(fit)


  # point estimate
  htilt.h<-pmin(e.h,1-e.h)
  mu1.h <- sum(z*y*htilt.h/e.h) / sum(z*htilt.h/e.h)
  mu0.h <- sum((1-z)*y*htilt.h/(1-e.h)) / sum((1-z)*htilt.h/(1-e.h))


  #augmentation and no bootstrap
  if(sum(colnames(out.estimate)%in%categoryz)<2){
    out.estimate<-as.matrix(out.estimate)
    m0.h<-out.estimate[,1]
    m1.h<-out.estimate[,2]
    warning("wrong column name set for out.estimate, treatment set as: ",oldlevel)
  }else{
    out.estimate<-as.matrix(out.estimate)
    out.estimate<-out.estimate[,match(oldlevel,colnames(out.estimate))]
    m0.h<-out.estimate[,1]
    m1.h<-out.estimate[,2]
  }


  #calculate the aumentation term
  aug0h.h<-sum(m0.h*htilt.h)/sum(htilt.h)
  aug0z.h<-sum((1-z)*htilt.h*m0.h/(1-e.h))/sum((1-z)*htilt.h/(1-e.h))
  aug1h.h<-sum(m1.h*htilt.h)/sum(htilt.h)
  aug1z.h<-sum(z*htilt.h*m1.h/e.h)/sum(z*htilt.h/e.h)

  mu0aug.h<-mu0.h+aug0h.h-aug0z.h
  mu1aug.h<-mu1.h+aug1h.h-aug1z.h

  #Estimate the sandwich variance after augmentation
  p<-ncol(W)
  #define the m estimator
  phi<-function(theta){
    mu1<-theta[1]
    aug1h<-theta[2]
    aug1z<-theta[3]
    mu0<-theta[4]
    aug0h<-theta[5]
    aug0z<-theta[6]

    #component in m estimator
    f1<-z*htilt.h/e.h*(y-mu1)
    f2<-(m1.h-aug1h)*htilt.h
    f3<-z*htilt.h/e.h*(m1.h-aug1z)
    f4<-(1-z)*htilt.h/(1-e.h)*(y-mu0)
    f5<-(m0.h-aug0h)*htilt.h
    f6<-(1-z)*htilt.h/(1-e.h)*(m0.h-aug0z)
    f<-rbind(f1,f2,f3,f4,f5,f6)
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

  theta.h<-c(mu1.h,aug1h.h,aug1z.h,mu0.h,aug0h.h,aug0z.h)
  Atheta<-jacobian(mphi,theta.h)
  invAtheta <- solve(Atheta)
  a0<-c(0,0,0,1,1,-1)
  a1<-c(1,1,-1,0,0,0)

  V1<-invAtheta%*%Omega(theta.h)%*%t(invAtheta)/n

  a<-rbind(a0,a1)
  covmu<-a%*%V1%*%t(a)
  muhat<-c(mu0aug.h, mu1aug.h)
  muboot<-NULL
  names(muhat)<-oldlevel
  colnames(covmu)<-rownames(covmu)<-oldlevel


  e.h<-cbind(1-e.h,e.h)
  colnames(e.h)<-oldlevel
  out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
  class(out)<-'PSweight'
  out
}





