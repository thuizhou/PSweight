# Binary ATO with outcome value supplied


ATObin_o<-function(ps.formula,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=200,out.formula=NULL,out.estimate=NULL,family=NULL,delta=0){

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
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1  # total sample size


  # estimate ps
  fit <- glm(ps.formula, family = binomial(link = "logit"),data=data1)
  W<- model.matrix(fit)                      # design matrix (including intercept)
  e.h <- as.numeric(fit$fitted.values)
  beta.h<-as.numeric(coef(fit))


  # point estimate
  mu1.h <- sum(z*y*(1-e.h)) / sum(z*(1-e.h))
  mu0.h <- sum((1-z)*y*e.h) / sum((1-z)*e.h)

  #augmentation and no bootstrap
  if(sum(colnames(out.estimate)%in%categoryz)<2){
    out.estimate<-as.matrix(out.estimate)
    m0.h<-out.estimate[,1]
    m1.h<-out.estimate[,2]
    warning("wrong column name set for out.estimate, treatment set as: ",oldlevel)
  }else{
    out.estimate<-out.estimate[,match(oldlevel,colnames(out.estimate))]
    m0.h<-out.estimate[,1]
    m1.h<-out.estimate[,2]
  }


  #calculate the augmentation term and updata tau
  aug0h.h<-sum(e.h*(1-e.h)*m0.h)/sum(e.h*(1-e.h))
  aug0z.h<-sum((1-z)*m0.h*e.h)/sum((1-z)*e.h)
  aug1h.h<-sum(e.h*(1-e.h)*m1.h)/sum(e.h*(1-e.h))
  aug1z.h<-sum(z*m1.h*(1-e.h))/sum(z*(1-e.h))

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
    beta<-theta[7:(p+6)]
    e<-plogis(c(W%*%beta))


    #component in m estimator
    f1<-(1-e)*z*(y-mu1)
    f2<-(1-e)*e*(m1.h-aug1h)
    f3<-(1-e)*z*(m1.h-aug1z)
    f4<-e*(1-z)*(y-mu0)
    f5<-(1-e)*e*(m0.h-aug0h)
    f6<-e*(1-z)*(m0.h-aug0z)
    f7<-W*(z-e)

    f<-rbind(f1,f2,f3,f4,f5,f6,t(f7))
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

  #when not pd
  tryCatch( {
    theta.h<-c(mu1.h,aug1h.h,aug1z.h,mu0.h,aug0h.h,aug0z.h,beta.h)
    Atheta<-jacobian(mphi,theta.h)
    invAtheta <- solve(Atheta)
    a0<-c(0,0,0,1,1,-1,rep(0,p))
    a1<-c(1,1,-1,0,0,0,rep(0,p))
    conser<-0
  },error = function(w) {
    warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
  })

  #if not pd than use conservative
  if(conser==1){
    #define the m estimator
    phi<-function(theta){
      mu1<-theta[1]
      aug1h<-theta[2]
      aug1z<-theta[3]
      mu0<-theta[4]
      aug0h<-theta[5]
      aug0z<-theta[6]

      #component in m estimator
      f1<-(1-e.h)*z*(y-mu1)
      f2<-(1-e.h)*e.h*(m1.h-aug1h)
      f3<-(1-e.h)*z*(m1.h-aug1z)
      f4<-e.h*(1-z)*(y-mu0)
      f5<-(1-e.h)*e.h*(m0.h-aug0h)
      f6<-e.h*(1-z)*(m0.h-aug0z)
      f<-rbind(f1,f2,f3,f4,f5,f6)
      return(f)
    }

    theta.h<-c(mu1.h,aug1h.h,aug1z.h,mu0.h,aug0h.h,aug0z.h)
    Atheta<-jacobian(mphi,theta.h)
    invAtheta <- solve(Atheta)
    a0<-c(0,0,0,1,1,-1)
    a1<-c(1,1,-1,0,0,0)
  }


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




