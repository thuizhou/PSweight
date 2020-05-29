# Binary ATT with outcome value supplied


ATTbin_o <- function(ps.formula,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=200,out.formula=NULL,out.estimate=NULL,family=NULL,delta=0){

  #extract z name
  if(typeof(ps.formula)!="character"){
    ps.formula=Reduce(paste0,deparse(ps.formula))
  }
  ps.formula=gsub(" ","",ps.formula)
  zname=unlist(strsplit(ps.formula,'~'))[[1]][1]



  #set the reference group
  categoryz<-unique(unlist(data[zname]))
  oldlevel0<-categoryz[order(unique(as.numeric(factor(unlist(data[zname])))))]
  if(is.null(trtgrp)){
    z<-as.numeric(factor(unlist(data[zname])))-1
    oldlevel<-categoryz[order(unique(z))]
  }else{
    categoryz<-categoryz[-which(categoryz==trtgrp)]
    categoryz<-c(categoryz,trtgrp)
    z<-as.numeric(factor(unlist(data[zname]),levels = categoryz))-1
    oldlevel<-categoryz
  }
  data[zname]<-z
  data1<-data[,colnames(data)!=yname]

  matchlevel<-match(oldlevel0,oldlevel)

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
  mu1.h <- sum(z*y) / sum(z)
  mu0.h <- sum((1-z)*y*e.h/(1-e.h)) / sum(e.h*(1-z)/(1-e.h))

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

  #calculate the aumentation term
  aug0h.h<-sum(m0.h*e.h)/sum(e.h)
  aug0z.h<-sum((1-z)*m0.h*e.h/(1-e.h))/sum((1-z)*e.h/(1-e.h))
  aug1h.h<-sum(m1.h*e.h)/sum(e.h)
  aug1z.h<-sum(z*m1.h)/sum(z)

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
    f1<-z*(y-mu1)
    f2<-(m1.h-aug1h)*e
    f3<-z*(m1.h-aug1z)
    f4<-(1-z)*e/(1-e)*(y-mu0)
    f5<-(m0.h-aug0h)*e
    f6<-(1-z)*e/(1-e)*(m0.h-aug0z)
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

  tryCatch({
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
      f1<-z*(y-mu1)
      f2<-(m1.h-aug1h)*e.h
      f3<-z*(m1.h-aug1z)
      f4<-(1-z)*e.h/(1-e.h)*(y-mu0)
      f5<-(m0.h-aug0h)*e.h
      f6<-(1-z)*e.h/(1-e.h)*(m0.h-aug0z)

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

  muhat<-muhat[matchlevel]
  covmu<-covmu[matchlevel,matchlevel]

  e.h<-cbind(1-e.h,e.h)
  colnames(e.h)<-oldlevel
  e.h<-e.h[,matchlevel]
  out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel0), trtgrp=oldlevel[2])
  class(out)<-'PSweight'
  out
}
