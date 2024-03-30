## Notation #########################
# single e stands for no augmentation
# ea stands for augmentation
# c represents conservative variance
#####################################

################## Sandwich Variance in binary group case####################################################
sand_bin<-function(z,y,n,ftilt,thetaest,W=NULL,XY=NULL,eest,m0est,m1est,family="gaussian",offset.e,type="e"){
  p<-ncol(W)
  q<-ncol(XY)

  #matrix to calucalte mu0,mu1 correlation
  a<-matrix(0,2,length(thetaest))

  ##score function
  if(type=="e"){
    phi<-function(theta){
      mu0<-theta[1]
      mu1<-theta[2]
      beta<-theta[3:(p+2)]
      e<-plogis(c(W%*%beta))
      f1<-(y-mu0)*(1-z)*ftilt(e)/(1-e)
      f2<-(y-mu1)*z*ftilt(e)/e
      f3<-W*(z-e)

      f<-rbind(f1,f2,t(f3))
      return(f)
    }
    a[1,1]<-1
    a[2,2]<-1
  }else if(type=="ec"){
    phi<-function(theta){
      mu0<-theta[1]
      mu1<-theta[2]

      f1<-(y-mu0)*(1-z)*ftilt(eest)/(1-eest)
      f2<-(y-mu1)*z*ftilt(eest)/eest

      f<-rbind(f1,f2)
      return(f)
    }
    a[1,1]<-1
    a[2,2]<-1
  }else if(type=="ea"){
    phi<-function(theta){
      mu0<-theta[1]
      aug0h<-theta[2]
      aug0z<-theta[3]
      mu1<-theta[4]
      aug1h<-theta[5]
      aug1z<-theta[6]
      beta<-theta[7:(p+6)]
      gamma0<-theta[(p+7):(p+6+q)]
      gamma1<-theta[(p+7+q):(p+2*q+6)]
      e<-plogis(c(W%*%beta))
      if (family=='gaussian'){
        m0<-c(XY%*%gamma0)
        m1<-c(XY%*%gamma1)
      }else if(family=='binomial'){
        m0<-plogis(c(XY%*%gamma0))
        m1<-plogis(c(XY%*%gamma1))
      }else if(family=='poisson'){
        m0<-exp(c(XY%*%gamma0))*offset.e
        m1<-exp(c(XY%*%gamma1))*offset.e
      }

      f1<-(1-z)*(y-mu0)*ftilt(e)/(1-e)
      f2<-ftilt(e)*(m0-aug0h)
      f3<-(1-z)*(m0-aug0z)*ftilt(e)/(1-e)
      f4<-z*(y-mu1)*ftilt(e)/e
      f5<-ftilt(e)*(m1-aug1h)
      f6<-z*(m1-aug1z)*ftilt(e)/e
      f7<-XY*(y-m0)*(1-z)
      f8<-XY*(y-m1)*z
      f9<-W*(z-e)

      f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8),t(f9))
      return(f)
    }
    a[1,1:3]<-c(1,1,-1)
    a[2,4:6]<-c(1,1,-1)
  }else if(type=='ecac'){
    phi<-function(theta){
      mu0<-theta[1]
      aug0h<-theta[2]
      aug0z<-theta[3]
      mu1<-theta[4]
      aug1h<-theta[5]
      aug1z<-theta[6]

      f1<-(1-z)*(y-mu0)*ftilt(eest)/(1-eest)
      f2<-ftilt(eest)*(m0est-aug0h)
      f3<-(1-z)*(m0est-aug0z)*ftilt(eest)/(1-eest)
      f4<-z*(y-mu1)*ftilt(eest)/eest
      f5<-ftilt(eest)*(m1est-aug1h)
      f6<-z*(m1est-aug1z)*ftilt(eest)/eest

      f<-rbind(f1,f2,f3,f4,f5,f6)
      return(f)
    }
    a[1,1:3]<-c(1,1,-1)
    a[2,4:6]<-c(1,1,-1)
  }else if(type=='eca'){
    #define the m estimator
    phi<-function(theta){
      mu0<-theta[1]
      aug0h<-theta[2]
      aug0z<-theta[3]
      mu1<-theta[4]
      aug1h<-theta[5]
      aug1z<-theta[6]
      gamma0<-theta[7:(6+q)]
      gamma1<-theta[(7+q):(2*q+6)]
      if (family=='gaussian'){
        m1<-c(XY%*%gamma1)
        m0<-c(XY%*%gamma0)
      }else if(family=='binomial'){
        m1<-plogis(c(XY%*%gamma1))
        m0<-plogis(c(XY%*%gamma0))
      }else if(family=='poisson'){
        m1<-exp(c(XY%*%gamma1))*offset.e
        m0<-exp(c(XY%*%gamma0))*offset.e
      }

      f1<-(1-z)*(y-mu0)*ftilt(eest)/(1-eest)
      f2<-ftilt(eest)*(m0-aug0h)
      f3<-(1-z)*(m0-aug0z)*ftilt(eest)/(1-eest)
      f4<-z*(y-mu1)*ftilt(eest)/eest
      f5<-ftilt(eest)*(m1-aug1h)
      f6<-z*(m1-aug1z)*ftilt(eest)/eest
      f7<-XY*(y-m0)*(1-z)
      f8<-XY*(y-m1)*z

      f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8))
      return(f)
    }
    a[1,1:3]<-c(1,1,-1)
    a[2,4:6]<-c(1,1,-1)
  }else if(type=="eac"){
      phi<-function(theta){
        mu0<-theta[1]
        aug0h<-theta[2]
        aug0z<-theta[3]
        mu1<-theta[4]
        aug1h<-theta[5]
        aug1z<-theta[6]
        beta<-theta[7:(p+6)]
        e<-plogis(c(W%*%beta))

        f1<-(1-z)*(y-mu0)*ftilt(e)/(1-e)
        f2<-ftilt(e)*(m0est-aug0h)
        f3<-(1-z)*(m0est-aug0z)*ftilt(e)/(1-e)
        f4<-z*(y-mu1)*ftilt(e)/e
        f5<-ftilt(e)*(m1est-aug1h)
        f6<-z*(m1est-aug1z)*ftilt(e)/e
        f7<-W*(z-e)

        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7))
        return(f)
      }
      a[1,1:3]<-c(1,1,-1)
      a[2,4:6]<-c(1,1,-1)
  }


  mphi<-function(theta){
    rowMeans(phi(theta))
  }

  #define the meat B, covariance operator
  Omega<-function(theta){
    phis<-phi(theta)
    return(tcrossprod(phis)/n)
  }

  Atheta<-jacobian(mphi,thetaest)
  invAtheta <- MASS::ginv(Atheta)

  #calculate the sandwich variance
  Vtmp<-invAtheta%*%Omega(thetaest)%*%t(invAtheta)/n

  covmu<-a%*%Vtmp%*%t(a)

  return(covmu)
}


################## Sandwich Variance in multiple group case####################################################

sand_mul<-function(z,y,n,ncate,ftilt,thetaest,W=NULL,XY=NULL,eest,mest,family="gaussian",offset.e,type="e"){
  p<-ncol(W)
  q<-ncol(XY)

  #matrix to calucalte correlation
  a<-matrix(0,ncate,length(thetaest))

  ##score function
  if(type=="e"){
    p<-ncol(W)
    phi<-function(theta){
      mu<-theta[1:ncate]
      beta<-theta[(ncate+1):(ncate+(ncate-1)*p)]
      e<-rep(1,n)
      for(i in 1:(ncate-1)){
        etmp<-exp(W%*%beta[((i-1)*p+1):(i*p)])
        e<-cbind(e,etmp)
      }
      e<-e/(apply(e,1,sum))
      wtt<-(1/e)*ftilt(e)
      fmu<-c()
      fbeta<-c()
      for(i in 1:ncate){
        fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
        fmu<-rbind(fmu,fmutmp)
      }
      for(i in 2:ncate){
        fbetatmp<-((z==i)-e[,i])*W
        fbeta<-cbind(fbeta,fbetatmp)
      }

      f<-rbind(fmu,t(fbeta))
      return(f)
    }
    diag(a)<-1
  }else if(type=="ec"){
    phi<-function(theta){
      mu<-theta[1:ncate]
      wtt<-(1/eest)*ftilt(eest)
      fmu<-c()
      fbeta<-c()
      for(i in 1:ncate){
        fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
        fmu<-rbind(fmu,fmutmp)
      }

      f<-rbind(fmu)
      return(f)
    }
    diag(a)<-1
  }else if(type=="ea"){
    phi<-function(theta){
      mu<-theta[1:ncate]
      augz<-theta[(ncate+1):(2*ncate)]
      augh<-theta[(2*ncate+1):(3*ncate)]
      beta<-theta[(3*ncate+1):(3*ncate+(ncate-1)*p)]
      gamma<-theta[(3*ncate+(ncate-1)*p+1):(3*ncate+(ncate-1)*p+ncate*q)]

      e<-rep(1,n)
      for(i in 1:(ncate-1)){
        etmp<-exp(W%*%beta[((i-1)*p+1):(i*p)])
        e<-cbind(e,etmp)
      }
      e<-e/(apply(e,1,sum))
      wtt<-(1/e)*ftilt(e)

      #outcome aug
      m<-c()
      if (family=='gaussian'){
        for(i in 1:ncate){
          mtmp<-c(XY%*%gamma[((i-1)*q+1):(i*q)])
          m<-cbind(m,mtmp)}
      }else if(family=='binomial'){
        for(i in 1:ncate){
          mtmp<-plogis(c(XY%*%gamma[((i-1)*q+1):(i*q)]))
          m<-cbind(m,mtmp)}
      }else if(family=='poisson'){
        for(i in 1:ncate){
          mtmp<-exp(c(XY%*%gamma[((i-1)*q+1):(i*q)]))*offset.e
          m<-cbind(m,mtmp)}
      }

      fmu<-c()
      fbeta<-c()
      fgamma<-c()
      faugz<-c()
      faugh<-c()

      for(i in 1:ncate){
        fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
        fmu<-rbind(fmu,fmutmp)
        fgammatmp<-XY*(y-m[,i])*(z==i)
        fgamma<-cbind(fgamma,fgammatmp)
        faugztmp<-(m[,i]-augz[i])*wtt[,i]*(z==i)
        faughtmp<-ftilt(e)*(m[,i]-augh[i])
        faugz<-rbind(faugz,faugztmp)
        faugh<-rbind(faugh,faughtmp)
      }
      for(i in 2:ncate){
        fbetatmp<-((z==i)-e[,i])*W
        fbeta<-cbind(fbeta,fbetatmp)
      }

      f<-rbind(fmu,faugz,faugh,t(fbeta),t(fgamma))
      return(f)
    }
    for(i in 1:ncate){
      a[i,c(i,2*ncate+i)]<-1
      a[i,ncate+i]<--1
    }
  }else if(type=='ecac'){
    phi<-function(theta){
      mu<-theta[1:ncate]
      augz<-theta[(ncate+1):(2*ncate)]
      augh<-theta[(2*ncate+1):(3*ncate)]
      htilt<-ftilt(eest)
      wtt<-(1/eest)*htilt

      fmu<-c()
      faugz<-c()
      faugh<-c()

      for(i in 1:ncate){
        fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
        fmu<-rbind(fmu,fmutmp)
        faugztmp<-(mest[,i]-augz[i])*wtt[,i]*(z==i)
        faughtmp<-htilt*(mest[,i]-augh[i])
        faugz<-rbind(faugz,faugztmp)
        faugh<-rbind(faugh,faughtmp)
      }

      f<-rbind(fmu,faugz,faugh)
      return(f)
    }
    for(i in 1:ncate){
      a[i,c(i,2*ncate+i)]<-1
      a[i,ncate+i]<--1
    }
  }else if(type=='eca'){
    phi<-function(theta){
      mu<-theta[1:ncate]
      augz<-theta[(ncate+1):(2*ncate)]
      augh<-theta[(2*ncate+1):(3*ncate)]
      gamma<-theta[(3*ncate+1):(3*ncate+ncate*q)]

      htilt<-ftilt(eest)
      wtt<-(1/eest)*htilt

      m<-c()

      if (family=='gaussian'){
        for(i in 1:ncate){
          mtmp<-c(XY%*%gamma[((i-1)*q+1):(i*q)])
          m<-cbind(m,mtmp)}
      }else if(family=='binomial'){
        for(i in 1:ncate){
          mtmp<-plogis(c(XY%*%gamma[((i-1)*q+1):(i*q)]))
          m<-cbind(m,mtmp)}
      }else if(family=='poisson'){
        for(i in 1:ncate){
          mtmp<-exp(c(XY%*%gamma[((i-1)*q+1):(i*q)]))*offset.e
          m<-cbind(m,mtmp)}
      }

      fmu<-c()
      fgamma<-c()
      faugz<-c()
      faugh<-c()

      for(i in 1:ncate){
        fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
        fmu<-rbind(fmu,fmutmp)
        fgammatmp<-XY*(y-m[,i])*(z==i)
        fgamma<-cbind(fgamma,fgammatmp)
        faugztmp<-(m[,i]-augz[i])*wtt[,i]*(z==i)
        faughtmp<-htilt*(m[,i]-augh[i])
        faugz<-rbind(faugz,faugztmp)
        faugh<-rbind(faugh,faughtmp)
      }

      f<-rbind(fmu,faugz,faugh,t(fgamma))
      return(f)
    }

    for(i in 1:ncate){
      a[i,c(i,2*ncate+i)]<-1
      a[i,ncate+i]<--1
    }
  }else if(type=="eac"){
    phi<-function(theta){
      mu<-theta[1:ncate]
      augz<-theta[(ncate+1):(2*ncate)]
      augh<-theta[(2*ncate+1):(3*ncate)]
      beta<-theta[(3*ncate+1):(3*ncate+(ncate-1)*p)]

      e<-rep(1,n)
      for(i in 1:(ncate-1)){
        etmp<-exp(W%*%beta[((i-1)*p+1):(i*p)])
        e<-cbind(e,etmp)
      }
      e<-e/(apply(e,1,sum))
      htilt<-ftilt(e)
      wtt<-(1/e)*htilt

      fmu<-c()
      fbeta<-c()
      faugz<-c()
      faugh<-c()

      for(i in 1:ncate){
        fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
        fmu<-rbind(fmu,fmutmp)
        faugztmp<-(mest[,i]-augz[i])*wtt[,i]*(z==i)
        faughtmp<-htilt*(mest[,i]-augh[i])
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
    for(i in 1:ncate){
      a[i,c(i,2*ncate+i)]<-1
      a[i,ncate+i]<--1
    }
  }

  mphi<-function(theta){
    rowMeans(phi(theta))
  }

  #define the meat B, covariance operator
  Omega<-function(theta){
    phis<-phi(theta)
    return(tcrossprod(phis)/n)
  }

  Atheta<-jacobian(mphi,thetaest)
  invAtheta <- MASS::ginv(Atheta)

  #calculate the sandwich variance
  Vtmp<-invAtheta%*%Omega(thetaest)%*%t(invAtheta)/n

  covmu<-a%*%Vtmp%*%t(a)

  return(covmu)
}
