# Multiple ATO with ps value supplied



ATOmul_p<-function(ps.formula,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=200,out.formula=NULL,out.estimate=NULL,family=NULL,delta=0){


  #extract y
  if(typeof(out.formula)!="character"){
    out.formula=Reduce(paste0,deparse(out.formula))
  }

  y=unlist(data[yname])



  #set ordered group
  categoryz<-unique(unlist(data[zname]))
  z<-as.numeric(factor(unlist(data[zname])))
  oldlevel<-categoryz[order(unique(z))]


  #set ordered group
  ncate<-length(categoryz)

  if(sum(colnames(ps.estimate)%in%categoryz)<ncate){
    z<-as.numeric(factor(unlist(data[zname])))
    oldlevel<-categoryz[order(unique(z))]
    ps.estimate<-as.matrix(ps.estimate)
    colnames(ps.estimate)<-oldlevel
    e.h<-ps.estimate
    warning("wrong column name set for ps.estimate, treatment set as: ",oldlevel)
  }else{
    z<-as.numeric(factor(unlist(data[zname])))
    oldlevel<-categoryz[order(unique(z))]
    ps.estimate<-as.matrix(ps.estimate)
    ps.estimate<-ps.estimate[,match(oldlevel,colnames(ps.estimate))]
    e.h<-ps.estimate
  }



  #dataframe used for propensity model excluding yname
  data1<-data[,!names(data)%in%yname]
  data1[zname]<-z

  # summary statistics
  n<-length(y)


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

  #no augmentation, no bootstrap
  if(augmentation==FALSE){
    #define the m estimator
    phi<-function(theta){
      mu<-theta[1:ncate]
      wtt<-(1/e.h)*(1/apply(1/e.h,1,sum))
      fmu<-c()
      fbeta<-c()
      for(i in 1:ncate){
        fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
        fmu<-rbind(fmu,fmutmp)
      }
      f<-rbind(fmu)
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

    theta.h<-c(mu.h)
    Atheta<-jacobian(mphi,theta.h)
    invAtheta <- solve(Atheta)
    V<-invAtheta%*%Omega(theta.h)%*%t(invAtheta)/n
    covmu<-V[1:ncate,1:ncate]
    muhat<-mu.h
    names(muhat)<-oldlevel
    colnames(covmu)<-rownames(covmu)<-oldlevel
    muboot<-NULL
  }

  #augmentation and no bootstrap
  if(augmentation==TRUE){
    offset.e<-rep(1,n)
    #compute the augmented tau and variance
    #fit outcome regression model for different treatment groups
    if(is.null(out.estimate)){
      dataaug<-data[,colnames(data)!=zname]
      XY<-model.matrix(formula(out.formula),data=dataaug)
      gamma.h<-c()
      if(family=='gaussian'){
        m.h<-c()
        for(i in 1:ncate){
          dataaugtmp<-dataaug[z==i,]
          outcomefittmp<-lm(out.formula,data=dataaug[z==i,])
          gamma.tmp<-as.numeric(coef(outcomefittmp))
          gamma.h<-c(gamma.h,gamma.tmp)
          #predict using outcome regression for two treatent groups
          m.h<-cbind(m.h,XY%*%gamma.tmp)
        }
      }else if(family=='binomial'){
        m.h<-c()
        for(i in 1:ncate){
          dataaugtmp<-dataaug[z==i,]
          outcomefittmp<-glm(out.formula,data=dataaug[z==i,],family =binomial(link = "logit"))
          gamma.tmp<-as.numeric(coef(outcomefittmp))
          gamma.h<-c(gamma.h,gamma.tmp)
          #predict using outcome regression for two treatent groups
          m.h<-cbind(m.h,plogis(XY%*%gamma.tmp))
        }
      }else{
        m.h<-c()
        for(i in 1:ncate){
          dataaugtmp<-dataaug[z==i,]
          outcomefittmp<-glm(out.formula,data=dataaug[z==i,],family = poisson())
          gamma.tmp<-as.numeric(coef(outcomefittmp))
          gamma.h<-c(gamma.h,gamma.tmp)
          #predict using outcome regression for two treatent groups
          m.h<-cbind(m.h,predict(outcomefittmp,dataaug,type='response'))
        }
        if(sum(grepl('offset',names(outcomefittmp$model)))>0){
          outcomefitall<-glm(out.formula,family = poisson(),data = dataaug)
          offset.e<-apply(outcomefitall$model[,grepl('offset',names(outcomefitall$model)),drop=FALSE],1,function(x) exp(sum(x)))
        }
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
      q<-ncol(XY)
      #Estimate the sandwich variance after augmentation
      #define the m estimator
      phi<-function(theta){
        mu<-theta[1:ncate]
        augz<-theta[(ncate+1):(2*ncate)]
        augh<-theta[(2*ncate+1):(3*ncate)]
        gamma<-theta[(3*ncate+1):(3*ncate+ncate*q)]

        htilt<-(1/apply(1/e.h,1,sum))
        wtt<-(1/e.h)*htilt

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
        theta.h<-c(mu.h,augz.h,augh.h,gamma.h)
        Atheta<-jacobian(mphi,theta.h)
        invAtheta <- solve(Atheta)
        conser<-0
      },error = function(w) {
        warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
       })

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
    }else{
        dataaug<-data[,colnames(data)!=zname]
        if(sum(colnames(out.estimate)%in%categoryz)<ncate){
          out.estimate<-as.matrix(out.estimate)
          m.h<-out.estimate
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
        #define the m estimator
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
        #define score function
        mphi<-function(theta){
          rowMeans(phi(theta))
        }

        #define the meat B, covariance operator
        Omega<-function(theta){
          phis<-phi(theta)
          return(tcrossprod(phis)/n)
        }

        theta.h<-c(mu.h,augz.h,augh.h)
        Atheta<-jacobian(mphi,theta.h)
        invAtheta <- solve(Atheta)
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
      }
    }


  colnames(e.h)<-oldlevel
  out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
  class(out)<-'PSweight'
  out
}






