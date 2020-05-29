# Multiple ATE with ps and outcome model specified


# module for calculating the estimates
ATEmul <- function(ps.formula,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=200,out.formula=NULL,out.estimate=NULL,family=NULL,delta=0){

  #extract the name for treatment
  if(typeof(ps.formula)!="character"){
    ps.formula=Reduce(paste0,deparse(ps.formula))
  }
  ps.formula=gsub(" ","",ps.formula)
  zname=unlist(strsplit(ps.formula,'~'))[[1]][1]



  #extract y
  if(typeof(out.formula)!="character"){
    out.formula=Reduce(paste0,deparse(out.formula))
  }
  y=unlist(data[yname])


  #set ordered group
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


  #generate weight for table and IPW
  allwt<-1/e.h
  wt<-rep(0,dim(allwt)[1])
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
  if(augmentation==FALSE&bootstrap==FALSE){
    p<-ncol(W)
    #define the m estimator
    phi<-function(theta){
      mu<-theta[1:ncate]
      beta<-theta[(ncate+1):(ncate+(ncate-1)*p)]
      e<-rep(1,n)
      for(i in 1:(ncate-1)){
        etmp=exp(W%*%beta[((i-1)*p+1):(i*p)])
        e<-cbind(e,etmp)
      }
      e<-e/(apply(e,1,sum))
      wtt<-1/e
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
      theta.h<-c(mu.h,beta.h)
      Atheta<-jacobian(mphi,theta.h)
      invAtheta <- solve(Atheta)
      conser<-0
    },error = function(w) {
      warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
      })

    #if not pd than use conservative
    if(conser==1){
      #define the m estimator
      phi<-function(theta){
        mu<-theta[1:ncate]
        wtt<-1/e.h
        fmu<-c()
        for(i in 1:ncate){
          fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
          fmu<-rbind(fmu,fmutmp)
        }
        f<-rbind(fmu)
        return(f)
      }
      theta.h<-c(mu.h)
      Atheta<-jacobian(mphi,theta.h)
      invAtheta <- solve(Atheta)
    }

    V<-invAtheta%*%Omega(theta.h)%*%t(invAtheta)/n
    covmu<-V[1:ncate,1:ncate]
    muhat<-mu.h
    names(muhat)<-oldlevel
    colnames(covmu)<-rownames(covmu)<-oldlevel
    muboot<-NULL
  }

  #no augmentation, with bootstrap
  if(augmentation==FALSE&bootstrap==TRUE){
    muboot<-c()
    for(i in 1:R){
      if(i %% 50==0){
        message("bootstrap ", i, " samples")
      }
      # estimate ps
      samp.b<-sample(n,n,replace = TRUE)
      data.b<-data1[samp.b,]
      y.b<-y[samp.b]
      z.b<-z[samp.b]
      fit.b <- multinom(formula = ps.formula, data=data.b,maxit = 500, Hess = TRUE, trace = FALSE)
      e.b <- fit.b$fitted.values
      # point estimate
      allwt.b<-(1/e.b)
      #the point estimate
      mu.tmp<-c()
      for(i in 1:ncate){
        ytmp.b<-y.b[z.b==i]
        wttmp.b<-allwt.b[z.b==i,i]
        mu.tmp<-c(mu.tmp,sum(ytmp.b*wttmp.b)/sum(wttmp.b))
      }
      muboot<-rbind(muboot,mu.tmp)
    }
    muhat<-apply(muboot,2,mean)
    covmu<-cov(muboot)
    names(muhat)<-oldlevel
    colnames(covmu)<-rownames(covmu)<-oldlevel
    colnames(muboot)<-oldlevel
    rownames(muboot)<-NULL
  }

  #augmentation and no bootstrap
  if(augmentation==TRUE&bootstrap==FALSE){
    offset.e<-rep(1,n)
    #compute the augmented tau and variance
    #fit outcome regression model for different treatment groups
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
        #predict using outcome regression
        m.h<-cbind(m.h,XY%*%gamma.tmp)
      }
    }else if(family=='binomial'){
      m.h<-c()
      for(i in 1:ncate){
        dataaugtmp<-dataaug[z==i,]
        outcomefittmp<-glm(out.formula,data=dataaug[z==i,],family =binomial(link = "logit"))
        gamma.tmp<-as.numeric(coef(outcomefittmp))
        gamma.h<-c(gamma.h,gamma.tmp)
        #predict using outcome regression
        m.h<-cbind(m.h,plogis(XY%*%gamma.tmp))
      }
    }else{
      m.h<-c()
      for(i in 1:ncate){
        dataaugtmp<-dataaug[z==i,]
        outcomefittmp<-glm(out.formula,data=dataaug[z==i,],family = poisson())
        gamma.tmp<-as.numeric(coef(outcomefittmp))
        gamma.h<-c(gamma.h,gamma.tmp)
        #predict using outcome regression
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
      aughtmp<-mean(mtmp)
      augz.h<-c(augz.h,augztmp)
      augh.h<-c(augh.h,aughtmp)
    }

    muhat<-mu.h+augh.h-augz.h
    names(muhat)<-oldlevel

    #Estimate the sandwich variance after augmentation
    p<-ncol(W)
    q<-ncol(XY)
    #define the m estimator
    phi<-function(theta){
      mu<-theta[1:ncate]
      augz<-theta[(ncate+1):(2*ncate)]
      augh<-theta[(2*ncate+1):(3*ncate)]
      beta<-theta[(3*ncate+1):(3*ncate+(ncate-1)*p)]
      gamma<-theta[(3*ncate+(ncate-1)*p+1):(3*ncate+(ncate-1)*p+ncate*q)]

      e<-rep(1,n)
      for(i in 1:(ncate-1)){
        etmp=exp(W%*%beta[((i-1)*p+1):(i*p)])
        e<-cbind(e,etmp)
      }
      e<-e/(apply(e,1,sum))
      wtt<-(1/e)

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
        faughtmp<-(m[,i]-augh[i])
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
      theta.h<-c(mu.h,augz.h,augh.h,beta.h,gamma.h)
      Atheta<-jacobian(mphi,theta.h)
      invAtheta <- solve(Atheta)
      conser<-0
      },error = function(w) {
        warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })

    #if not pd than use conservative
    if(conser==1){
      #define the m estimator
      phi<-function(theta){
        mu<-theta[1:ncate]
        augz<-theta[(ncate+1):(2*ncate)]
        augh<-theta[(2*ncate+1):(3*ncate)]

        wtt<-(1/e.h)
        fmu<-c()
        faugz<-c()
        faugh<-c()

        for(i in 1:ncate){
          fmutmp<-(y-mu[i])*wtt[,i]*(z==i)
          fmu<-rbind(fmu,fmutmp)
          faugztmp<-(m.h[,i]-augz[i])*wtt[,i]*(z==i)
          faughtmp<-(m.h[,i]-augh[i])
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
  }

  #augmentation and  bootstrap
  if(augmentation==TRUE&bootstrap==TRUE){
    dataaug<-data[,colnames(data)!=zname]
    muboot<-c()
    if(family=='gaussian'){
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples", "\n")
        }
        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data1[samp.b,]
        dataaug.b<-dataaug[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]
        fit.b <- multinom(formula = ps.formula, data=data.b,maxit = 500, Hess = TRUE, trace = FALSE)
        e.b <- fit.b$fitted.values
        # point estimate
        allwt.b<-(1/e.b)
        #the point estimate
        mu.tmp<-c()
        augz.tmp<-c()
        augh.tmp<-c()
        for(i in 1:ncate){
          dataaug.tmp<-dataaug.b[z==i,]
          outcomefit.b<-lm(out.formula,data=dataaug.tmp)
          m.b<-predict(outcomefit.b,dataaug.b)
          mu.tmp<-c(mu.tmp,sum(y.b*allwt.b[,i]*(z.b==i))/sum(allwt.b[,i]*(z.b==i)))
          augz.tmp<-c(augz.tmp,sum(m.b*allwt.b[,i]*(z.b==i))/sum(allwt.b[,i]*(z.b==i)))
          augh.tmp<-c(augh.tmp,mean(m.b))
        }
        muboot<-rbind(muboot,mu.tmp-augz.tmp+augh.tmp)
      }
    }else if(family=='binomial'){
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }
        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data1[samp.b,]
        dataaug.b<-dataaug[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]
        fit.b <- multinom(formula = ps.formula, data=data.b,maxit = 500, Hess = TRUE, trace = FALSE)
        e.b <- fit.b$fitted.values
        # point estimate
        allwt.b<-(1/e.b)
        #the point estimate
        mu.tmp<-c()
        augz.tmp<-c()
        augh.tmp<-c()
        for(i in 1:ncate){
          dataaug.tmp<-dataaug.b[z==i,]
          outcomefit.b<-glm(out.formula,family = binomial(link = "logit"),data=dataaug.tmp)
          m.b<-predict(outcomefit.b,type = "response",dataaug.b)
          mu.tmp<-c(mu.tmp,sum(y.b*allwt.b[,i]*(z.b==i))/sum(allwt.b[,i]*(z.b==i)))
          augz.tmp<-c(augz.tmp,sum(m.b*allwt.b[,i]*(z.b==i))/sum(allwt.b[,i]*(z.b==i)))
          augh.tmp<-c(augh.tmp,mean(m.b))
        }
        muboot<-rbind(muboot,mu.tmp-augz.tmp+augh.tmp)
      }
    }else{
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }
        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data1[samp.b,]
        dataaug.b<-dataaug[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]
        fit.b <- multinom(formula = ps.formula, data=data.b,maxit = 500, Hess = TRUE, trace = FALSE)
        e.b <- fit.b$fitted.values
        # point estimate
        allwt.b<-(1/e.b)
        #the point estimate
        mu.tmp<-c()
        augz.tmp<-c()
        augh.tmp<-c()
        for(i in 1:ncate){
          dataaug.tmp<-dataaug.b[z==i,]
          outcomefit.b<-glm(out.formula,family = poisson(),data=dataaug.tmp)
          m.b<-predict(outcomefit.b,type='response',dataaug.b)
          mu.tmp<-c(mu.tmp,sum(y.b*allwt.b[,i]*(z.b==i))/sum(allwt.b[,i]*(z.b==i)))
          augz.tmp<-c(augz.tmp,sum(m.b*allwt.b[,i]*(z.b==i))/sum(allwt.b[,i]*(z.b==i)))
          augh.tmp<-c(augh.tmp,mean(m.b))
        }
        muboot<-rbind(muboot,mu.tmp-augz.tmp+augh.tmp)
      }
    }
    muhat<-apply(muboot,2,mean)
    covmu<-cov(muboot)
    names(muhat)<-oldlevel
    colnames(covmu)<-rownames(covmu)<-oldlevel
    colnames(muboot)<-oldlevel
    rownames(muboot)<-NULL
  }



  colnames(e.h)<-oldlevel
  out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
  class(out)<-'PSweight'
  out
}
