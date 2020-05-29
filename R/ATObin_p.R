# Binary ATO with ps value supplied



ATObin_p<- function(ps.formula,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=200,out.formula=NULL,out.estimate=NULL,family=NULL,delta=0){



  #extract y
  if(typeof(out.formula)!="character"){
    out.formula=Reduce(paste0,deparse(out.formula))
  }
  y=unlist(data[yname])



  #set ordered group
  categoryz<-unique(unlist(data[zname]))
  if(sum(colnames(ps.estimate)%in%categoryz)<2){
    z<-as.numeric(factor(unlist(data[zname])))-1
    oldlevel<-categoryz[order(unique(z))]
    ps.estimate<-as.matrix(ps.estimate)
    colnames(ps.estimate)<-oldlevel
    e.h<-ps.estimate[,2]
    warning("wrong column name set for ps.estimate, treatment set as: ",oldlevel)
  }else{
    z<-as.numeric(factor(unlist(data[zname])))-1
    oldlevel<-categoryz[order(unique(z))]
    ps.estimate<-as.matrix(ps.estimate)
    ps.estimate<-ps.estimate[,match(oldlevel,colnames(ps.estimate))]
    e.h<-ps.estimate[,2]
  }


  data[zname]<-z
  data1<-data[,colnames(data)!=yname]


  # summary statistics
  n1 <- sum(z)             # number of treated
  n0 <- sum(1-z)           # number of untreated
  n <- n0 + n1  # total sample size


  # point estimate
  mu1.h <- sum(z*y*(1-e.h)) / sum(z*(1-e.h))
  mu0.h <- sum((1-z)*y*e.h) / sum((1-z)*e.h)


  #Compute the variance by using numeric gradient
  if(augmentation==FALSE){
    #define the m estimator
    phi<-function(theta){
      mu1<-theta[1]
      mu0<-theta[2]
      f1<-(y-mu1)*z*(1-e.h)
      f2<-(y-mu0)*(1-z)*e.h
      f<-rbind(f1,f2)
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

    theta.h<-c(mu1.h,mu0.h)
    Atheta<-jacobian(mphi,theta.h)
    invAtheta <- solve(Atheta)

    V<-invAtheta%*%Omega(theta.h)%*%t(invAtheta)/n
    covmu<-V[2:1,2:1]
    muhat<-c(mu0.h,mu1.h)
    names(muhat)<-oldlevel
    colnames(covmu)<-rownames(covmu)<-oldlevel
    muboot<-NULL
  }


  if(augmentation==TRUE){
    offset.e<-rep(1,n)
    if(is.null(out.estimate)){
      dataaug<-data[,colnames(data)!=zname]
      #fit two outcome regression model for different treatment groups
      dataaug0<-dataaug[z==0,]
      dataaug1<-dataaug[z==1,]
      if(family=='gaussian'){
        outcomefit0<-lm(out.formula,data=dataaug0)
        outcomefit1<-lm(out.formula,data=dataaug1)
        XY<-model.matrix(formula(out.formula),data=dataaug)
        gamma0.h<-as.numeric(coef(outcomefit0))
        gamma1.h<-as.numeric(coef(outcomefit1))
        #predict using outcome regression for two treatent groups
        m0.h<-c(XY%*%gamma0.h)
        m1.h<-c(XY%*%gamma1.h)
      }else if(family=='binomial'){
        outcomefit0<-glm(out.formula,family = binomial(link = "logit"),data=dataaug0)
        outcomefit1<-glm(out.formula,family = binomial(link = "logit"),data = dataaug1)
        XY<-model.matrix(formula(out.formula),data=dataaug)
        gamma0.h<-as.numeric(coef(outcomefit0))
        gamma1.h<-as.numeric(coef(outcomefit1))
        #predict using outcome regression for two treatent groups
        m0.h<-c(plogis(XY%*%gamma0.h))
        m1.h<-c(plogis(XY%*%gamma1.h))
      }else{
        outcomefit1<-glm(out.formula,family = poisson(),data = dataaug1)
        outcomefit0<-glm(out.formula,family = poisson(),data=dataaug0)
        XY<-model.matrix(formula(out.formula),data=dataaug)
        gamma0.h<-as.numeric(coef(outcomefit0))
        gamma1.h<-as.numeric(coef(outcomefit1))
        #predict using outcome regression for two treatent groups
        m0.h<-predict(outcomefit0,dataaug,type='response')
        m1.h<-predict(outcomefit1,dataaug,type='response')
        #if offset exits
        if(sum(grepl('offset',names(outcomefit0$model)))>0){
          outcomefitall<-glm(out.formula,family = poisson(),data = dataaug)
          offset.e<-apply(outcomefitall$model[,grepl('offset',names(outcomefitall$model)),drop=FALSE],1,function(x) exp(sum(x)))
        }
      }

      #calculate the augmentation term and updata tau
      aug0h.h<-sum(e.h*(1-e.h)*m0.h)/sum(e.h*(1-e.h))
      aug0z.h<-sum((1-z)*m0.h*e.h)/sum((1-z)*e.h)
      aug1h.h<-sum(e.h*(1-e.h)*m1.h)/sum(e.h*(1-e.h))
      aug1z.h<-sum(z*m1.h*(1-e.h))/sum(z*(1-e.h))

      mu0aug.h<-mu0.h+aug0h.h-aug0z.h
      mu1aug.h<-mu1.h+aug1h.h-aug1z.h

      #Estimate the sandwich variance after augmentation
      q<-ncol(XY)
      #define the m estimator
      phi<-function(theta){
        mu1<-theta[1]
        aug1h<-theta[2]
        aug1z<-theta[3]
        mu0<-theta[4]
        aug0h<-theta[5]
        aug0z<-theta[6]
        gamma1<-theta[7:(6+q)]
        gamma0<-theta[(7+q):(2*q+6)]
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

        #component in m estimator
        f1<-(1-e.h)*z*(y-mu1)
        f2<-(1-e.h)*e.h*(m1-aug1h)
        f3<-(1-e.h)*z*(m1-aug1z)
        f4<-e.h*(1-z)*(y-mu0)
        f5<-(1-e.h)*e.h*(m0-aug0h)
        f6<-e.h*(1-z)*(m0-aug0z)
        f7<-XY*(y-m1)*z
        f8<-XY*(y-m0)*(1-z)

        f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8))
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
        theta.h<-c(mu1.h,aug1h.h,aug1z.h,mu0.h,aug0h.h,aug0z.h,gamma1.h,gamma0.h)
        Atheta<-jacobian(mphi,theta.h)
        invAtheta <- solve(Atheta)
        a0<-c(0,0,0,1,1,-1,rep(0,2*q))
        a1<-c(1,1,-1,0,0,0,rep(0,2*q))
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


      a<-rbind(a0,a1)
      V1<-invAtheta%*%Omega(theta.h)%*%t(invAtheta)/n
      covmu<-a%*%V1%*%t(a)
      muhat<-c(mu0aug.h, mu1aug.h)
      muboot<-NULL
      names(muhat)<-oldlevel
      colnames(covmu)<-rownames(covmu)<-oldlevel

    }else{
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
      invAtheta <- ginv(Atheta)
      tryCatch( { invAtheta <- solve(Atheta) }
                ,error = function(w) {
                  warning("The sandwich matrix not pd, therefore not invertable, use g-inverse instead, please double check") })

      V1<-invAtheta%*%Omega(theta.h)%*%t(invAtheta)/n
      a0<-c(0,0,0,1,1,-1)
      a1<-c(1,1,-1,0,0,0)
      a<-rbind(a0,a1)
      covmu<-a%*%V1%*%t(a)
      muhat<-c(mu0aug.h, mu1aug.h)
      muboot<-NULL
      names(muhat)<-oldlevel
      colnames(covmu)<-rownames(covmu)<-oldlevel
    }
  }



  e.h<-cbind(1-e.h,e.h)
  colnames(e.h)<-oldlevel
  out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
  class(out)<-'PSweight'
  out
}




