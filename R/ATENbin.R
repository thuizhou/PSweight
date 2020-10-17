# Binary ATEN(entropy) with ps and outcome model specified


# module for calculating the estimates
ATENbin <- function(ps.formula,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=200,out.formula=NULL,out.estimate=NULL,family=NULL,delta=0){

  #extract z name
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
  z<-as.numeric(factor(unlist(data[zname])))-1
  oldlevel<-categoryz[order(unique(z))]

  data[zname]<-z
  data1<-data[,colnames(data)!=yname]


  # summary statistics
  n1 <- sum(z)
  n0 <- sum(1-z)
  n <- n0 + n1


  #estimated ps
  fit <- glm(ps.formula, family = binomial(link = "logit"),data=data1)
  e.h <- as.numeric(fit$fitted.values)
  beta.h<-as.numeric(coef(fit))
  W<- model.matrix(fit)


  ##calculate tilting function, clip propensity to certain range to stablize
  e.hclip<-as.matrix(pmax(cbind(e.h,1-e.h),1e-6))
  htilt.h<-c(-rowSums(log(e.hclip)*e.hclip))

  # point estimate
  mu1.h <- sum(z*y*htilt.h/e.h) / sum(z*htilt.h/e.h)
  mu0.h <- sum((1-z)*y*htilt.h/(1-e.h)) / sum((1-z)*htilt.h/(1-e.h))


  #no augmentation, no bootstrap compute the variance by using numeric gradient
  if(augmentation==FALSE&bootstrap==FALSE){
      p<-ncol(W)
      #define the m estimator
      phi<-function(theta){
        mu1<-theta[1]
        mu0<-theta[2]
        beta<-theta[3:(p+2)]
        e<-plogis(c(W%*%beta))
        htilt<--((log(e)*e)+log(1-e)*(1-e))
        f1<-(y-mu1)*htilt*z/e
        f2<-(y-mu0)*htilt*(1-z)/(1-e)
        f3<-W*(z-e)
        f<-rbind(f1,f2,t(f3))
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
        theta.h<-c(mu1.h,mu0.h,beta.h)
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
          mu1<-theta[1]
          mu0<-theta[2]

          f1<-(y-mu1)*htilt.h*z/e.h
          f2<-(y-mu0)*htilt.h*(1-z)/(1-e.h)
          f<-rbind(f1,f2)
          return(f)
        }

        theta.h<-c(mu1.h,mu0.h)
        Atheta<-jacobian(mphi,theta.h)
        invAtheta <- solve(Atheta)
      }

      V<-invAtheta%*%Omega(theta.h)%*%t(invAtheta)/n
      covmu<-V[2:1,2:1]
      muhat<-c(mu0.h,mu1.h)
      names(muhat)<-oldlevel
      colnames(covmu)<-rownames(covmu)<-oldlevel
      muboot<-NULL
    }

  #no augmentation, with bootstrap
  if(augmentation==FALSE&bootstrap==TRUE){
      mu1.boot<-c()
      mu0.boot<-c()
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }
        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data1[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]
        fit.b <- glm(ps.formula, family = binomial(link = "logit"),data=data.b)
        e.b <- as.numeric(fit.b$fitted.values)

        e.clip.b<-as.matrix(pmax(cbind(e.b,1-e.b),1e-6))
        htilt.b<-c(-rowSums(log(e.clip.b)*e.clip.b))


        # point estimate
        mu1.b <- sum(z.b*htilt.b*y.b/e.b) / sum(z.b*htilt.b/e.b)
        mu0.b <- sum((1-z.b)*y.b*htilt.b/(1-e.b)) / sum((1-z.b)*htilt.b/(1-e.b))
        mu1.boot<-c(mu1.boot,mu1.b)
        mu0.boot<-c(mu0.boot,mu0.b)
      }
      muhat<-c(mu0.h,mu1.h)
      covmu<-cov(cbind(mu0.boot,mu1.boot))
      muboot<-cbind(mu0.boot,mu1.boot)
      names(muhat)<-oldlevel
      colnames(covmu)<-rownames(covmu)<-oldlevel
      colnames(muboot)<-oldlevel
    }

  #augmentation
  if(augmentation==TRUE){
      #fit two outcome regression model for different treatment groups
      dataaug<-data[,colnames(data)!=zname]
      dataaug0<-dataaug[z==0,]
      dataaug1<-dataaug[z==1,]
      offset.e<-rep(1,n)
      if(family=='gaussian'){
        outcomefit1<-lm(out.formula,data = dataaug1)
        outcomefit0<-lm(out.formula,data=dataaug0)
        XY<-model.matrix(formula(out.formula),data=dataaug)
        gamma0.h<-as.numeric(coef(outcomefit0))
        gamma1.h<-as.numeric(coef(outcomefit1))
        #predict using outcome regression for two treatent groups
        m0.h<-c(XY%*%gamma0.h)
        m1.h<-c(XY%*%gamma1.h)
      }else if(family=='binomial'){
        outcomefit1<-glm(out.formula,family = binomial(link = "logit"),data = dataaug1)
        outcomefit0<-glm(out.formula,family = binomial(link = "logit"),data=dataaug0)
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
        m0.h<-c(predict(outcomefit0,dataaug,type='response'))
        m1.h<-c(predict(outcomefit1,dataaug,type='response'))
        #if offset exits
        if(sum(grepl('offset',names(outcomefit0$model)))>0){
          outcomefitall<-glm(out.formula,family = poisson(),data = dataaug)
          offset.e<-apply(outcomefitall$model[,grepl('offset',names(outcomefitall$model)),drop=FALSE],1,function(x) exp(sum(x)))
        }
      }

      #calculate the aumentation term
      aug0h.h<-sum(m0.h*htilt.h)/sum(htilt.h)
      aug0z.h<-sum((1-z)*htilt.h*m0.h/(1-e.h))/sum((1-z)*htilt.h/(1-e.h))
      aug1h.h<-sum(m1.h*htilt.h)/sum(htilt.h)
      aug1z.h<-sum(z*htilt.h*m1.h/e.h)/sum(z*htilt.h/e.h)

      mu0aug.h<-mu0.h+aug0h.h-aug0z.h
      mu1aug.h<-mu1.h+aug1h.h-aug1z.h

      if(bootstrap==FALSE){
        #Estimate the sandwich variance after augmentation
        p<-ncol(W)
        q<-ncol(XY)
        #define the m estimator
        phi<-function(theta){
          mu1<-theta[1]
          aug1h<-theta[2]
          aug1z<-theta[3]
          mu0<-theta[4]
          aug0h<-theta[5]
          aug0z<-theta[6]
          beta<-theta[7:(p+6)]
          gamma1<-theta[(p+7):(p+6+q)]
          gamma0<-theta[(p+7+q):(p+2*q+6)]
          e<-plogis(c(W%*%beta))
          htilt<--((log(e)*e)+log(1-e)*(1-e))
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
          f1<-z*htilt/e*(y-mu1)
          f2<-(m1-aug1h)*htilt
          f3<-z*htilt/e*(m1-aug1z)
          f4<-(1-z)*htilt/(1-e)*(y-mu0)
          f5<-(m0-aug0h)*htilt
          f6<-(1-z)*htilt/(1-e)*(m0-aug0z)
          f7<-XY*(y-m1)*z
          f8<-XY*(y-m0)*(1-z)
          f9<-W*(z-e)

          f<-rbind(f1,f2,f3,f4,f5,f6,t(f7),t(f8),t(f9))
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
          theta.h<-c(mu1.h,aug1h.h,aug1z.h,mu0.h,aug0h.h,aug0z.h,beta.h,gamma1.h,gamma0.h)
          Atheta<-jacobian(mphi,theta.h)
          invAtheta <- solve(Atheta)
          a0<-c(0,0,0,1,1,-1,rep(0,p+2*q))
          a1<-c(1,1,-1,0,0,0,rep(0,p+2*q))
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
            f1<-z*htilt.h/e.h*(y-mu1)
            f2<-htilt.h*(m1.h-aug1h)
            f3<-z*htilt.h/e.h*(m1.h-aug1z)
            f4<-(1-z)*htilt.h/(1-e.h)*(y-mu0)
            f5<-htilt.h*(m0.h-aug0h)
            f6<-(1-z)*htilt.h/(1-e.h)*(m0.h-aug0z)
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
      }else{

        #augmentation and bootstrap
        dataaug<-data[,colnames(data)!=zname]
        mu1aug.boot<-c()
        mu0aug.boot<-c()
        if(family=='gaussian'){
          for(i in 1:R){
            if(i %% 50==0){
              message("bootstrap ", i, " samples")
            }
            # estimate ps
            samp.b<-sample(n,n,replace = TRUE)
            data.b<-data1[samp.b,]
            y.b<-y[samp.b]
            z.b<-z[samp.b]
            fit.b <- glm(ps.formula, family = binomial(link = "logit"),data=data.b)
            e.b <- as.numeric(fit.b$fitted.values)

            e.clip.b<-as.matrix(pmax(cbind(e.b,1-e.b),1e-6))
            htilt.b<-c(-rowSums(log(e.clip.b)*e.clip.b))

            # point estimate
            mu1.b <- sum(z.b*htilt.b*y.b/e.b) / sum(z.b*htilt.b/e.b)
            mu0.b <- sum((1-z.b)*y.b*htilt.b/(1-e.b)) / sum((1-z.b)*htilt.b/(1-e.b))

            #calculate the augmentation term
            dataaug.b<-dataaug[samp.b,]
            dataaug0.b<-dataaug.b[z.b==0,]
            dataaug1.b<-dataaug.b[z.b==1,]
            outcomefit0.b<-lm(out.formula,data = dataaug0.b)
            outcomefit1.b<-lm(out.formula,data = dataaug1.b)
            #predict using outcome regression for two treatent groups
            m0.b<-predict(outcomefit0.b,dataaug.b)
            m1.b<-predict(outcomefit1.b,dataaug.b)

            aug0h.b<-sum(m0.b*htilt.b)/sum(htilt.b)
            aug0z.b<-sum((1-z.b)*m0.b*htilt.b/(1-e.b))/sum((1-z.b)*htilt.b/(1-e.b))
            aug1h.b<-sum(m1.b*htilt.b)/sum(htilt.b)
            aug1z.b<-sum(z.b*m1.b*htilt.b/e.b)/sum(z.b*htilt.b/e.b)
            mu1aug.boot<-c(mu1aug.boot,mu1.b+aug1h.b-aug1z.b)
            mu0aug.boot<-c(mu0aug.boot,mu0.b+aug0h.b-aug0z.b)
          }
        }else if(family=='binomial'){
          for(i in 1:R){
            if(i %% 50==0){
              message("bootstrap", i, "samples")
            }
            # estimate ps
            samp.b<-sample(n,n,replace = TRUE)
            data.b<-data1[samp.b,]
            y.b<-y[samp.b]
            z.b<-z[samp.b]
            fit.b <- glm(ps.formula, family = binomial(link = "logit"),data=data.b)
            e.b <- as.numeric(fit.b$fitted.values)

            e.clip.b<-pmax(cbind(e.b,1-e.b),1e-6)
            htilt.b<-c(-colSums(log(e.clip.b)*e.clip.b))

            # point estimate
            mu1.b <- sum(z.b*htilt.b*y.b/e.b) / sum(z.b*htilt.b/e.b)
            mu0.b <- sum((1-z.b)*y.b*htilt.b/(1-e.b)) / sum((1-z.b)*htilt.b/(1-e.b))

            #calculate the augmentation term
            dataaug.b<-dataaug[samp.b,]
            dataaug0.b<-dataaug.b[z.b==0,]
            dataaug1.b<-dataaug.b[z.b==1,]
            outcomefit0.b<-glm(out.formula,family = binomial(link = "logit"),data = dataaug0.b)
            outcomefit1.b<-glm(out.formula,family = binomial(link = "logit"),data = dataaug1.b)
            #predict using outcome regression for two treatent groups
            m0.b<-predict(outcomefit0.b,type = "response",dataaug.b)
            m1.b<-predict(outcomefit1.b,type = "response",dataaug.b)

            aug0h.b<-sum(m0.b*htilt.b)/sum(htilt.b)
            aug0z.b<-sum((1-z.b)*m0.b*htilt.b/(1-e.b))/sum((1-z.b)*htilt.b/(1-e.b))
            aug1h.b<-sum(m1.b*htilt.b)/sum(htilt.b)
            aug1z.b<-sum(z.b*m1.b*htilt.b/e.b)/sum(z.b*htilt.b/e.b)
            mu1aug.boot<-c(mu1aug.boot,mu1.b+aug1h.b-aug1z.b)
            mu0aug.boot<-c(mu0aug.boot,mu0.b+aug0h.b-aug0z.b)
          }
        }else{
          for(i in 1:R){
            if(i %% 50==0){
              message("bootstrap ", i, " samples")
            }
            # estimate ps
            samp.b<-sample(n,n,replace = TRUE)
            data.b<-data1[samp.b,]
            y.b<-y[samp.b]
            z.b<-z[samp.b]
            fit.b <- glm(ps.formula, family = binomial(link = "logit"),data=data.b)
            e.b <- as.numeric(fit.b$fitted.values)

            e.clip.b<-pmax(cbind(e.b,1-e.b),1e-6)
            htilt.b<-c(-colSums(log(e.clip.b)*e.clip.b))

            # point estimate
            mu1.b <- sum(z.b*htilt.b*y.b/e.b) / sum(z.b*htilt.b/e.b)
            mu0.b <- sum((1-z.b)*y.b*htilt.b/(1-e.b)) / sum((1-z.b)*htilt.b/(1-e.b))

            #calculate the augmentation term
            dataaug.b<-dataaug[samp.b,]
            dataaug0.b<-dataaug.b[z.b==0,]
            dataaug1.b<-dataaug.b[z.b==1,]
            outcomefit0.b<-glm(out.formula,family = poisson(),data = dataaug0.b)
            outcomefit1.b<-glm(out.formula,family = poisson(),data = dataaug1.b)

            #predict using outcome regression for two treatent groups
            m0.b<-predict(outcomefit0.b,dataaug.b,type='response')
            m1.b<-predict(outcomefit1.b,dataaug.b,type='response')

            aug0h.b<-sum(m0.b*htilt.b)/sum(htilt.b)
            aug0z.b<-sum((1-z.b)*m0.b*htilt.b/(1-e.b))/sum((1-z.b)*htilt.b/(1-e.b))
            aug1h.b<-sum(m1.b*htilt.b)/sum(htilt.b)
            aug1z.b<-sum(z.b*m1.b*htilt.b/e.b)/sum(z.b*htilt.b/e.b)
            mu1aug.boot<-c(mu1aug.boot,mu1.b+aug1h.b-aug1z.b)
            mu0aug.boot<-c(mu0aug.boot,mu0.b+aug0h.b-aug0z.b)
          }
        }
        muhat<-c(mu0aug.h, mu1aug.h)
        covmu<-cov(cbind(mu0aug.boot,mu1aug.boot))
        muboot<-cbind(mu0aug.boot,mu1aug.boot)
        names(muhat)<-oldlevel
        colnames(covmu)<-rownames(covmu)<-oldlevel
        colnames(muboot)<-oldlevel
      }
    }

    e.h<-cbind(1-e.h,e.h)
    colnames(e.h)<-oldlevel
    out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
    class(out)<-'PSweight'
    out
}





