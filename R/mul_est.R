#### Point and variance for multiple treatment ############################################
###########################################################################################

mulest<-function(ps.formula=NULL,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=50,out.formula=NULL,out.estimate=NULL,family=NULL,weight="overlap",ps.method='glm',ps.control=list(),out.method='glm',out.control=list()){

  #preprocess formula and extract y
  out.formula<-as.formula(out.formula)
  y<-unlist(data[yname])


  #set ordered group
  facz<-as.factor(unlist(data[zname]))

  #set the treatment label
  oldlevel<-levels(facz)
  newlevel<-oldlevel

  z<-as.numeric(facz)
  n<-length(z) #total obs
  ncate<-length(oldlevel) #number of categories

  #set group for ATT
  if(weight=="treated"){
    if(is.null(trtgrp)){
      trtgrp<-tail(oldlevel,1)
    }else{
      newlevel<-rev(unique(c(trtgrp,oldlevel)))
      facz<-factor(unlist(data[zname]),levels = newlevel)
      z<-as.numeric(facz)
    }
  }

  matchlevel<-match(oldlevel,newlevel) #match level for ATT


  #reassign value as numeric
  data[zname]<-z
  data_p<-data[,colnames(data)!=yname]



  # obtain ps estimation
  # estimate with formula
  if(is.null(ps.estimate)){
    fit<-do.call(PSmethod, list(ps.formula=ps.formula, method=ps.method, data=data_p, ncate=ncate, ps.control=ps.control))
    W<- model.matrix(formula(ps.formula),data_p)                      # design matrix
    e.h <- fit$e.h
    beta.h<-as.numeric(fit$beta.h)
    if(ps.method!='glm'){
      ps.estimate<-fit$e.h
    }
  }else{
    #the name for the propensity score
    if(!setequal(colnames(ps.estimate),newlevel)){
      ps.estimate<-as.matrix(ps.estimate)
      e.h<-ps.estimate
      warning("wrong column name set for ps.estimate, treatment set as: ",paste(newlevel,collapse = ', '))
    }else{
      ps.estimate<-ps.estimate[,match(newlevel,colnames(ps.estimate))]
      e.h<-ps.estimate
    }
  }


  if(weight=="entropy"){
    e.h<-pmax(e.h,1e-6)
  }


  #compute outcome regression for augmentation
  if(augmentation){
    #no outcome estimation provided
    if(is.null(out.estimate)){
      offset.e<-rep(1,n)
      #fit outcome regression model for different treatment groups
      dataaug<-data[,colnames(data)!=zname]
      XY<-model.matrix(formula(out.formula),data=dataaug)
      gamma.h<-c()
      m.h<-c()
      for(i in 1:ncate){
        dataaugtmp<-dataaug[z==i,]
        outcomefittmp<-do.call(OUTmethod,list(out.formula=out.formula,y=y[z==i], out.method=out.method, family=family, datain=dataaug[z==i,], dataout=dataaug,out.control=out.control))
        gamma.tmp<-outcomefittmp$gamma.h
        gamma.h<-c(gamma.h,gamma.tmp)
        #predict using outcome regression for two treatent groups
        m.h<-cbind(m.h,outcomefittmp$m.est)
      }
      if(family=='poisson'){
        offsetlog<-model.extract(model.frame(out.formula,data = dataaug),'offset')
        if(!is.null(offsetlog)){offset.e<-exp(offsetlog)}
      }

      if(out.method!='glm'){
        out.estimate<-m.h
        colnames(out.estimate)<-newlevel
      }

    }else{
      #the name for the outcome regression
      if(!setequal(colnames(out.estimate),newlevel)){
        out.estimate<-as.matrix(out.estimate)
        m.h<-out.estimate
        warning("wrong column name set for out.estimate, treatment set as: ",paste(newlevel,collapse = ', '))
      }else{
        out.estimate<-out.estimate[,match(newlevel,colnames(out.estimate))]
        m.h<-out.estimate
      }
    }
  }


  #tilting function
  ftilt<-tiltmul(weight = weight)


  ##No Augmentation###############################################################
  if(augmentation==FALSE){
    muhat<-ptmul(e.h,z,y,ftilt,ncate,n)
    ##No bootstrap###############################################################
    if(bootstrap==FALSE){

      #choose conservative or not
      conser<-1

      #use ps formula
      if(is.null(ps.estimate)){
        tryCatch( {
          theta.h<-c(muhat,beta.h)
          covmu<-sand_mul(z,y,n,ncate,ftilt,theta.h,W,type='e')
          conser<-0
          },error = function(w) {
            warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
             })
      }

      #use conservative
      if(conser==1){
        theta.h<-muhat
        covmu<-sand_mul(z,y,n,ncate,ftilt,theta.h,eest=e.h,type='ec')
      }

      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
      muboot<-NULL
    }else{
      ##With bootstrap###############################################################

      muboot<-c()
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }
        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data_p[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]

        fit.b<-do.call(PSmethod,list(ps.formula = ps.formula, method=ps.method, data=data.b,ncate=ncate,ps.control=ps.control))
        e.b <- fit.b$e.h

        mu.b<-ptmul(e.b,z.b,y.b,ftilt,ncate,n)

        muboot<-rbind(muboot,mu.b)
      }
      covmu<-cov(muboot)
      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
      colnames(muboot)<-newlevel
      rownames(muboot)<-NULL
    }
  }

  ##Augmentation###############################################################
  if(augmentation==TRUE){

    #calculate point estimate
    augest <- ptmul(e.h,z,y,ftilt,ncate,n,m.h)
    muhat <- tail(augest,ncate)

    ##No bootstrap###############################################################
    if(bootstrap==FALSE){

      #choose conservative or not
      conser<-1

      if(is.null(ps.estimate) & is.null(out.estimate)){
        #both formula
        tryCatch( {
          theta.h<-c(augest[1:(length(augest)-ncate)],beta.h,gamma.h)
          covmu<-sand_mul(z,y,n,ncate,ftilt,theta.h,W,XY,family=family,offset.e=offset.e,type='ea')
          conser<-0
          },error = function(w) {
            warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
           })

      }else if(!is.null(ps.estimate) & is.null(out.estimate)){
        ## formula on outcome regression only
        tryCatch( {
          theta.h<-c(augest[1:(length(augest)-ncate)],gamma.h)
          covmu<-sand_mul(z,y,n,ncate,ftilt,theta.h,XY=XY,eest=e.h,family=family,offset.e=offset.e,type='eca')
          conser<-0
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })

      }else if(is.null(ps.estimate) & !is.null(out.estimate)){
        ## formula on ps only
        tryCatch( {
          theta.h<-c(augest[1:(length(augest)-ncate)],beta.h)
          covmu<-sand_mul(z,y,n,ncate,ftilt,theta.h,W,mest=m.h,type='eac')
          conser<-0
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })

      }

      #if not pd than use conservative
      if(conser==1){
        theta.h<-augest[1:(length(augest)-ncate)]
        covmu<-sand_mul(z,y,n,ncate,ftilt,theta.h,eest=e.h,mest=m.h,type='ecac')

      }

      names(muhat)<-newlevel
      muboot<-NULL
      colnames(covmu)<-rownames(covmu)<-newlevel
    }else{
      #augmentation and  bootstrap
      dataaug<-data[,colnames(data)!=zname]
      muboot<-NULL

      #bootstrap runs
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }
        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data_p[samp.b,]
        dataaug.b<-dataaug[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]
        fit.b<-do.call(PSmethod,list(ps.formula = ps.formula, method=ps.method, data=data.b,ncate=ncate,ps.control=ps.control))
        e.b <- fit.b$e.h

        m.b<-NULL
        for(i in 1:ncate){
          dataaug.tmp<-dataaug.b[z.b==i,]

          dataaugtmp<-dataaug[z==i,]
          outcomefit.b<-do.call(OUTmethod,list(out.formula=out.formula,y=y.b[z.b==i], out.method=out.method, family=family, datain=dataaug.tmp, dataout=dataaug.b,out.control=out.control))
          #predict using outcome regression for two treatent groups
          m.b<-cbind(m.b,outcomefit.b$m.est)
        }
        mu.b <- tail(ptmul(e.b,z.b,y.b,ftilt,ncate,n,m.b),ncate)
        muboot<-rbind(muboot,mu.b)
      }

      covmu<-cov(muboot)
      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
      colnames(muboot)<-newlevel
      rownames(muboot)<-NULL
    }
  }

  colnames(e.h)<-newlevel

  muhat<-muhat[matchlevel]
  covmu<-covmu[matchlevel,matchlevel]
  e.h<-e.h[,matchlevel]
  muboot<-muboot[,matchlevel]


  out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
  class(out)<-'PSweight'

  out
}




