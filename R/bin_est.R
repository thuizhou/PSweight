#### Point and variance for binary treatment ##############################################
###########################################################################################

binest<-function(ps.formula=NULL,ps.estimate=NULL,zname=NULL,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=50,out.formula=NULL,out.estimate=NULL,family=NULL,weight="overlap",ps.method='glm',ps.control=list(),out.method='glm',out.control=list()){

  #preprocess formula and extract y
  out.formula<-as.formula(out.formula)
  y<-unlist(data[yname])


  #set ordered group
  facz<-as.factor(unlist(data[zname]))

  #set the treatment label
  oldlevel<-levels(facz)
  newlevel<-oldlevel

  z<-as.numeric(facz)-1
  n<-length(z) #total obs

  #set group for ATT
  if(weight=="treated"){
    if(is.null(trtgrp)){
      trtgrp<-oldlevel[2]
    }else{
      newlevel<-rev(unique(c(trtgrp,oldlevel)))
      facz<-factor(unlist(data[zname]),levels = newlevel)
      z<-as.numeric(facz)-1
    }
  }

  matchlevel<-match(oldlevel,newlevel) #match level for ATT


  #reassign value as numeric
  data[zname]<-z
  data_p<-data[,colnames(data)!=yname] #data without outcome



  # obtain ps estimation
  # estimate with formula
  if(is.null(ps.estimate)){
    fit<-do.call(PSmethod,list(ps.formula = ps.formula, method=ps.method, data=data_p,ncate=2,ps.control=ps.control))
    W<- model.matrix(ps.formula,data_p)                      # design matrix

    e.h <- as.numeric(fit$e.h[,2])
    beta.h<-as.numeric(fit$beta.h)


    if(ps.method!='glm'){
      ps.estimate<-fit$e.h
    }
  }else{
    #the name for the propensity score
    if(length(ps.estimate)==n){
      e.h<-c(ps.estimate)
    }else if(!setequal(colnames(ps.estimate),newlevel)){
      ps.estimate<-as.matrix(ps.estimate)
      e.h<-c(ps.estimate[,2])
      warning("wrong column name set for ps.estimate, treatment set as: ",newlevel[1], " , ", newlevel[2])
    }else{
      ps.estimate<-ps.estimate[,match(newlevel,colnames(ps.estimate))]
      e.h<-c(ps.estimate[,2])
    }
  }

  #weight entropy needs extra clipping
  if(weight=="entropy"){
    e.h<-pmax(e.h,1e-6)
    e.h<-pmin(e.h,1-1e-6)
  }


  #compute outcome regression for augmentation
  if(augmentation){
    #no outcome estimation provided
    if(is.null(out.estimate)){
      #fit two outcome regression model for different treatment groups
      offset.e<-rep(1,n) #for poisson regression
      dataaug<-data[,colnames(data)!=zname]
      dataaug0<-dataaug[z==0,]
      dataaug1<-dataaug[z==1,]
      XY<-model.matrix(formula(out.formula),data=dataaug)

      #predict outcome
      fitout0<-do.call(OUTmethod,list(out.formula=out.formula,y=y[z==0], out.method=out.method, family=family, datain=dataaug0, dataout=dataaug,out.control=out.control))
      m0.h<-fitout0$m.est
      gamma0.h<-fitout0$gamma.h
      fitout1<-do.call(OUTmethod,list(out.formula=out.formula,y=y[z==1], out.method=out.method, family=family, datain=dataaug1, dataout=dataaug,out.control=out.control))
      m1.h<-fitout1$m.est
      gamma1.h<-fitout1$gamma.h

      if(family=='poisson'){
        offsetlog<-model.extract(model.frame(out.formula,data = dataaug),'offset')
        if(!is.null(offsetlog)){offset.e<-exp(offsetlog)}
      }
      if(out.method!='glm'){
        out.estimate<-cbind(m0.h,m1.h)
        colnames(out.estimate)<-newlevel
      }

    }else{
      #the name for the outcome regression
      if(!setequal(colnames(out.estimate),newlevel)){
        out.estimate<-as.matrix(out.estimate)
        m0.h<-out.estimate[,1]
        m1.h<-out.estimate[,2]
        warning("wrong column name set for out.estimate, treatment set as: ",newlevel[1], " , ", newlevel[2])
      }else{
        out.estimate<-out.estimate[,match(newlevel,colnames(out.estimate))]
        m0.h<-out.estimate[,1]
        m1.h<-out.estimate[,2]
      }
    }
  }


  #tilting function
  ftilt<-tiltbin(weight = weight)

  ##No Augmentation###############################################################
  if(augmentation==FALSE){
    muhat<-ptbin(e.h,z,y,ftilt)
    ##No bootstrap###############################################################
    if(bootstrap==FALSE){

      conser<-1   #choose conservative or not

      #use ps formula
      if(is.null(ps.estimate)){
        tryCatch( {
          theta.h<-c(muhat,beta.h)
          covmu<-sand_bin(z,y,n,ftilt,theta.h,W,type='e')
          conser<-0 #is pd
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })
      }

      #use conservative
      if(conser==1){
        theta.h<-muhat
        covmu<-sand_bin(z,y,n,ftilt,theta.h,eest=e.h,type='ec')
      }

      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
      muboot<-NULL
    }else{
      ##With bootstrap###############################################################

      muboot<-NULL
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }
        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data_p[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]

        fit.b<-do.call(PSmethod,list(ps.formula = ps.formula, method=ps.method, data=data.b,ncate=2,ps.control=ps.control))
        e.b <- as.numeric(fit.b$e.h[,2])

        if(weight=="entropy"){
          e.b<-pmax(e.b,1e-6)
          e.b<-pmin(e.b,1-1e-6)
        }
        # point estimate
        mu.b <- ptbin(e.b,z.b,y.b,ftilt)
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
    augest <- ptbin(e.h,z,y,ftilt,m0.h,m1.h)
    muhat <- tail(augest,2)

    ##No bootstrap###############################################################
    if(bootstrap==FALSE){

      conser<-1       #choose conservative or not

      if(is.null(ps.estimate) & is.null(out.estimate)){

        ## both with formula
        tryCatch({
          theta.h<-c(augest[1:6],beta.h,gamma0.h,gamma1.h)
          covmu<-sand_bin(z,y,n,ftilt,theta.h,W,XY,family=family,offset.e=offset.e,type='ea')
          conser<-0 #is pd
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })

      }else if(!is.null(ps.estimate) & is.null(out.estimate)){

        ## formula on outcome regression only
        tryCatch( {
          theta.h<-c(augest[1:6],gamma0.h,gamma1.h)
          covmu<-sand_bin(z,y,n,ftilt,theta.h,XY=XY,eest=e.h,family=family,offset.e=offset.e,type='eca')
          conser<-0 #is pd
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })

      }else if(is.null(ps.estimate) & !is.null(out.estimate)){

        ## formula on ps only
        #when not pd
        tryCatch( {
          theta.h<-c(augest[1:6],beta.h)
          covmu<-sand_bin(z,y,n,ftilt,theta.h,W,m0est=m0.h,m1est=m1.h,type='eac')
          conser<-0 #is pd
        },error = function(w) {
          warning("The sandwich matrix not pd, therefore not invertable, use conservative variance instead, please double check")
        })
      }



      #if not pd then use conservative
      if(conser==1){
        theta.h<-augest[1:6]
        covmu<-sand_bin(z,y,n,ftilt,theta.h,eest=e.h,m0est=m0.h,m1est=m1.h,type='ecac')
      }

      muboot<-NULL
      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
    }else{
      #use bootstrap for aumentation
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
        y.b<-y[samp.b]
        z.b<-z[samp.b]
        fit.b<-do.call(PSmethod,list(ps.formula = ps.formula, method=ps.method, data=data.b,ncate=2,ps.control=ps.control))
        e.b <- as.numeric(fit.b$e.h[,2])

        if(weight=="entropy"){
          e.b<-pmax(e.b,1e-6)
          e.b<-pmin(e.b,1-1e-6)
        }

        #calculate the augmentation term and updata tau
        dataaug.b<-dataaug[samp.b,]
        dataaug0.b<-dataaug.b[z.b==0,]
        dataaug1.b<-dataaug.b[z.b==1,]

        #predict outcome
        fitout0.b<-do.call(OUTmethod,list(out.formula=out.formula,y=y.b[z.b==0], out.method=out.method, family=family, datain=dataaug0.b, dataout=dataaug.b,out.control=out.control))
        m0.b<-fitout0.b$m.est
        fitout1.b<-do.call(OUTmethod,list(out.formula=out.formula,y=y.b[z.b==1], out.method=out.method, family=family, datain=dataaug1.b, dataout=dataaug.b,out.control=out.control))
        m1.b<-fitout1.b$m.est

        mu.b <- ptbin(e.b,z.b,y.b,ftilt,m0.b,m1.b)
        muboot <- rbind(muboot,tail(mu.b,2))
      }

      covmu<-cov(muboot)
      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
      colnames(muboot)<-newlevel
      rownames(muboot)<-NULL
    }
  }
  e.h<-cbind(1-e.h,e.h)
  colnames(e.h)<-newlevel

  #match back to original levels
  e.h<-e.h[,matchlevel]
  muhat<-muhat[matchlevel]
  covmu<-covmu[matchlevel,matchlevel]
  muboot<-muboot[,matchlevel]

  out<-list(propensity=e.h, muhat=muhat, covmu=covmu, muboot=muboot, group=c(oldlevel), trtgrp=trtgrp)
  class(out)<-'PSweight'
  out
}




