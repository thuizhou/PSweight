#### Point and variance for binary treatment  cluster level(cl)############################
###########################################################################################
binest_cl<-function(ps.formula=NULL,zname,yname,data,trtgrp=NULL,augmentation=FALSE,bootstrap=FALSE,R=50,bs_level=NULL,out.formula=NULL,family=NULL,weight="overlap",nAGQ=1L){

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

  # obtain ps estimation
  # estimate with formula
  fit.e<-lme4::glmer(formula = ps.formula, data=data,nAGQ=nAGQ,family="binomial")
  e.h <- as.numeric(predict(fit.e,newdata=data,type="response"))
  ps.estimate<-cbind(1-e.h,e.h)


  #weight entropy needs extra clipping
  if(weight=="entropy"){
    e.h<-pmax(e.h,1e-6)
    e.h<-pmin(e.h,1-1e-6)
  }

  #tilting function
  ftilt<-tiltbin(weight = weight)

  #compute outcome regression for augmentation
  if(augmentation){

    #fit two outcome regression model for different treatment groups
    data0<-data1<-data
    data0[,zname]<-0
    data1[,zname]<-1

    if(family=="gaussian"){
      #predict outcome
      fitout<-lme4::lmer(formula = out.formula, data=data)
      m0.h<-as.numeric(predict(fitout,newdata=data0,type="response"))
      m1.h<-as.numeric(predict(fitout,newdata=data1,type="response"))
    }else{
      fitout<-lme4::glmer(formula = out.formula, data=data,nAGQ=nAGQ,family=family)
      m0.h<-as.numeric(predict(fitout,newdata=data0,type="response"))
      m1.h<-as.numeric(predict(fitout,newdata=data1,type="response"))
    }
  }



  ##No Augmentation###############################################################
  if(augmentation==FALSE){
    muhat<-ptbin(e.h,z,y,ftilt)
    ##No bootstrap###############################################################
    if(bootstrap==FALSE){
      #use conservative
      theta.h<-muhat
      covmu<-sand_bin(z=z,y=y,n=n,ftilt=ftilt,thetaest=theta.h,eest=e.h,type='ec')

      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
      muboot<-NULL
    }else{
      ##With bootstrap###############################################################
      muboot<-NULL
      #create dictionary for cluster level
      if (!is.null(bs_level)){
        cl_bs<-c(data[,bs_level])
        unique_cl<-unique(cl_bs)
        cl_dict<-list()
        for(i_c in 1:n){
          cl_dict[[cl_bs[i_c]]]<-c(cl_dict[[cl_bs[i_c]]],i_c)
        }
      }

      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }


        # estimate ps
        if (is.null(bs_level)){
          samp.b<-sample(n,n,replace = TRUE)
        }else{
          samp.l<-sample(unique_cl,replace = TRUE)
          samp.b<-c()
          for (tmp.l in samp.l){
            samp.b<-c(samp.b,cl_dict[[tmp.l]])
          }
        }

        data.b<-data[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]


        tryCatch( {
          fit.b<-lme4::glmer(formula = ps.formula, data=data.b,nAGQ=nAGQ,family="binomial")
          e.b <- as.numeric(predict(fit.b,newdata=data.b,type="response"))
          if(weight=="entropy"){
            e.b<-pmax(e.b,1e-6)
            e.b<-pmin(e.b,1-1e-6)
          }
          # point estimate
          mu.b <- ptbin(e.b,z.b,y.b,ftilt)
          muboot<-rbind(muboot,mu.b)
        },error = function(w) {
          warning("Bootstrap sample does not converge, please double check")
        })
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

      theta.h<-augest[1:6]
      covmu<-sand_bin(z,y,n,ftilt,theta.h,eest=e.h,m0est=m0.h,m1est=m1.h,type='ecac')

      muboot<-NULL
      names(muhat)<-newlevel
      colnames(covmu)<-rownames(covmu)<-newlevel
    }else{
      #use bootstrap for aumentation
      muboot<-NULL

      #create dictionary for cluster level
      if (!is.null(bs_level)){
        cl_bs<-c(data[,bs_level])
        unique_cl<-unique(cl_bs)
        cl_dict<-list()
        for(i_c in 1:n){
          cl_dict[[cl_bs[i_c]]]<-c(cl_dict[[cl_bs[i_c]]],i_c)
        }
      }

      #bootstrap runs
      for(i in 1:R){
        if(i %% 50==0){
          message("bootstrap ", i, " samples")
        }


        # estimate ps
        if (is.null(bs_level)){
          samp.b<-sample(n,n,replace = TRUE)
        }else{
          samp.l<-sample(unique_cl,replace = TRUE)
          samp.b<-c()
          for (tmp.l in samp.l){
            samp.b<-c(samp.b,cl_dict[[tmp.l]])
          }
        }

        # estimate ps
        samp.b<-sample(n,n,replace = TRUE)
        data.b<-data[samp.b,]
        y.b<-y[samp.b]
        z.b<-z[samp.b]

        #the augmentation data

        data.b<-data[samp.b,]
        data0.b<-data1.b<-data.b
        data0.b[,zname]<-0
        data1.b[,zname]<-1

        tryCatch( {
          fit.b<-lme4::glmer(formula = ps.formula, data=data.b,nAGQ=nAGQ,family="binomial")
          e.b <- as.numeric(predict(fit.b,newdata=data.b,type="response"))
          if(weight=="entropy"){
            e.b<-pmax(e.b,1e-6)
            e.b<-pmin(e.b,1-1e-6)
          }

          if(family=="gaussian"){
            #predict outcome
            fitout.b<-lme4::lmer(formula = ps.formula, data=data.b)
            m0.b<-as.numeric(predict(fitout.b,newdata=data0.b,type="response"))
            m1.b<-as.numeric(predict(fitout.b,newdata=data1.b,type="response"))
          }else{
            fitout.b<-lme4::glmer(formula = ps.formula, data=data.b,nAGQ=nAGQ,family=family)
            m0.b<-as.numeric(predict(fitout.b,newdata=data0.b,type="response"))
            m1.b<-as.numeric(predict(fitout.b,newdata=data1.b,type="response"))
          }

          mu.b <- ptbin(e.b,z.b,y.b,ftilt,m0.b,m1.b)
          muboot <- rbind(muboot,tail(mu.b,2))


        },error = function(w) {
          warning("Bootstrap sample does not converge, please double check")
        })

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




