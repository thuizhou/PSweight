# The point estimating function for binary group
ptbin<-function(psest,z,y,ftilt,m0est=NULL,m1est=NULL){

  # point estimate
  mu1est <- sum(z*y*ftilt(psest)/psest) / sum(z*ftilt(psest)/psest)
  mu0est <- sum((1-z)*y*ftilt(psest)/(1-psest)) / sum((1-z)*ftilt(psest)/(1-psest))

  if(is.null(m0est)){
    return(c(mu0est,mu1est))
  }else{
    aug0hest<-sum(ftilt(psest)*m0est)/sum(ftilt(psest))
    aug0zest<-sum((1-z)*m0est*ftilt(psest)/(1-psest))/sum((1-z)*ftilt(psest)/(1-psest))
    aug1hest<-sum(ftilt(psest)*m1est)/sum(ftilt(psest))
    aug1zest<-sum(z*m1est*ftilt(psest)/psest)/sum(z*ftilt(psest)/psest)

    mu0estaug<-mu0est+aug0hest-aug0zest
    mu1estaug<-mu1est+aug1hest-aug1zest

    return(c(mu0est,aug0hest,aug0zest,mu1est,aug1hest,aug1zest,mu0estaug,mu1estaug))
  }
}



# The point estimating function for multiple group
ptmul<-function(psest,z,y,ftilt,ncate,n,mest=NULL){

  allwt<-(1/psest)*ftilt(psest)
  wt<-rep(0,n)
  for(i in 1:ncate){
    wt[z==i]<-allwt[z==i,i]
  }

  #point estimate
  muest<-c()
  for(i in 1:ncate){
    ytmp<-y[z==i]
    wttmp<-wt[z==i]
    muest<-c(muest,sum(ytmp*wttmp)/sum(wttmp))
  }

  if(is.null(mest)){
    return(muest)
  }else{
    #calculate the augmentation term and updata mu
    augzest<-c()
    aughest<-c()
    for(i in 1:ncate){
      mtmp<-mest[,i]
      wttmp<-wt[z==i]

      augztmp<-sum(mtmp[z==i]*wttmp)/sum(wttmp)
      aughtmp<-sum(ftilt(psest)*mtmp)/sum(ftilt(psest))
      augzest<-c(augzest,augztmp)
      aughest<-c(aughest,aughtmp)
    }
    muestaug<-muest+aughest-augzest

    return(c(muest,augzest,aughest,muestaug))
  }
}

