# The tilting function for binary treatment
tiltbin<-function(weight="overlap"){
  if(weight=="overlap"){
    return(function(x){x*(1-x)})
  }else if(weight=="IPW"){
    return(function(x){rep(1,length(x))})
  }else if(weight=="matching"){
    return(function(x){pmin(x,1-x)})
  }else if(weight=="entropy"){
    return(function(x){-(x*log(x)+(1-x)*log(1-x))})
  }else if(weight=="treated"){
    return(function(x){x})
  }
}


# The tilting function for multiple treatment
tiltmul<-function(weight="overlap"){
  if(weight=="overlap"){
    return(function(x){(1/apply(1/x,1,sum))})
  }else if(weight=="IPW"){
    return(function(x){rep(1,dim(x)[1])})
  }else if(weight=="matching"){
    return(function(x){apply(x, 1, min)})
  }else if(weight=="entropy"){
    return(function(x){ -(apply(x*log(x),1,sum))})
  }else if(weight=="treated"){
    return(function(x){x[,dim(x)[2]]})
  }
}
