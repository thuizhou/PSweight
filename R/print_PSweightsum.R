#' Print the results of Summary.PSweight
#'
#' The \code{\link{print}} method for class "PSweightsum"
#'
#' @param x an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return The output from \code{\link{print}}
#' @export
#'
#'
print.PSweightsum <- function(x,...) {
  group<-x$group
  trtgrp<-x$trtgrp
  contrast<-x$contrast
  estimates<-x$estimates
  type<-x$type
  CI<-x$CI
  closedform<-is.null(x$bootestimates)

  cat('\n')

  if(closedform){
    cat('Closed-form inference:','\n')
  }else{
      cat('Use Bootstrap sample for inference:','\n')
    }
  cat('\n')

  if(type!="DIF"){
    cat('Inference in log scale:','\n')
  }

  prtgroup<-paste(x$group,collapse = ', ')
  cat('Original group value: ',prtgroup, "\n")
  if(!is.null(x$trtgrp)){
    cat(paste('Treatment group value: ',trtgrp),'\n')
  }
  cat('\n')
  cat('Contrast:','\n')
  print(contrast)
  cat('\n')


  if(CI){
    printCoefmat(estimates[,c(-3),drop=FALSE],has.Pvalue = TRUE,tst.ind=0)
  }else{
    printCoefmat(estimates[,c(-4,-5),drop=FALSE],has.Pvalue = TRUE,tst.ind=0)
    }

}

