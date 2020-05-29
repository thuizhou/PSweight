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
  bt<-is.null(x$contrastboot)

  prtgroup<-paste(x$group,collapse = ', ')
  cat('Original group value: ',prtgroup, "\n")
  if(!is.null(x$trtgrp)){
    cat(paste('Treatment group value: ',trtgrp),'\n')
  }
  cat('\n')
  cat('contrast','\n')
  print(contrast)
  cat('\n')
  if(!bt){
    cat('Use bootstrap sample for inference','\n')
    cat('\n')
  }
  print(estimates)
}

