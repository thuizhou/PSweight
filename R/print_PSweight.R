#' Print the results of PSweight
#'
#' The \code{\link{print}} method for class "PSweight"
#'
#' @param x an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#'
print.PSweight <- function(x,...) {
  group<-paste(x$group,collapse = ', ')
  cat('Original group value: ',group, "\n")
  if(!is.null(x$trtgrp)){
    cat(paste('Treatment group value: ',x$trtgrp),'\n')
  }
  cat('\n')
  cat('Point estimate:','\n')
  ptest<-paste(round(x$muhat,4),collapse = ', ')
  cat(ptest, "\n")
}



