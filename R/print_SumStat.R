#' Print the results of SumStat
#'
#' The \code{\link{print}} method for class "SumStat"
#'
#' @param x an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#' @return The output from \code{\link{print}}
#' @export
#'

print.SumStat<-function(x,...) {
  #weights<-paste(x$group,collapse = ', ')
  L<-length(x)
  if (names(x)[L]=="trim"){
    cat(as.numeric(rowSums(x$trim))[1],"cases trimmed, ",as.numeric(rowSums(x$trim))[2],"cases remained","\n")
    cat("\n")
    cat("trimmed result by trt group:","\n")
    cat(capture.output(x$trim),sep="\n")
    cat("\n")
  }
  weights<-names(x$ps.weights)[-(1:2)]
  if ("ATT" %in% weights){
    cat(paste('trt group for PS model is: ',x$trtgrp), "\n")
  }
  cat(paste('weights estimated for: ',paste(weights,collapse=" ")), "\n")
}
