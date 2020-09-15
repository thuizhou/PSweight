#' Print the results of Summary.SumStat
#'
#' The \code{\link{print}} method for class "SumSumStat"
#'
#' @param x an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'

print.SumSumStat<- function(x,...) {
  nobj<-length(x)
  L<-dim(x$unweighted)[2]
  ncate<-(L-1)/2
  for(i in 2:nobj)
  { cat(paste0(names(x)[i], " result"), sep="\n")
    cat(capture.output(round(x[[i]][,c(1:ncate,L),drop=F],3)),sep="\n")
    cat("\n")
  }
}
