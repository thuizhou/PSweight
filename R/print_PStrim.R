#' Print the results of PStrim
#'
#' The \code{\link{print}} method for class "PStrim"
#'
#' @param x an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{print}}
#' @export
#'

print.PStrim<- function(x,...) {
  cat("Summary of the trimming:")
  cat("\n")
  print(x$trim_sum)
  cat("\n")
}
