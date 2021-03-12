#' Variance of PSweight
#'
#' The \code{\link{vcov}} method for class "PSweight"
#'
#' @param object an object used to select a method.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The output from \code{\link{vcov}}
#' @export
#'
#'
vcov.PSweight <- function(object, ...)
{
  mu_vcov <- object$covmu
  mu_vcov
}



