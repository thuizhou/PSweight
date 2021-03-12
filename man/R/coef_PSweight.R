#' Point estimates of PSweight
#'
#' The \code{\link{coef}} method for class "PSweight"
#'
#' @param object an object for which the extraction of model coefficients is meaningful.
#' @param ... other arguments.
#'
#' @return The output from \code{\link{coef}}
#' @export
#'
#'
coef.PSweight <- function(object, ...)
{
  mu_coef <- object$muhat
  mu_coef
}


