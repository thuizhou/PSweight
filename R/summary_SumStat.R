#' Summarize a SumStat object.
#'
#' \code{summary.SumStat} is used to summarize results obtained from function
#' \code{\link{SumStat}}. The output includes effective sample sizes and tables for balance statistics.
#'
#' @param object a \code{SumStat} object obtained with the \code{\link{SumStat}} function.
#' @param weighted.var logical. Indicate whether the propensity score weighted variance should be used in calculating the balance metrics. Default is \code{TRUE}.
#' @param metric a chatacter indicating the type of balance metrics. \code{"ASD"} refers to the pairwise absolute standardized difference and \code{"PSD"} refers to the population standardized difference. Default is \code{"ASD"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details  For \code{metric}, the two options \code{"ASD"} and \code{"PSD"} are defined in Li and Li (2019)
#' for the general family of balancing weights. Similar definitions are also given in McCaffrey et al. (2013)
#' for inverse probability weighting. \code{weighted.var} specifies whether weighted or unweighted variance
#' should be used in calculating ASD or PSD. An example of weighted variance with two treatment groups is given in
#' Austin and Stuart (2015). For more than two treatment groups, the maximum of ASD (across all pairs of treatments)
#' and maximum of PSD (across all treatments) are calcualted, as explained in Li and Li (2019).
#'
#' @return A list of tables containing effective sample sizes and balance statistics on covariates
#' for specified propensity score weighting schemes.
#'
#' \describe{
#' \item{\code{ effective.sample.size}}{a table of effective sample sizes. This serves as a conservative measure to
#' characterize the variance inflation or precision loss due to weighting, see Li and Li (2019).}
#'
#' \item{\code{ unweighted}}{A table summarizing mean, variance by treatment groups, and standardized mean difference.}
#'
#' \item{\code{ IPW}}{If \code{"IPW"} is specified, this is a data table summarizing mean, variance by treatment groups,
#' and standardized mean difference under inverse probability of treatment weighting.}
#'
#' \item{\code{ treated}}{If \code{"treated"} is specified, this is a data table summarizing mean, variance by treatment groups,
#' and standardized mean difference under the ATT weights.}
#'
#' \item{\code{ overlap}}{If \code{"overlap"} is specified, this is a data table summarizing mean, variance by treatment groups,
#' and standardized mean difference under the (generalized) overlap weights.}
#'
#' \item{\code{ matching}}{If \code{"matching"} is specified, this is a data table summarizing mean, variance by treatment groups,
#' and standardized mean difference under the (generalized) matching weights.}
#'
#' \item{\code{ entropy}}{If \code{"entropy"} is specified, this is a data table summarizing mean, variance by treatment groups,
#' and standardized mean difference under the (generalized) entropy weights.}
#' }
#'
#'
#' @references
#' Crump, R. K., Hotz, V. J., Imbens, G. W., Mitnik, O. A. (2009).
#' Dealing with limited overlap in estimation of average treatment effects. Biometrika, 96(1), 187-199.
#'
#' Li, L., Greene, T. (2013).
#' A weighting analogue to pair matching in propensity score analysis. The International Journal of Biostatistics, 9(2), 215-234.
#'
#' Austin, P.C. and Stuart, E.A. (2015). Moving towards best practice when using inverse probability of treatment weighting (IPTW) using the propensity score to estimate causal treatment effects in observational studies.
#' Statistics in Medicine, 34(28), 3661-3679.
#'
#'
#' Li, F., Li, F. (2019). Propensity score weighting for causal inference with multiple treatments.
#' The Annals of Applied Statistics, 13(4), 2389-2415.
#'
#' Zhou, Y., Matsouaka, R. A., Thomas, L. (2020).
#' Propensity score weighting under limited overlap and model misspecification. Statistical Methods in Medical Research (Online)
#'
#'
#' @export
#'
#' @examples
#' ## For examples, run: example(SumStat).
#' @importFrom  stats binomial coef cov formula glm lm model.matrix plogis poisson predict qnorm quantile sd
#' @importFrom  utils capture.output combn
#' @importFrom  graphics hist legend
summary.SumStat<-function(object,weighted.var=TRUE,metric="ASD",...){

  #extract object info
  wt_list<-colnames(object$ess)[-1]
  ncate<-length(unique(object$ps.weights$zindex))
  metric<-toupper(metric)

  if (metric!="ASD" && metric!="PSD"){
    cat("metric argument unrecognized, 'ASD' calculated instead.","\n")
    metric="ASD"
  }

  output<-function(target){
    if (metric=="ASD"){
        if (weighted.var==TRUE)
        { vm<-target[,(ncate*1+1):(ncate*2)]
          SMD<-apply(abs(target[,(ncate*3+1):(ncate*4)]),1,max)
        }else{
          vm<-target[,(ncate*2+1):(ncate*3)]
          SMD<-apply(abs(target[,(ncate*4+1):(ncate*5)]),1,max)
        }
    } else{
        if (weighted.var==TRUE)
        { vm<-target[,(ncate*1+1):(ncate*2)]
          SMD<-apply(abs(target[,(ncate*5+1):(ncate*6)]),1,max)
        }else{
          vm<-target[,(ncate*2+1):(ncate*3)]
          SMD<-apply(abs(target[,(ncate*6+1):(ncate*7)]),1,max)
        }
    }
    return(cbind(target[,(ncate*0+1):(ncate*2)],SMD))
  }

  output_sumsumstat<-list(effective.sample.size=object$ess,unweighted=output(object$unweighted.sumstat))

  for (i in 1:length(wt_list)){
      wt<-paste0(wt_list[i],".sumstat",collapse = "")
      target<-object[[wt]]
      output_sumsumstat[[wt_list[i]]]<-output(target)
  }

  class(output_sumsumstat)<-'SumSumStat'
  output_sumsumstat
}


