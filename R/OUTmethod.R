#' Fitting potential outcome regression with different methods
#'
#' The function \code{OUTmethod} is an internal function to estimate the potential outcomes given a specified model through formula.
#' It is built into function \code{PSweight}, and is used for constructing the augmented estimators.
#'
#' @param out.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the outcome model to be fitted.
#' @param y a vector of the observed outcome in the training data (\code{datain}).
#' @param out.method a character to specify the method for estimating the outcome regression model. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param family a description of the error distribution and link function to be used in the outcome model. Supported distributional families include
#' \code{"gaussian" (link = identity)}, \code{"binomial" (link = logit)} and \code{"poisson" (link = log)}. Default is \code{"gaussian"}.
#' @param datain The training data for the outcome model. In the context of \code{PSweight}, it refers to the data observed for each treatment group.
#' @param dataout The prediction data for the outcome model. In the context of \code{PSweight}, it refers to the full data.
#' @param ... further arguments passed to or from other methods.
#'
#' @details  A typical form for \code{out.formula} is \code{y ~ terms} where \code{y} is the outcome
#' variable and \code{terms} is a series of terms which specifies linear predictors (on the link function scale). \code{out.formula} by default specifies generalized
#' linear model given the gaussian error through the default arguments \code{method = "glm"} and \code{family = "gaussian"}.  It fits the logistic regression when \code{family = "binomal"},and poisson
#' regression when \code{family = "poisson"}. The argument \code{out.method} allows user to choose
#' model other than glm to fit the outcome regression models for constructing the augmented estimator. We have included \code{gbm} and \code{SuperLearner} as alternative machine learning estimators.
#' Additional argument in them can be supplied through the \code{...} argument. Please refer to the user manual of the \code{gbm} and \code{SuperLearner} packages for all the
#' allowed arguments.
#'
#' @return
#'
#' \describe{
#'
#' \item{\code{ m.est}}{a vector of predicted outcome on the \code{dataout}.}
#'
#' \item{\code{ gamma.h}}{estimated coefficient of the outcome model when \code{method = "glm"}.}
#'
#' }
#'
#' @export
#'
#' @examples
#'
#' #' the outcome model
#' out.formula <- Y~cov1+cov2+cov3+cov4+cov5+cov6
#' y <- psdata$Y
#' #train on model with treatment group 1
#' datain <- psdata[psdata$trt==1, ]
#' outfit <- OUTmethod(out.formula = out.formula, y=y, datain = datain, dataout = psdata)
#'
OUTmethod<-function(out.formula=out.formula, y=y, out.method="glm", family='gaussian', datain=datain, dataout=dataout, out.control){

################################################################################################
  gamma.h<-NULL #only useful in glm
  if(out.method=='glm'){

    if(family=='gaussian'){
      fitglm<-lm(out.formula,data=datain)
    }else if(family=='binomial'){
      fitglm<-glm(out.formula,family = binomial(link = "logit"),data=datain)
    }else if(family=='poisson'){
      fitglm<-glm(out.formula,family = poisson(),data=datain)
    }
    m.est<-predict(fitglm,type = "response",dataout)
    gamma.h<-as.numeric(coef(fitglm))

  }else if(out.method=='gbm'){
    if ("distribution" %in% names(out.control)){
      if (!out.control$distribution %in% c("bernoulli","adaboost","gaussian","poisson")){
        stop("only bernoulli, adaboost, gaussian, poisson, supported for augmentation")
      }
    }

    if ("formula" %in% names(out.control)) {
      stop("the out formula should not be specified in 'out.control'; use argument 'out.formula'")
    } else {
      out.control$formula <- out.formula
    }

    if (!("var.monotone" %in% names(out.control))) out.control$var.monotone=NULL
    if (!("weights" %in% names(out.control))) out.control$weights=NULL
    if (!("n.trees" %in% names(out.control))) out.control$n.trees=100
    if (!("interaction.depth" %in% names(out.control))) out.control$interaction.depth=1
    if (!("n.minobsinnode" %in% names(out.control))) out.control$n.minobsinnode=10
    if (!("shrinkage" %in% names(out.control))) out.control$shrinkage=0.1
    if (!("bag.fraction" %in% names(out.control))) out.control$bag.fraction=0.5
    if (!("train.fraction" %in% names(out.control))) out.control$train.fraction=1
    if (!("cv.folds" %in% names(out.control))) out.control$cv.folds=0
    if (!("class.stratify.cv" %in% names(out.control))) out.control$class.stratify.cv=NULL
    if (!("n.cores" %in% names(out.control))) out.control$n.cores=NULL
    if ("verbose" %in% names(out.control)) warning("verbose argument set to F for gbm in PSweight")
    out.control$verbose <- FALSE
    
    if(family=='gaussian'){
      out.control$distribution <- "gaussian"
    }else if(family=='binomial'){
      if (!"distribution" %in% out.control){
        out.control$distribution <- "bernoulli"
      }
    }else if(family=='poisson'){
      out.control$distribution <- "poisson"
    }

    out.control$data <- datain

    fitgbm <- do.call(gbm::gbm, out.control)
    
    m.est<-predict(fitgbm, type="response", newdata=dataout, n.trees=m$n.trees) # stop warning with 'n.trees'
  }else if(out.method=='SuperLearner'){
      
    if ("newX" %in% names(out.control)) warning("newX argument set to NULL for SuperLearner in PSweight; please use out.method argument")
    out.control$newX <- NULL
    if ("id" %in% names(out.control)) warning("id argument set to NULL for SuperLearner in PSweight")
    out.control$id <- NULL
    if ("verbose" %in% names(out.control)) warning("verbose argument set to F for SuperLearner in PSweight")
    out.control$verbose <- FALSE
    if (!("obsWeights" %in% names(out.control))) out.control$obsWeights<-NULL
    if (!("control" %in% names(out.control))) out.control$control<-list()
    if (!("cvControl" %in% names(out.control))) out.control$cvControl<-list()
    if (!("env" %in% names(out.control))) out.control$env<-parent.frame()

    if(family=='poisson'){
      warning("poisson regression not supported in SuperLearner")
      family<-'gaussian'
    }

    SL.all<-c("SL.bartMachine","SL.bayesglm","SL.biglasso","SL.caret","SL.caret.rpart","SL.cforest",
              "SL.earth","SL.extraTrees", "SL.gam","SL.gbm","SL.glm","SL.glm.interaction","SL.glmnet","SL.ipredbagg",
              "SL.kernelKnn","SL.knn","SL.ksvm","SL.lda","SL.leekasso","SL.lm", "SL.loess","SL.logreg","SL.mean","SL.nnet",
              "SL.nnls","SL.polymars","SL.qda","SL.randomForest","SL.ranger","SL.ridge","SL.rpart","SL.rpartPrune","SL.speedglm",
              "SL.speedlm","SL.step","SL.step.forward","SL.step.interaction","SL.stepAIC","SL.svm","SL.template","SL.xgboost")


    if ("SL.library" %in% names(out.control)){
      if (length(unlist(out.control$SL.library))>1) {
        out.control$SL.library<-unlist(out.control$SL.library)[1]
        warning("only one method allowed in SL.library argument of SuperLearner in PSweight; the first element in SL.library is taken")
      }
      if (!out.control$SL.library %in% SL.all) {
        stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values")
      }
    }else{ #no SL.library specified
      if(family=='binomial'){
        out.control$SL.library="SL.glm"
      }else{
        out.control$SL.library="SL.lm"
      }
    }

    covM<-model.matrix(out.formula, datain)
    #out.method is fixed to "method.NNLS"
    fitSL <- do.call(SuperLearner::SuperLearner,
                     c(list(Y=zvalue, X=data.frame(covM), family = family, method = "method.NNLS"), out.control=out.control))

    covout<-model.matrix(out.formula,dataout)

    m.est <- predict(fitSL, newdata=data.frame(covout))$pred
  }



  return(list(m.est=m.est,gamma.h=gamma.h))
}





















