#' Fitting propensity score models with different methods
#'
#' The function \code{PSmethod} is an internal function to estimate the propensity scores given a specified model through formula.
#' It is built into function \code{Sumstat}, \code{PStrim} and \code{PSweight}.
#'
#' @param ps.formula an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the propensity score model to be fitted. Additional details of model specification
#' are given under "Details". This argument is optional if \code{ps.estimate} is not \code{NULL}.
#' @param method a character to specify the method for estimating propensity scores. \code{"glm"} is default, and \code{"gbm"} and \code{"SuperLearner"} are also allowed.
#' @param data an optional data frame containing the variables in the propensity score model.
#' @param ncate a numeric to specify the number of treatment groups present in the given data.
#' @param ... further arguments passed to or from other methods.
#'
#' @details  A typical form for \code{ps.formula} is \code{treatment ~ terms} where \code{treatment} is the treatment
#' variable and \code{terms} is a series of terms which specifies a linear predictor. \code{ps.formula} by default specifies generalized
#' linear models given the default argument \code{method = "glm"}.  It fits the logistic regression when \code{ncate = 2},and multinomial
#' logistic regression when \code{ncate > 2}. The argument \code{method} allows user to choose
#' model other than glm to fit the propensity score models. We have included \code{gbm} and \code{SuperLearner} as two alternative machine learning methods.
#' Additional arguments of the machine learning estimators can be supplied through the \code{...} argument. Note that SuperLearner does not handle multiple groups and the current version of multinomial
#' logistic regression is not supported by gbm. We suggest user to use them with extra caution. Please refer to the user manual of the \code{gbm} and \code{SuperLearner} packages for all the
#' allowed arguments.
#'
#'
#' @return
#'
#' \describe{
#'
#' \item{\code{ e.h}}{a data frame of estimated propensity scores.}
#'
#' \item{\code{ beta.h}}{estimated coefficient of the propensity model when \code{method = "glm"}.}
#'
#' }
#'
#' @export
#'
#' @examples
#' # the propensity model
#' ps.formula <- trt~cov1+cov2+cov3+cov4+cov5+cov6
#' psfit <- PSmethod(ps.formula = ps.formula,data = psdata,ncate=3)
#'
#'
PSmethod<-function(ps.formula=ps.formula, method="glm", data=data, ncate=ncate, ps.control){

  zname<-all.vars(ps.formula)[1]

  facz<-as.factor(data[,zname])
  #creat a dictionary for the original and recoded values in Z
  dic<-levels(facz)

  ############## logistic #############################################################
  beta.h<-NULL #only return coefficient when gbm
  if (method=="glm"){
    if (exists("distribution")){
      warning("distribution argument only necessary for gbm method; set to bernoulli in glm method with logit link function by default")
    }

    if(ncate==2){
      #change z to 0/1 if
      dataps<-data

      dataps[,zname]<- as.numeric(facz)-1

      fitglm <- glm(formula = ps.formula, data=dataps,family = binomial(link = "logit"))
      e.h <- fitglm$fitted.values
      e.h <- cbind(1-e.h,e.h)
    }else{
      fitglm <- nnet::multinom(formula = ps.formula, data=data, maxit = 500, Hess = TRUE, trace = FALSE)
      e.h <- fitglm$fitted.values
    }
    beta.h<-as.numeric(t(coef(fitglm)))
  }
 
  ############## gbm ###################################################################
  if (method=="gbm"){
    if ("distribution" %in% names(ps.control)){
      if (!ps.control$distribution %in% c("bernoulli","adaboost","multinomial")){
        stop("only bernoulli, adaboost, or multinomial distributions in 'gbm' are supported in propensity score models of PSweight")
      }
    }

    if (!("var.monotone" %in% names(ps.control))) ps.control$var.monotone=NULL
    if (!("weights" %in% names(ps.control))) ps.control$weights=NULL
    if (!("n.trees" %in% names(ps.control))) ps.control$n.trees=100
    if (!("interaction.depth" %in% names(ps.control))) ps.control$interaction.depth=1
    if (!("n.minobsinnode" %in% names(ps.control))) ps.control$n.minobsinnode=10
    if (!("shrinkage" %in% names(ps.control))) ps.control$shrinkage=0.1
    if (!("bag.fraction" %in% names(ps.control))) ps.control$bag.fraction=0.5
    if (!("train.fraction" %in% names(ps.control))) ps.control$train.fraction=1
    if (!("cv.folds" %in% names(ps.control))) ps.control$cv.folds=0
    if (!("class.stratify.cv" %in% names(ps.control))) ps.control$class.stratify.cv=NULL
    if (!("n.cores" %in% names(ps.control))) ps.control$n.cores=NULL
    if ("verbose" %in% names(ps.control)) warning("verbose argument set to F for gbm in PSweight")
    ps.control$verbose <- FALSE

    if(ncate==2){
      #change z to 0/1
      dataps<-data
      dataps[,zname]<- as.numeric(facz)-1

      if ("distribution" %in% names(ps.control)) {
        if (!ps.control$distribution %in% c("adaboost","bernoulli")) {
          ps.control$distribution<-"bernoulli"
          warning("supplied unsupported distribution for binary outcome in gbm; reset to bernoulli")
        }
      }else{
        ps.control$distribution<-"bernoulli"
      }

      fitgbm <- do.call(gbm::gbm, c(list(formula = ps.formula, data=dataps), ps.control))

      e.h<-exp(fitgbm$fit)/(1+exp(fitgbm$fit))
      e.h<-cbind(1-e.h,e.h)
    }else if (ncate>2){
      if ("distribution" %in% names(ps.control)) {
        if (!ps.control$distribution=="multinomial"){
          warning("distribution for multi-category outcome reset to multinomial")
        }
      }

      ps.control$distribution<-"multinomial"

      warning("current multinomial distribution is broken in gbm; fitted results are rescaled to have rowsums of 1")

      fitgbm <- do.call(gbm::gbm, c(list(formula = ps.formula, data=data), ps.control))

      #standardize
      e.h<-predict(fitgbm, newdata = data, type = "response")[,,1]

    }
  }

############## super learner #############################################################
  if (method=="SuperLearner"){
    if ("distribution" %in% names(ps.control)){
      warning("distribution argument not supported by SuperLearner; only family argument is supported")
    }

    if ("newX" %in% names(ps.control)) warning("newX argument set to NULL for SuperLearner in PSweight; please use method argument")
    ps.control$newX <- NULL
    if ("id" %in% names(ps.control)) warning("id argument set to NULL for SuperLearner in PSweight")
    ps.control$id <- NULL
    if ("verbose" %in% names(ps.control)) warning("verbose argument set to F for SuperLearner in PSweight")
    ps.control$verbose <- FALSE
    if (!("obsWeights" %in% names(ps.control))) ps.control$obsWeights<-NULL
    if (!("control" %in% names(ps.control))) ps.control$control<-list()
    if (!("cvControl" %in% names(ps.control))) ps.control$cvControl<-list()
    if (!("env" %in% names(ps.control))) ps.control$env<-parent.frame()


    family="binomial"
    zvalue<-as.numeric(facz)-1

    if(ncate>2){
      stop("only binary outcomes are supported in SuperLearner for propensity score model in PSweight")
    }else{
      SL.all<-c("SL.bartMachine","SL.bayesglm","SL.biglasso","SL.caret","SL.caret.rpart","SL.cforest",
                "SL.earth","SL.extraTrees", "SL.gam","SL.gbm","SL.glm","SL.glm.interaction","SL.glmnet","SL.ipredbagg",
                "SL.kernelKnn","SL.knn","SL.ksvm","SL.lda","SL.leekasso","SL.lm", "SL.loess","SL.logreg","SL.mean","SL.nnet",
                "SL.nnls","SL.polymars","SL.qda","SL.randomForest","SL.ranger","SL.ridge","SL.rpart","SL.rpartPrune","SL.speedglm",
                "SL.speedlm","SL.step","SL.step.forward","SL.step.interaction","SL.stepAIC","SL.svm","SL.template","SL.xgboost")

      if ("SL.library" %in% names(ps.control)){
        if (length(unlist(ps.control$SL.library))>1) {
          ps.control$SL.library<-unlist(ps.control$SL.library)[1]
          warning("only one method allowed in SL.library argument of SuperLearner in PSweight; the first element in SL.library is taken")
        }
        if (!ps.control$SL.library %in% SL.all) {
          stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values")
        }
      }else{ #no SL.library specified
        ps.control$SL.library="SL.glm"
      }

      covM<-model.matrix(ps.formula, data)
      #method is fixed to "method.NNLS"
      fitsl <- do.call(SuperLearner::SuperLearner,
                       c(list(Y=zvalue, X=data.frame(covM), family = family, method = "method.NNLS"), ps.control))

      e.h<-cbind((1-fitsl$SL.predict),fitsl$SL.predict)
    }

  }

  #relabel the propensity name
  colnames(e.h)<-dic

  return(list(e.h=e.h,beta.h=beta.h))
}





















