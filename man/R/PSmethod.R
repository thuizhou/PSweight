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
PSmethod<-function(ps.formula=ps.formula, method="glm", data=data, ncate=ncate,...){

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
    if (exists("distribution")){
      if (!distribution %in% c("bernoulli","adaboost","multinomial")){
        stop("only bernoulli, adaboost, or multinomial distributions in 'gbm' are supported in propensity score models of PSweight")
      }
    }

    if (exists("var.monotone")) var.monotone=var.monotone else var.monotone=NULL
    if (exists("weights")) weights=weights else weights=NULL
    if (exists("n.trees")) n.trees=n.trees else n.trees=100
    if (exists("interaction.depth")) interaction.depth=interaction.depth else interaction.depth=1
    if (exists("n.minobsinnode")) n.minobsinnode=n.minobsinnode else n.minobsinnode=10
    if (exists("shrinkage")) shrinkage=shrinkage else shrinkage=0.1
    if (exists("bag.fraction")) bag.fraction=bag.fraction else bag.fraction=0.5
    if (exists("train.fraction")) train.fraction=train.fraction else train.fraction=1
    if (exists("cv.folds")) cv.folds=cv.folds else cv.folds=0
    if (exists("class.stratify.cv ")) class.stratify.cv =class.stratify.cv  else class.stratify.cv=NULL
    if (exists("n.cores ")) n.cores =n.cores  else n.cores=NULL
    if (exists("verbose")) warning("verbose argument set to F for SuperLearner in PSweight")

    if(ncate==2){
      #change z to 0/1
      dataps<-data
      dataps[,zname]<- as.numeric(facz)-1

      if (exists("distribution")) {
        if (!distribution %in% c("adaboost","bernoulli")) {
          distribution<-"bernoulli"
          warning("supplied unsupported distribution for binary outcome in gbm; reset to bernoulli")
        }
      }else{
        distribution<-"bernoulli"
      }

      fitgbm <- gbm::gbm(formula = ps.formula, data=dataps,distribution=distribution, var.monotone=var.monotone, n.trees=n.trees,
                       interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, bag.fraction = bag.fraction,
                       train.fraction = train.fraction, cv.folds = cv.folds, keep.data = T, verbose = F,
                       class.stratify.cv = class.stratify.cv, n.cores = n.cores)

      e.h<-exp(fitgbm$fit)/(1+exp(fitgbm$fit))
      e.h<-cbind(1-e.h,e.h)
    }else if (ncate>2){
      if (exists("distribution")) {
        if (!distribution=="multinomial"){
          warning("distribution for multi-category outcome reset to multinomial")
        }
      }

      distribution<-"multinomial"

      warning("current multinomial distribution is broken in gbm; fitted results are rescaled to have rowsums of 1")

      fitgbm <- gbm::gbm(formula = ps.formula, data=data,distribution=distribution, var.monotone=var.monotone, n.trees=n.trees,
                       interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, bag.fraction = bag.fraction,
                       train.fraction = train.fraction, cv.folds = cv.folds, keep.data = T, verbose = F,
                       class.stratify.cv = class.stratify.cv, n.cores = n.cores)
      #standardize
      e.h<-predict(fitgbm, newdata = data, type = "response")[,,1]

    }
  }

############## super learner #############################################################
  if (method=="SuperLearner"){
    if (exists("distribution")){
      warning("distribution argument not supported by SuperLearner; only family argument is supported")
    }

    if (exists("newX")) warning("newX argument set to NULL for SuperLearner in PSweight; please use method argument")
    if (exists("id")) warning("id argument set to NULL for SuperLearner in PSweight")
    if (exists("verbose")) warning("verbose argument set to F for SuperLearner in PSweight")
    if (exists("obsWeights")) obsWeights<-obsWeights else obsWeights=NULL
    if (exists("control")) control<-control else control<-list()
    if (exists("cvControl")) cvControl<-cvControl else cvControl<-list()
    if (exists("env")) env<-env else env<-parent.frame()


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

      if (exists("SL.library")){
        if (length(unlist(SL.library))>1) {
          SL.library<-unlist(SL.library)[1]
          warning("only one method allowed in SL.library argument of SuperLearner in PSweight; the first element in SL.library is taken")
        }
        if (!SL.library %in% SL.all) stop("SL.library argument unrecgonized; please use listWrappers() in SuperLearner to find the list of supported values")
      }else{ #no SL.library specified
        SL.library="SL.glm"
      }

      covM<-model.matrix(ps.formula, data)
      #method is fixed to "method.NNLS"
      fitsl <-SuperLearner::SuperLearner(Y=zvalue, X=data.frame(covM), newX = NULL, family = family, SL.library=SL.library,
                                         method = "method.NNLS", id = NULL, verbose = FALSE, control = control, cvControl = cvControl,
                                         obsWeights = obsWeights, env = env)

      e.h<-cbind((1-fitsl$SL.predict),fitsl$SL.predict)
    }

  }

  #relabel the propensity name
  colnames(e.h)<-dic

  return(list(e.h=e.h,beta.h=beta.h))
}





















