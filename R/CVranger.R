

#' Title K-fold cross validation of a ranger random forest model.
#'
#' @param x data.frame of all training instances for model
#' @param nfolds integer Number of folds to use for validation
#' @param fm Formula to use for random forest model building
#' @param os string. Operating system for implementation. Options are "windows" or "linux".
#' @param train.params List of training parameters for building random forest models.
#' @param nthreads Integer number of logical cores to use for parallelization of function.
#' @param casewts Character name of field with case weights values
#'
#' @return data.frame with original fields along with CV predictions
#' @export
#'
CVranger <- function(x, nfolds = 10, fm, os = "windows", train.params, nthreads = detectCores() - 1, casewts = "tot_wts"){

  if(!validObject(train.params)){
    message('No train.params given, using default RF params: ntree=100, min.node.size=1')
    train.params <- list(ntrees = 100, min.node.size = 1)
  }

  if(! "tot_wts" %in% colnames(x)){
    message('No case weights given, all case weights are being set to 1')
    x$casewts <- 1
  } else { x$casewts <- x[,casewts] }

  if(nthreads == 0){nthreads = 1}

  ## Linux CV function implementation in forked parallel list apply
  if(os == "linux"){
    pts.extcvm <- x
    set.seed(420)
    pts.extcvm$folds <- sample.int(nfolds,size =length(pts.extcvm[,1]),replace=T)
    fnthreads <- ifelse(nthreads > nfolds, floor((nthreads - nfolds)/nfolds), 1)
    lappthreads <- ifelse(nthreads >= nfolds, nfolds, nthreads)
    CV_factorRF <- function(g){#,pts.extcvm, formulaStringCVm){
      traindf <- subset(pts.extcvm, pts.extcvm$folds != g)
      testdf <- subset(pts.extcvm, pts.extcvm$folds == g)
      set.seed(420)
      rf.pcv <- ranger(formulaStringRF, data=traindf, num.trees = train.params$ntrees, quantreg = T, num.threads = fnthreads,
                       min.node.size = train.params$min.node.size,case.weights = traindf$casewts)
      traindf$pcvpredpre <- predict(rf.pcv,data=traindf, num.threads = fnthreads)$predictions
      testdf$pcvpredpre <- predict(rf.pcv, data=testdf, num.threads = fnthreads)$predictions
      testdf$pcvpredpre.025 <- predict(rf.pcv, data=testdf, type = "quantiles", quantiles = c(0.025), num.threads = fnthreads)$predictions
      testdf$pcvpredpre.975 <- predict(rf.pcv, data=testdf, type = "quantiles", quantiles = c(0.975), num.threads = fnthreads)$predictions
      attach(traindf)
      lm.pcv <- lm(prop_t~pcvpredpre)
      detach(traindf)
      testdf$pcvpred <- predict(lm.pcv, newdata=testdf)
      testdf$foldRsq <- 1-var(testdf$prop_t - testdf$pcvpred, na.rm=TRUE)/var(testdf$prop_t, na.rm=TRUE)
      return(testdf)
    }
    ## Linux parallel list apply
    cpus <- lappthreads
    cl <- makeCluster(cpus, type="FORK")
    registerDoParallel(cl)
    pts.extpcv.lst <- parLapply(cl,1:nfolds,try(CV_factorRF))
    stopCluster(cl)
    pts.extpcv <- plyr::rbind.fill(pts.extpcv.lst)
    pts.extpcv$pcvpred <- as.numeric(pts.extpcv$pcvpred)
    pts.extpcv$valtype <- paste("cv",nfolds,"f",sep="")
    return(pts.extpcv)
  }

  ## Windows (or linux) list apply implementation with ranger steps parallelized
  if(os == "windows"){
    pts.extcvm <- x
    set.seed(420)
    pts.extcvm$folds <- sample.int(nfolds,size =length(pts.extcvm[,1]),replace=T)
    CV_factorRF <- function(g){
      traindf <- subset(pts.extcvm, pts.extcvm$folds != g)
      testdf <- subset(pts.extcvm, pts.extcvm$folds == g)
      set.seed(420)
      rf.pcv <- ranger(formulaStringRF, data=traindf, num.trees = train.params$ntrees, quantreg = T, num.threads = nthreads,
                       min.node.size = train.params$min.node.size,case.weights = traindf$casewts)
      traindf$pcvpredpre <- predict(rf.pcv,data=traindf, num.threads = nthreads)$predictions
      testdf$pcvpredpre <- predict(rf.pcv, data=testdf, num.threads = nthreads)$predictions
      testdf$pcvpredpre.025 <- predict(rf.pcv, data=testdf, type = "quantiles", quantiles = c(0.025), num.threads = nthreads)$predictions
      testdf$pcvpredpre.975 <- predict(rf.pcv, data=testdf, type = "quantiles", quantiles = c(0.975), num.threads = nthreads)$predictions
      attach(traindf)
      lm.pcv <- lm(prop_t~pcvpredpre)
      detach(traindf)
      testdf$pcvpred <- predict(lm.pcv, newdata=testdf)
      testdf$foldRsq <- 1-var(testdf$prop_t - testdf$pcvpred, na.rm=TRUE)/var(testdf$prop_t, na.rm=TRUE)
      return(testdf)
    }
    ## list apply of function
    pts.extpcv.lst <- lapply(1:nfolds,try(CV_factorRF))
    pts.extpcv <- plyr::rbind.fill(pts.extpcv.lst)
    pts.extpcv$pcvpred <- as.numeric(pts.extpcv$pcvpred)
    pts.extpcv$valtype <- paste("cv",nfolds,"f",sep="")
    return(pts.extpcv)
  }

}
