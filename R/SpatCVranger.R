

#' Title Spatial K-fold cross validation of a ranger random forest model.
#'
#' @param sp SpatialPointsDataFrame of all training instances for model
#' @param nfolds integer Number of folds to use for validation
#' @param rast raster object to be used for builing spatial validation blocks. Must be in projected CRS in units of meters.
#' @param resol integer resolution in kilometers for spatial validation block size
#' @param fm Formula to use for random forest model building
#' @param os string. Operating system for implementation. Options are "windows" or "linux".
#' @param quants numeric vector with lower and upper quantiles for uncertainty prediction intervals
#' @param train.params List of training parameters for building random forest models.
#' @param nthreads Integer number of logical cores to use for parallelization of function.
#' @param casewts Character name of field with case weights values
#'
#' @return data.frame with original fields along with CV predictions
#' @export
#'
SpatCVranger <- function(sp, nfolds = 10, fm, rast, resol, os = "windows", quants, train.params, nthreads = detectCores() - 1, casewts = "tot_wts"){

  ## Check for train.params object
  if(!validObject(train.params)){
    message('No train.params given, using default RF params: ntree=100, min.node.size=1')
    train.params <- list(ntrees = 100, min.node.size = 1)
  }

  if(! "tot_wts" %in% colnames(sp@data)){
    message('No case weights given, all case weights are being set to 1')
    sp@data$casewts <- 1
  } else { sp@data$casewts <- sp@data[,casewts] }

  ## Check number of threads
  if(nthreads == 0){nthreads = 1}

  ## Check for matching projections
  # TODO update to proj6
  if(projection(sp) != projection(rast)){
    message('reprojecting points to match grid projections')
    sp <- spTransform(sp, CRS(projection(rast)))
  }

  ## Spatial Block creation and allocation to points
  blocksz <- resol * 1000 # Convert km to meters
  crs.rast <- CRS(projection(rast))
  img10kcv <- projectRaster(rast,res=blocksz,method='ngb',crs=crs.rast)
  set.seed(420)
  values(img10kcv) <- sample.int(nfolds, size = ncell(img10kcv), replace=T)
  names(img10kcv) <- "folds"
  pts.extcvm <- sp
  pts.extcvm <- extract(img10kcv, pts.extcvm, df=T, sp=T)
  pts.extcvm <- pts.extcvm@data

  ## Linux CV function implementation in forked parallel list apply
  if(os == "linux"){
    fnthreads <- ifelse(nthreads > nfolds, floor((nthreads - nfolds)/nfolds), 1)
    lappthreads <- ifelse(nthreads >= nfolds, nfolds, nthreads)
    CV_factorRF <- function(g){#,pts.extcvm, formulaStringCVm){
      traindf <- subset(pts.extcvm, pts.extcvm$folds != g)
      testdf <- subset(pts.extcvm, pts.extcvm$folds == g)
      set.seed(420)
      rf.pcv <- ranger(fm, data=traindf, num.trees = train.params$ntrees, quantreg = T, num.threads = fnthreads,
                       min.node.size = train.params$min.node.size,case.weights = traindf$casewts)
      traindf$pcvpredpre <- predict(rf.pcv,data=traindf, num.threads = fnthreads)$predictions
      testdf$pcvpredpre <- predict(rf.pcv, data=testdf, num.threads = fnthreads)$predictions
      testdf$pcvpredpre.025 <- predict(rf.pcv, data=testdf, type = "quantiles", quantiles = c(quants[1]), num.threads = fnthreads)$predictions
      testdf$pcvpredpre.975 <- predict(rf.pcv, data=testdf, type = "quantiles", quantiles = c(quants[2]), num.threads = fnthreads)$predictions
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
    pts.extpcv$valtype <- paste("s",resol,"cv",nfolds,"f",sep="")
    return(pts.extpcv)
  }

  ## Windows (or linux) list apply implementation with ranger steps parallelized
  if(os == "windows"){
    CV_factorRF <- function(g){
      traindf <- subset(pts.extcvm, pts.extcvm$folds != g)
      testdf <- subset(pts.extcvm, pts.extcvm$folds == g)
      set.seed(420)
      rf.pcv <- ranger(fm, data=traindf, num.trees = train.params$ntrees, quantreg = T, num.threads = nthreads,
                       min.node.size = train.params$min.node.size,case.weights = traindf$casewts)
      traindf$pcvpredpre <- predict(rf.pcv,data=traindf, num.threads = nthreads)$predictions
      testdf$pcvpredpre <- predict(rf.pcv, data=testdf, num.threads = nthreads)$predictions
      testdf$pcvpredpre.025 <- predict(rf.pcv, data=testdf, type = "quantiles", quantiles = c(quants[1]), num.threads = nthreads)$predictions
      testdf$pcvpredpre.975 <- predict(rf.pcv, data=testdf, type = "quantiles", quantiles = c(quants[2]), num.threads = nthreads)$predictions
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
    pts.extpcv$valtype <- paste("s",resol,"cv",nfolds,"f",sep="")
    return(pts.extpcv)
  }

}
