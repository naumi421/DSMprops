

#' Global relative prediction interval estimation for given sample sets over a chosen spatial prediction domain.
#'
#' @param x data.frame of all potential training instances for model
#' @param gsamp a random sample of the prediction domain with covariates pre-extracted
#' @param fm Formula to use for random forest model building
#' @param griddf A data.frame with geographic 'geo' and source 'srce' field expanded in a factorial to select
#' sample subsets for x.
#' @param os string. Operating system for implementation. Options are "windows" or "linux".
#' @param train.params List of training parameters for building random forest models.
#' @param nthreads Integer number of logical cores to use for parallelization of function.
#' @param casewts Character name of field with case weights values
#'
#' @return data.frame with original fields along with CV predictions
#' @export
#'
gRPI_estim_ranger <- function(x, gsamp, fm, os = "windows", grid, train.params, nthreads = detectCores() - 1, casewts = "tot_wts"){

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
    set.seed(420)
    fnthreads <- ifelse(nthreads > nfolds, floor((nthreads - nfolds)/nfolds), 1)
    lappthreads <- ifelse(nthreads >= nfolds, nfolds, nthreads)
    gRPI_RF <- function(g){#,pts.extcvm, formulaStringCVm){
      levs <- data.frame(grid_vec[x,])
      colnames(levs) <- colnames(grid_vec)
      geo_levs <- str_split(levs$geo, "_")[[1]]
      srce_levs <- str_split(levs$srce, "_")[[1]]
      ptseval <- x
      ptseval <- ptseval[ptseval$mtchtype %in% srce_levs,]
      # subset by geolevels, but leave in all scd since 2000 as the baseline training/evaluation data
      ptseval <- ptseval[(ptseval$mtchtype=="scd"&ptseval$geo_wt<8)|(ptseval$geo_cls %in% geo_levs),]
      #### summarize RPI in full raster prediction using sample for speed (tested against full average)
      ## Determine 95% interquartile range for relative prediction interval
      varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.975), na.rm=T)-quantile(ptseval$prop, probs=c(0.025),na.rm=T)) ## TRANSFORM IF NEEDED!
      if(varrange_gRPI != 0) {quants <- c(0.025,0.975)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.995), na.rm=T)-quantile(ptseval$prop, probs=c(0.005),na.rm=T))
      quants <- c(0.005,0.995)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.999), na.rm=T)-quantile(ptseval$prop, probs=c(0.001),na.rm=T))
      quants <- c(0.001,0.999)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.9995), na.rm=T)-quantile(ptseval$prop, probs=c(0.0005),na.rm=T))
      quants <- c(0.0005,0.9995)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.9999), na.rm=T)-quantile(ptseval$prop, probs=c(0.0001),na.rm=T))
      quants <- c(0.0001,0.9999)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(max(ptseval$prop, na.rm=T)-min(ptseval$prop,na.rm=T))
      quants <- c(0,1)}
      #gRPI_rf <- quantregForest(x=xtrain, y=ytrain, importance=TRUE, ntree=100, keep.forest=TRUE, nthreads = 60,nodesize = 1)
      gRPI_rf <- ranger(fm, data=ptseval, num.trees = trn.params$ntrees, quantreg = T, num.threads = fnthreads,
                        min.node.size = trn.params$min.node.size, importance = "impurity")
      ## Predict onto random gRPI sample pts
      gRPI_rf_quans <- predict(gRPI_rf, data=gsamp, num.threads = fnthreads, type = "quantiles", quantiles = quants)
      lowpred <- gRPI_rf_quans$predictions[,1]
      highpred <- gRPI_rf_quans$predictions[,2]
      ## Back transformation if necessary
      if(trans=="log10") {lowpred <- (10^(lowpred) - 0.1)
      highpred <- (10^(highpred) - 0.1)
      }
      if(trans=="log") {lowpred <- (exp(lowpred) - 1)
      highpred <- (exp(highpred) - 1)
      }
      if(trans=="sqrt") {lowpred <- ((lowpred)^2)
      highpred <- ((highpred)^2)
      }
      RPI <- (highpred - lowpred) / varrange_gRPI
      gRPI.ave <- mean(RPI)
      gRPI.med <- median(RPI)
      gRPI.n <- length(RPI)
      RPIg_df <- data.frame(gRPI.ave,gRPI.med,gRPI.n)
      RPIg_df$quant_l <- quants[1]
      RPIg_df$quant_h <- quants[2]
      RPIg_df$datagrid <- paste(colnames(levs),levs[1,],collapse="_",sep="_")
      RPIg_df$Rsq <- gRPI_rf$r.squared
      return(RPIg_df)
    }

    ## Linux parallel list apply
    cpus <- lappthreads
    cl <- makeCluster(cpus, type="FORK")
    registerDoParallel(cl)
    pts.extpcv.lst <- parLapply(cl,1:nrow(griddf),try(gRPI_RF))
    stopCluster(cl)
    pts.extpcv <- plyr::rbind.fill(pts.extpcv.lst)
    return(pts.extpcv)
  }

  ## Windows (or linux) list apply implementation with ranger steps parallelized
  if(os == "windows"){
    set.seed(420)
    gRPI_RF <- function(g){#,pts.extcvm, formulaStringCVm){
      levs <- data.frame(grid_vec[x,])
      colnames(levs) <- colnames(grid_vec)
      geo_levs <- str_split(levs$geo, "_")[[1]]
      srce_levs <- str_split(levs$srce, "_")[[1]]
      ptseval <- x
      ptseval <- ptseval[ptseval$mtchtype %in% srce_levs,]
      # subset by geolevels, but leave in all scd since 2000 as the baseline training/evaluation data
      ptseval <- ptseval[(ptseval$mtchtype=="scd"&ptseval$geo_wt<8)|(ptseval$geo_cls %in% geo_levs),]
      #### summarize RPI in full raster prediction using sample for speed (tested against full average)
      ## Determine 95% interquartile range for relative prediction interval
      varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.975), na.rm=T)-quantile(ptseval$prop, probs=c(0.025),na.rm=T)) ## TRANSFORM IF NEEDED!
      if(varrange_gRPI != 0) {quants <- c(0.025,0.975)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.995), na.rm=T)-quantile(ptseval$prop, probs=c(0.005),na.rm=T))
      quants <- c(0.005,0.995)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.999), na.rm=T)-quantile(ptseval$prop, probs=c(0.001),na.rm=T))
      quants <- c(0.001,0.999)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.9995), na.rm=T)-quantile(ptseval$prop, probs=c(0.0005),na.rm=T))
      quants <- c(0.0005,0.9995)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(quantile(ptseval$prop, probs=c(0.9999), na.rm=T)-quantile(ptseval$prop, probs=c(0.0001),na.rm=T))
      quants <- c(0.0001,0.9999)}
      if(varrange_gRPI == 0) {varrange_gRPI <- as.numeric(max(ptseval$prop, na.rm=T)-min(ptseval$prop,na.rm=T))
      quants <- c(0,1)}
      #gRPI_rf <- quantregForest(x=xtrain, y=ytrain, importance=TRUE, ntree=100, keep.forest=TRUE, nthreads = 60,nodesize = 1)
      gRPI_rf <- ranger(fm, data=ptseval, num.trees = trn.params$ntrees, quantreg = T, num.threads = fnthreads,
                        min.node.size = trn.params$min.node.size, importance = "impurity")
      ## Predict onto random gRPI sample pts
      gRPI_rf_quans <- predict(gRPI_rf, data=gsamp, num.threads = fnthreads, type = "quantiles", quantiles = quants)
      lowpred <- gRPI_rf_quans$predictions[,1]
      highpred <- gRPI_rf_quans$predictions[,2]
      ## Back transformation if necessary
      if(trans=="log10") {lowpred <- (10^(lowpred) - 0.1)
      highpred <- (10^(highpred) - 0.1)
      }
      if(trans=="log") {lowpred <- (exp(lowpred) - 1)
      highpred <- (exp(highpred) - 1)
      }
      if(trans=="sqrt") {lowpred <- ((lowpred)^2)
      highpred <- ((highpred)^2)
      }
      RPI <- (highpred - lowpred) / varrange_gRPI
      gRPI.ave <- mean(RPI)
      gRPI.med <- median(RPI)
      gRPI.n <- length(RPI)
      RPIg_df <- data.frame(gRPI.ave,gRPI.med,gRPI.n)
      RPIg_df$quant_l <- quants[1]
      RPIg_df$quant_h <- quants[2]
      RPIg_df$datagrid <- paste(colnames(levs),levs[1,],collapse="_",sep="_")
      RPIg_df$Rsq <- gRPI_rf$r.squared
      return(RPIg_df)
    }
    ## list apply of function
    pts.extpcv.lst <- lapply(1:nrow(griddf),try(gRPI_RF))
    pts.extpcv <- plyr::rbind.fill(pts.extpcv.lst)
    return(pts.extpcv)
  }
}
