

#' Title Tabulate validation statistics
#'
#' @param xlst list of data.frames or a dataframe of cross validation predictions
#' @param trans character string data transformation type: "none", "log10", "log", or "sqrt"
#' @param varrange numeric range of values representing the quantile interval used for uncertainty prediction intervals
#' @param prop character string name of property being modeled
#' @param depth integer depth in cm for model being evaluated
#'
#' @return
#' @export
#'
valmetrics <- function(xlst, trans, varrange, prop, depth){

  if(is.data.frame(xlst)){ xlst <- list(xlst)}

  valseq <- seq.int(length(xlst))

  valfn <- function(v){
    pts.extpcv <- xlst[[v]]
    model <- paste(prop,depth,"cm",sep="_")
    valtype <- pts.extpcv$valtype[1]
    varrange <- as.numeric(quantile(pts.extpcv$prop, probs=c(0.975), na.rm=T)-quantile(pts.extpcv$prop, probs=c(0.025),na.rm=T))
    qrtrange <- as.numeric(quantile(pts.extpcv$prop_t, probs=c(0.75), na.rm=T)-quantile(pts.extpcv$prop_t, probs=c(0.25),na.rm=T))
    ## CV statistics: all data
    n <- length(pts.extpcv[,1])
    RMSE <- sqrt(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2, na.rm=TRUE))
    Rsq <- 1-var(pts.extpcv$prop_t - pts.extpcv$pcvpred, na.rm=TRUE)/var(pts.extpcv$prop_t, na.rm=TRUE)
    if("foldRsq"%in%colnames(pts.extpcv)){Rsqvar <- var(unique(pts.extpcv$foldRsq), na.rm=TRUE)}else{Rsqvar <- 'NA'}
    Rsqpre <- 1-var(pts.extpcv$prop_t - pts.extpcv$pcvpredpre, na.rm=TRUE)/var(pts.extpcv$prop_t, na.rm=TRUE)
    Bias <- mean(pts.extpcv$prop_t - pts.extpcv$pcvpred, na.rm=TRUE)/mean(pts.extpcv$prop_t, na.rm=T)
    QRMSE <- qrtrange / sqrt(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2, na.rm=TRUE))
    QMedAE <- qrtrange / median(abs(pts.extpcv$prop_t - pts.extpcv$pcvpred), na.rm=TRUE)
    ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function: Using Duan's smearing estimator
    if(trans=="log10") {pts.extpcv$pcvpred_bt <- (10^(pts.extpcv$pcvpred) - 0.1)*(mean(10^(pts.extpcv$prop_t - pts.extpcv$pcvpred)))}
    if(trans=="log") {pts.extpcv$pcvpred_bt <- (exp(pts.extpcv$pcvpred) - 1)*(mean(exp(pts.extpcv$prop_t - pts.extpcv$pcvpred)))}
    if(trans=="sqrt") {pts.extpcv$pcvpred_bt <- ((pts.extpcv$pcvpred)^2)*(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2))}
    if(trans=="none") {pts.extpcv$pcvpred_bt <- pts.extpcv$pcvpred}
    ## Untransformed calcs
    RMSE_bt <- sqrt(mean((pts.extpcv$prop - pts.extpcv$pcvpred_bt)^2, na.rm=TRUE))
    Rsq_bt <- 1-var(pts.extpcv$prop - pts.extpcv$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv$prop, na.rm=TRUE)
    MAE_bt <- mean(abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    MedAE_bt <- median(abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    Bias_bt <- mean(pts.extpcv$prop - pts.extpcv$pcvpred_bt, na.rm=TRUE)/mean(pts.extpcv$prop_t, na.rm=T)
    ## RPI
    ## Back transform low PI
    # TODO Not sure if smearing estimator should be used on PIs? Probably not since the linear
    # adjustment does not apply.
    if(trans=="log10") {pts.extpcv$pcvpredpre.025_bt <- 10^(pts.extpcv$pcvpredpre.025) - 0.1}
    if(trans=="log") {pts.extpcv$pcvpredpre.025_bt <- exp(pts.extpcv$pcvpredpre.025) - 1}
    if(trans=="sqrt") {pts.extpcv$pcvpredpre.025_bt <- (pts.extpcv$pcvpredpre.025)^2}
    if(trans=="none") {pts.extpcv$pcvpredpre.025_bt <- pts.extpcv$pcvpredpre.025}
    ## Back transform Upper PI
    if(trans=="log10") {pts.extpcv$pcvpredpre.975_bt <- 10^(pts.extpcv$pcvpredpre.975) - 0.1}
    if(trans=="log") {pts.extpcv$pcvpredpre.975_bt <- exp(pts.extpcv$pcvpredpre.975) - 1}
    if(trans=="sqrt") {pts.extpcv$pcvpredpre.975_bt <- (pts.extpcv$pcvpredpre.975)^2}
    if(trans=="none") {pts.extpcv$pcvpredpre.975_bt <- pts.extpcv$pcvpredpre.975}
    ## Calculate different metrics
    pts.extpcv$abs.resid <- abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt)
    pts.extpcv$RPI <- (pts.extpcv$pcvpredpre.975_bt - pts.extpcv$pcvpredpre.025_bt)/varrange
    ## Summarize RPI and residuals
    ## Back transform original property to avoid bias in PICP
    if(trans=="log10") {pts.extpcv$prop_bt <- 10^(pts.extpcv$prop_t) - 0.1}
    if(trans=="log") {pts.extpcv$prop_bt <- exp(pts.extpcv$prop_t) - 1}
    if(trans=="sqrt") {pts.extpcv$prop_bt <- (pts.extpcv$prop_t)^2}
    if(trans=="none") {pts.extpcv$prop_bt <- pts.extpcv$prop_t}
    RPI.cvave <- mean(pts.extpcv$RPI)
    RPI.cvmed <- median(pts.extpcv$RPI)
    PICP <- sum(ifelse(pts.extpcv$prop_bt <= pts.extpcv$pcvpredpre.975_bt & pts.extpcv$prop_bt >= pts.extpcv$pcvpredpre.025_bt,1,0))/length(pts.extpcv[,1])
    ## Create CV statistics table
    CVdf <- data.frame(model, valtype, n, RMSE, Rsq, Rsqpre, Bias, QRMSE, QMedAE, RMSE_bt, Rsq_bt, MAE_bt, MedAE_bt, Bias_bt, RPI.cvave,RPI.cvmed,PICP)
    if("cvgrid" %in% colnames(pts.extpcv)){CVdf$cvgrid <- pts.extpcv$cvgrid[1]}
    return(CVdf)
  }

  ## list apply of function
  CVdf.lst <- lapply(valseq,valfn)
  CVdfall <- do.call("rbind", CVdf.lst)
  return(CVdfall)

}

