

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

  if(!is.list(xlst)){ xlst <- list(xlst)}

  valseq <- seq.int(length(xlst))

  valfn <- function(v){
    pts.extpcv <- xlst[[v]]
    model <- paste(prop,depth,"cm",sep="_")
    valtype <- pts.extpcv$valtype[1]
    ## CV statistics: all data
    RMSE = sqrt(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2, na.rm=TRUE))
    Rsq = 1-var(pts.extpcv$prop_t - pts.extpcv$pcvpred, na.rm=TRUE)/var(pts.extpcv$prop_t, na.rm=TRUE)
    Rsqpre <- 1-var(pts.extpcv$prop_t - pts.extpcv$pcvpredpre, na.rm=TRUE)/var(pts.extpcv$prop_t, na.rm=TRUE)
    Bias <- mean(pts.extpcv$prop_t - pts.extpcv$pcvpred, na.rm=TRUE)/mean(pts.extpcv$prop_t, na.rm=T)
    ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function: Using Duan's smearing estimator
    if(trans=="log10") {pts.extpcv$pcvpred_bt <- (10^(pts.extpcv$pcvpred) - 0.1)*(mean(10^(pts.extpcv$prop_t - pts.extpcv$trainpredsadj)))}
    if(trans=="log") {pts.extpcv$pcvpred_bt <- (exp(pts.extpcv$pcvpred) - 1)*(mean(exp(pts.extpcv$prop_t - pts.extpcv$trainpredsadj)))}
    if(trans=="sqrt") {pts.extpcv$pcvpred_bt <- ((pts.extpcv$pcvpred)^2)*(mean((pts.extpcv$prop_t - pts.extpcv$trainpredsadj)^2))}
    if(trans=="none") {pts.extpcv$pcvpred_bt <- pts.extpcv$pcvpred}
    ## Untransformed calcs
    RMSE_bt = sqrt(mean((pts.extpcv$prop - pts.extpcv$pcvpred_bt)^2, na.rm=TRUE))
    Rsq_bt = 1-var(pts.extpcv$prop - pts.extpcv$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv$prop, na.rm=TRUE)
    MAE_bt <- mean(abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    MedAE_bt <- median(abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    Bias_bt <- mean(pts.extpcv$prop - pts.extpcv$pcvpred_bt, na.rm=TRUE)/mean(pts.extpcv$prop_t, na.rm=T)
    ## PCV stats for scd points
    pts.extpcv.scd <- subset(pts.extpcv, pts.extpcv$tid == "scd")
    RMSE.scd <- sqrt(mean((pts.extpcv.scd$prop_t - pts.extpcv.scd$pcvpred)^2, na.rm=TRUE))
    Rsq.scd <- 1-var(pts.extpcv.scd$prop_t - pts.extpcv.scd$pcvpred, na.rm=TRUE)/var(pts.extpcv.scd$prop_t, na.rm=TRUE)
    ## PCV stats for scd points: backtransformed
    RMSE.scd_bt <- sqrt(mean((pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt)^2, na.rm=TRUE))
    Rsq.scd_bt <- 1-var(pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv.scd$prop, na.rm=TRUE)
    MAE.scd_bt <- mean(abs(pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    MedAE.scd_bt <- median(abs(pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    Bias.scd_bt <- mean(pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt, na.rm=TRUE)/mean(pts.extpcv.scd$prop_t, na.rm=T)
    ## Number of SCD samples
    n_scd <- length(pts.extpcv.scd[,1])
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
    pts.extpcv$rel.abs.resid <- pts.extpcv$abs.resid/varrange
    RPI.cvave <- mean(pts.extpcv$RPI)
    RPI.cvmed <- median(pts.extpcv$RPI)
    rel.abs.res.ave <- mean(pts.extpcv$rel.abs.resid)
    rel.abs.res.med <- median(pts.extpcv$rel.abs.resid)
    pts.extpcv$BTbias <- pts.extpcv$prop_bt - pts.extpcv$prop
    BTbias.abs.max <- max(abs(pts.extpcv$BTbias))
    BTbias.ave <- mean(pts.extpcv$BTbias)
    PICP <- sum(ifelse(pts.extpcv$prop_bt <= pts.extpcv$pcvpredpre.975_bt & pts.extpcv$prop_bt >= pts.extpcv$pcvpredpre.025_bt,1,0))/length(pts.extpcv[,1])
    ## Create CV statistics table
    CVdf <- data.frame(model, valtype, RMSE, Rsq, Rsqpre, Bias, RMSE_bt, Rsq_bt, MAE_bt, MedAE_bt, Bias_bt, RMSE.scd, Rsq.scd, RMSE.scd_bt, Rsq.scd_bt, MAE.scd_bt, MedAE.scd_bt, Bias.scd_bt, n_scd,RPI.cvave,RPI.cvmed,PICP,rel.abs.res.ave,rel.abs.res.med,BTbias.abs.max,BTbias.ave)
    names(CVdf) <- c("model","valtype","RMSE","Rsq", "Rsqpre","Bias","RMSE_bt", "Rsq_bt","MAE_bt","MedAE_bt","Bias_bt", "RMSE.scd", "Rsq.scd", "RMSE.scd_bt", "Rsq.scd_bt","MAE.scd_bt","MedAE.scd_bt","Bias.scd_bt","n_scd","RPI.CVave","RPI.CVmed","PICP","rel.abs.res.ave","rel.abs.res.med","BTbias.abs.max","BTbias.ave")
    return(CVdf)
  }

  ## list apply of function
  CVdf.lst <- lapply(valseq,valfn)
  CVdfall <- plyr::rbind.fill(CVdf.lst)
  return(CVdfall)

}

