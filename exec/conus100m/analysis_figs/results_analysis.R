
# Workspace setup
# Install packages if not already installed
required.packages <- c("raster","rgdal", "rasterVis","maptools","RColorBrewer","ggplot2","gridExtra","classInt","RStoolbox","hexbin", "ranger", "parallel", "doParallel" ,"dplyr","Hmisc","viridisLite","DSMprops")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
# memory.limit(100000)
rasterOptions(maxmemory = 5e+09, chunksize = 4e+08, memfrac = 0.9)

cpus <- min(124,detectCores() - 2)

## Key Folder Locations
topfolder <- "/mnt/covs/solus_preds/v2tst_gnat_pts"
resultfolder <- "/mnt/covs/solus_preds/v2tst_gnat_pts/analysis"

## Create lists of folders an files to cycle through
topdirs <-  list.dirs(topfolder,recursive = F, full.names = T)
topdirs
#topdirs <- c("/mnt/covs/solus_preds/v2tst_gnat_pts/Rock_tot_gRPI_250k")
## Subsetting
topdirs <- topdirs[grepl("_250k",topdirs)]
topdirs <- topdirs[!grepl("first_test",topdirs)]
topdirs <- topdirs[!grepl("analysis",topdirs)]
topdirs <- topdirs[!grepl("other_test",topdirs)]
topdirs <- topdirs[!grepl("AWC_",topdirs)]
topdirsnm <- basename(topdirs)
#topdirsnm <- list.dirs(topfolder,recursive = F, full.names = F)
## Pedon level predictions
topdirsnmdepth <- topdirsnm[grep("RestrDepth_gRPI",topdirsnm)]
topdirsdepth <- topdirs[grep("RestrDepth_gRPI",topdirs)]
topdirs <- subset(topdirs, !topdirs %in% topdirsdepth)
topdirsnm <- subset(topdirsnm, !topdirsnm %in% topdirsnmdepth)


## Index for list apply
flist <- 1:length(topdirs)

cvsum_fn <- function(f) {
#for(f in 1:length(topdirs)){
  dirfiles <- list.files(topdirs[f], full.names = T, recursive = T)
  prop_fldr <- topdirsnm[f]
  newfile <- data.frame(folder = "", model = "", prop="", propsym="", depth="", valtype = "", data="", n = "", nt = "", n_l = "", nt_l = "", n_n = "", nt_n ="", n_g = "",
                        min = "", max = "", mean = "", median = "", cvRPIave="",
                        cvRPImed="",cvRPIall="",PICP="",  Rsq="", RMSE="", MAE="", MedAE="", Bias="", Rsq_l="", RMSE_l="", MAE_l="", MedAE_l="", Bias_l="",
                        Rsq_n="", RMSE_n="", MAE_n="", MedAE_n="", Bias_n="", Rsq_g="", RMSE_g="", MAE_g="", MedAE_g="", Bias_g="",
                        sRsq="", sRMSE="", sMAE="", sMedAE="", sBias="", sRsq_l="", sRMSE_l="", sMAE_l="", sMedAE_l="", sBias_l="",
                        sRsq_n="", sRMSE_n="", sMAE_n="", sMedAE_n="", sBias_n="", sRsq_g="", sRMSE_g="", sMAE_g="", sMedAE_g="", sBias_g="",
                        bestval="", gRPIave="", gRPImed="", gRPIall="")
  newfile <- newfile[-1,] # empty df
  for(d in c("_0_","_5_","_15_","_30_","_60_","_100_","_150_")){
    depthfiles <- dirfiles[grep(d, dirfiles,ignore.case = TRUE)]
    depthfiles <- depthfiles[!grepl("first_test", depthfiles)]
    depthfiles <- depthfiles[!grepl("results_oldSCD", depthfiles)]
    valfile <- read.delim(depthfiles[grep("sCVstatsV2", depthfiles,ignore.case = TRUE)])[1,] ## [1,] Just grabs the first row (CV)
    cvfile <- readRDS(depthfiles[grep("CVlist", depthfiles,ignore.case = TRUE)])
    avlfile <- readRDS(depthfiles[grep("AvailPTS", depthfiles,ignore.case = TRUE)])
    propsym <- strsplit(basename(depthfiles[grep("CVlist", depthfiles,ignore.case = TRUE)]),"_")[[1]][2]
    ## Pull out CV prediction matrix to compute metrics
    cv10f <- cvfile[[1]]
    ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function:
    if(valfile$Rsq != valfile$Rsq_bt) {cv10f$pcvpred_bt <- (exp(cv10f$pcvpred) - 1)}
    if(valfile$Rsq == valfile$Rsq_bt) {cv10f$pcvpred_bt <- cv10f$pcvpred}
    #### Normal Cross validation metrics
    Rsq <- 1-var(cv10f$prop_t - cv10f$pcvpred, na.rm=TRUE)/var(cv10f$prop_t, na.rm=TRUE)
    RMSE <- sqrt(mean((cv10f$prop - cv10f$pcvpred_bt)^2, na.rm=TRUE))
    MAE <- mean(abs(cv10f$prop - cv10f$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    MedAE <- median(abs(cv10f$prop - cv10f$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    Bias <- mean(cv10f$prop - cv10f$pcvpred_bt, na.rm=TRUE)/mean(cv10f$prop, na.rm=T)
    n_a <- nrow(avlfile@data)
    min <- min(cv10f$prop, na.rm = T)
    max <- max(cv10f$prop, na.rm = T)
    mean <- mean(cv10f$prop, na.rm = T)
    median <- median(cv10f$prop, na.rm = T)
    ## Calcs just with scd lab data
    Rsq_l <- 1-var(cv10f[cv10f$tid == "scd",]$prop_t - cv10f[cv10f$tid == "scd",]$pcvpred, na.rm=TRUE)/var(cv10f[cv10f$tid == "scd",]$prop_t, na.rm=TRUE)
    RMSE_l <- sqrt(mean((cv10f[cv10f$tid == "scd",]$prop - cv10f[cv10f$tid == "scd",]$pcvpred_bt)^2, na.rm=TRUE))
    MAE_l <- mean(abs(cv10f[cv10f$tid == "scd",]$prop - cv10f[cv10f$tid == "scd",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    MedAE_l <- median(abs(cv10f[cv10f$tid == "scd",]$prop - cv10f[cv10f$tid == "scd",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    Bias_l <- mean(cv10f[cv10f$tid == "scd",]$prop - cv10f[cv10f$tid == "scd",]$pcvpred_bt, na.rm=TRUE)/mean(cv10f[cv10f$tid == "scd",]$prop, na.rm=T)
    n_l <- nrow(avlfile@data[avlfile@data$tid=="scd",])
    nt_l <- nrow(cv10f[cv10f$tid == "scd",])
    ## Calcs just with nasis pedon data
    Rsq_n <- 1-var(cv10f[cv10f$tid == "nasis",]$prop_t - cv10f[cv10f$tid == "nasis",]$pcvpred, na.rm=TRUE)/var(cv10f[cv10f$tid == "nasis",]$prop_t, na.rm=TRUE)
    RMSE_n <- sqrt(mean((cv10f[cv10f$tid == "nasis",]$prop - cv10f[cv10f$tid == "nasis",]$pcvpred_bt)^2, na.rm=TRUE))
    MAE_n <- mean(abs(cv10f[cv10f$tid == "nasis",]$prop - cv10f[cv10f$tid == "nasis",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    MedAE_n <- median(abs(cv10f[cv10f$tid == "nasis",]$prop - cv10f[cv10f$tid == "nasis",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    Bias_n <- mean(cv10f[cv10f$tid == "nasis",]$prop - cv10f[cv10f$tid == "nasis",]$pcvpred_bt, na.rm=TRUE)/mean(cv10f[cv10f$tid == "nasis",]$prop, na.rm=T)
    n_n <- nrow(avlfile@data[avlfile@data$tid=="nasis",])
    nt_n <- nrow(cv10f[cv10f$tid == "nasis",])
    ## Calcs just with gNATSGO random point data
    Rsq_g <- 1-var(cv10f[cv10f$tid == "gNAT",]$prop_t - cv10f[cv10f$tid == "gNAT",]$pcvpred, na.rm=TRUE)/var(cv10f[cv10f$tid == "gNAT",]$prop_t, na.rm=TRUE)
    RMSE_g <- sqrt(mean((cv10f[cv10f$tid == "gNAT",]$prop - cv10f[cv10f$tid == "gNAT",]$pcvpred_bt)^2, na.rm=TRUE))
    MAE_g <- mean(abs(cv10f[cv10f$tid == "gNAT",]$prop - cv10f[cv10f$tid == "gNAT",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    MedAE_g <- median(abs(cv10f[cv10f$tid == "gNAT",]$prop - cv10f[cv10f$tid == "gNAT",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    Bias_g <- mean(cv10f[cv10f$tid == "gNAT",]$prop - cv10f[cv10f$tid == "gNAT",]$pcvpred_bt, na.rm=TRUE)/mean(cv10f[cv10f$tid == "gNAT",]$prop, na.rm=T)
    n_g <- nrow(avlfile@data[avlfile@data$tid=="gNAT",])
    #### Spatial Cross validation metrics
    if(valfile$bestval == "cv10f"){sCV_idx <- 1}
    if(valfile$bestval == "s1cv10f"){sCV_idx <- 2}
    if(valfile$bestval == "s10cv10f"){sCV_idx <- 3}
    if(valfile$bestval == "s50cv10f"){sCV_idx <- 4}
    if(valfile$bestval == "s100cv10f"){sCV_idx <- 5}
    sCV <- cvfile[[sCV_idx]]
    ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function:
    if(valfile$Rsq != valfile$Rsq_bt) {sCV$pcvpred_bt <- (exp(sCV$pcvpred) - 1)}
    if(valfile$Rsq == valfile$Rsq_bt) {sCV$pcvpred_bt <- sCV$pcvpred}
    ## Spatial CV stats
    sRsq <- 1-var(sCV$prop_t - sCV$pcvpred, na.rm=TRUE)/var(sCV$prop_t, na.rm=TRUE)
    sRMSE <- sqrt(mean((sCV$prop - sCV$pcvpred_bt)^2, na.rm=TRUE))
    sMAE <- mean(abs(sCV$prop - sCV$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    sMedAE <- median(abs(sCV$prop - sCV$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    sBias <- mean(sCV$prop - sCV$pcvpred_bt, na.rm=TRUE)/mean(sCV$prop, na.rm=T)
    ## Calcs just with scd lab data
    sRsq_l <- 1-var(sCV[sCV$tid == "scd",]$prop_t - sCV[sCV$tid == "scd",]$pcvpred, na.rm=TRUE)/var(sCV[sCV$tid == "scd",]$prop_t, na.rm=TRUE)
    sRMSE_l <- sqrt(mean((sCV[sCV$tid == "scd",]$prop - sCV[sCV$tid == "scd",]$pcvpred_bt)^2, na.rm=TRUE))
    sMAE_l <- mean(abs(sCV[sCV$tid == "scd",]$prop - sCV[sCV$tid == "scd",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    sMedAE_l <- median(abs(sCV[sCV$tid == "scd",]$prop - sCV[sCV$tid == "scd",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    sBias_l <- mean(sCV[sCV$tid == "scd",]$prop - sCV[sCV$tid == "scd",]$pcvpred_bt, na.rm=TRUE)/mean(sCV[sCV$tid == "scd",]$prop, na.rm=T)
    ## Calcs just with nasis pedon data
    sRsq_n <- 1-var(sCV[sCV$tid == "nasis",]$prop_t - sCV[sCV$tid == "nasis",]$pcvpred, na.rm=TRUE)/var(sCV[sCV$tid == "nasis",]$prop_t, na.rm=TRUE)
    sRMSE_n <- sqrt(mean((sCV[sCV$tid == "nasis",]$prop - sCV[sCV$tid == "nasis",]$pcvpred_bt)^2, na.rm=TRUE))
    sMAE_n <- mean(abs(sCV[sCV$tid == "nasis",]$prop - sCV[sCV$tid == "nasis",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    sMedAE_n <- median(abs(sCV[sCV$tid == "nasis",]$prop - sCV[sCV$tid == "nasis",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    sBias_n <- mean(sCV[sCV$tid == "nasis",]$prop - sCV[sCV$tid == "nasis",]$pcvpred_bt, na.rm=TRUE)/mean(sCV[sCV$tid == "nasis",]$prop, na.rm=T)
    ## Calcs just with gNATSGO random point data
    sRsq_g <- 1-var(sCV[sCV$tid == "gNAT",]$prop_t - sCV[sCV$tid == "gNAT",]$pcvpred, na.rm=TRUE)/var(sCV[sCV$tid == "gNAT",]$prop_t, na.rm=TRUE)
    sRMSE_g <- sqrt(mean((sCV[sCV$tid == "gNAT",]$prop - sCV[sCV$tid == "gNAT",]$pcvpred_bt)^2, na.rm=TRUE))
    sMAE_g <- mean(abs(sCV[sCV$tid == "gNAT",]$prop - sCV[sCV$tid == "gNAT",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    sMedAE_g <- median(abs(sCV[sCV$tid == "gNAT",]$prop - sCV[sCV$tid == "gNAT",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    sBias_g <- mean(sCV[sCV$tid == "gNAT",]$prop - sCV[sCV$tid == "gNAT",]$pcvpred_bt, na.rm=TRUE)/mean(sCV[sCV$tid == "gNAT",]$prop, na.rm=T)
    ## Now grap gRPI values
    grpifile <- read.delim(depthfiles[grep("gRPIs_", depthfiles,ignore.case = TRUE)])
    grpifile <- grpifile[grpifile$rank_final==1,] # Just grab row of top ranked model
    gRPIave <- grpifile$gRPI.ave
    gRPImed <- grpifile$gRPI.med
    gRPIall <- grpifile$gRPIfinal
    newrow <- data.frame(folder = prop_fldr,model = valfile$model, prop=gsub("_gRPI","",prop_fldr), propsym=propsym, depth=gsub("_","",d), valtype = valfile$valtype, data=grpifile$datagrid,
                         n = n_a, n_t = valfile$n, n_l = n_l, nt_l = nt_l, n_n = n_n, nt_n = nt_n, n_g = n_g,
                         min = min, max = max, mean = mean, median = median,
                         cvRPIave=valfile$RPI.cvave,cvRPImed=valfile$RPI.cvmed,cvRPIall=valfile$RPIall, PICP=valfile$PICP,
                         Rsq = Rsq, RMSE=RMSE, MAE=MAE, MedAE=MedAE, Bias=Bias,
                         Rsq_l = Rsq_l, RMSE_l=RMSE_l, MAE_l=MAE_l, MedAE_l=MedAE_l, Bias_l=Bias_l,
                         Rsq_n=Rsq_n, RMSE_n=RMSE_n, MAE_n=MAE_n, MedAE_n=MedAE_n, Bias_n=Bias_n,
                         Rsq_g=Rsq_g, RMSE_g=RMSE_g, MAE_g=MAE_g, MedAE_g=MedAE_g, Bias_g=Bias_g,
                         sRsq=sRsq, sRMSE=sRMSE, sMAE=sMAE, sMedAE=sMedAE, sBias=sBias,
                         sRsq_l=sRsq_l, sRMSE_l=sRMSE_l, sMAE_l=sMAE_l, sMedAE_l=sMedAE_l, sBias_l=sBias_l,
                         sRsq_n=sRsq_n, sRMSE_n=sRMSE_n, sMAE_n=sMAE_n, sMedAE_n=sMedAE_n, sBias_n=sBias_n,
                         sRsq_g=sRsq_g, sRMSE_g=sRMSE_g, sMAE_g=sMAE_g, sMedAE_g=sMedAE_g, sBias_g=sBias_g,
                         bestval=valfile$bestval, gRPIave=gRPIave, gRPImed=gRPImed, gRPIall=gRPIall)
    newfile <- rbind(newfile, newrow)
  }
  return(newfile)
}

## Linux parallel list apply
cpus <- min(length(flist),detectCores() - 2)
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
cvsum_fn.lst <- parLapply(cl,flist,try(cvsum_fn)) ## Only returning 0 cm depth, need 'return' command outside of loop. Need rbind in loop
stopCluster(cl)
cv_sum_tab <- plyr::rbind.fill(cvsum_fn.lst)


## Now add depth to restriction models
depthtypes <- c("anylithic","resdept")
for(t in depthtypes){
  dirfiles <- list.files(topdirsdepth, full.names = T, recursive = T)
  dirfiles <- dirfiles[grep(t, dirfiles,ignore.case = TRUE)]
  dirfiles <- dirfiles[!grepl("results_oldSCD", dirfiles)]
  prop_fldr <- t
  valfile <- read.delim(dirfiles[grep("sCVstatsV2", dirfiles,ignore.case = TRUE)])[1,] ## [1,] Just grabs the first row (CV)
  cvfile <- readRDS(dirfiles[grep("CVlist", dirfiles,ignore.case = TRUE)])
  avlfile <- readRDS(dirfiles[grep("AvailPTS", dirfiles,ignore.case = TRUE)])
  propsym <- strsplit(basename(dirfiles[grep("CVlist", dirfiles,ignore.case = TRUE)]),"_")[[1]][2]
  ## Pull out CV prediction matrix to recompute metrics without RF bias adjustments for the normal 10-fold cross-val
  cv10f <- cvfile[[1]]
  ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function:
  if(valfile$Rsq != valfile$Rsq_bt) {cv10f$pcvpred_bt <- (exp(cv10f$pcvpred) - 1)}
  if(valfile$Rsq == valfile$Rsq_bt) {cv10f$pcvpred_bt <- cv10f$pcvpred}
  #### Normal Cross validation metrics
  Rsq <- 1-var(cv10f$prop_t - cv10f$pcvpred, na.rm=TRUE)/var(cv10f$prop_t, na.rm=TRUE)
  RMSE <- sqrt(mean((cv10f$prop - cv10f$pcvpred_bt)^2, na.rm=TRUE))
  MAE <- mean(abs(cv10f$prop - cv10f$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  MedAE <- median(abs(cv10f$prop - cv10f$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  Bias <- mean(cv10f$prop - cv10f$pcvpred_bt, na.rm=TRUE)/mean(cv10f$prop, na.rm=T)
  n_a <- nrow(avlfile@data)
  min <- min(cv10f$prop, na.rm = T)
  max <- max(cv10f$prop, na.rm = T)
  mean <- mean(cv10f$prop, na.rm = T)
  median <- median(cv10f$prop, na.rm = T)
  ## Calcs just with scd lab data
  Rsq_l <- 1-var(cv10f[cv10f$tid == "scd",]$prop_t - cv10f[cv10f$tid == "scd",]$pcvpred, na.rm=TRUE)/var(cv10f[cv10f$tid == "scd",]$prop_t, na.rm=TRUE)
  RMSE_l <- sqrt(mean((cv10f[cv10f$tid == "scd",]$prop - cv10f[cv10f$tid == "scd",]$pcvpred_bt)^2, na.rm=TRUE))
  MAE_l <- mean(abs(cv10f[cv10f$tid == "scd",]$prop - cv10f[cv10f$tid == "scd",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  MedAE_l <- median(abs(cv10f[cv10f$tid == "scd",]$prop - cv10f[cv10f$tid == "scd",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  Bias_l <- mean(cv10f[cv10f$tid == "scd",]$prop - cv10f[cv10f$tid == "scd",]$pcvpred_bt, na.rm=TRUE)/mean(cv10f[cv10f$tid == "scd",]$prop, na.rm=T)
  n_l <- nrow(avlfile@data[avlfile@data$tid=="scd",])
  nt_l <- nrow(cv10f[cv10f$tid == "scd",])
  ## Calcs just with nasis pedon data
  Rsq_n <- 1-var(cv10f[cv10f$tid == "nasis",]$prop_t - cv10f[cv10f$tid == "nasis",]$pcvpred, na.rm=TRUE)/var(cv10f[cv10f$tid == "nasis",]$prop_t, na.rm=TRUE)
  RMSE_n <- sqrt(mean((cv10f[cv10f$tid == "nasis",]$prop - cv10f[cv10f$tid == "nasis",]$pcvpred_bt)^2, na.rm=TRUE))
  MAE_n <- mean(abs(cv10f[cv10f$tid == "nasis",]$prop - cv10f[cv10f$tid == "nasis",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  MedAE_n <- median(abs(cv10f[cv10f$tid == "nasis",]$prop - cv10f[cv10f$tid == "nasis",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  Bias_n <- mean(cv10f[cv10f$tid == "nasis",]$prop - cv10f[cv10f$tid == "nasis",]$pcvpred_bt, na.rm=TRUE)/mean(cv10f[cv10f$tid == "nasis",]$prop, na.rm=T)
  n_n <- nrow(avlfile@data[avlfile@data$tid=="nasis",])
  nt_n <- nrow(cv10f[cv10f$tid == "nasis",])
  ## Calcs just with gNATSGO random point data
  Rsq_g <- 1-var(cv10f[cv10f$tid == "gNAT",]$prop_t - cv10f[cv10f$tid == "gNAT",]$pcvpred, na.rm=TRUE)/var(cv10f[cv10f$tid == "gNAT",]$prop_t, na.rm=TRUE)
  RMSE_g <- sqrt(mean((cv10f[cv10f$tid == "gNAT",]$prop - cv10f[cv10f$tid == "gNAT",]$pcvpred_bt)^2, na.rm=TRUE))
  MAE_g <- mean(abs(cv10f[cv10f$tid == "gNAT",]$prop - cv10f[cv10f$tid == "gNAT",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  MedAE_g <- median(abs(cv10f[cv10f$tid == "gNAT",]$prop - cv10f[cv10f$tid == "gNAT",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  Bias_g <- mean(cv10f[cv10f$tid == "gNAT",]$prop - cv10f[cv10f$tid == "gNAT",]$pcvpred_bt, na.rm=TRUE)/mean(cv10f[cv10f$tid == "gNAT",]$prop, na.rm=T)
  n_g <- nrow(avlfile@data[avlfile@data$tid=="gNAT",])
  #### Spatial Cross validation metrics
  if(valfile$bestval == "cv10f"){sCV_idx <- 1}
  if(valfile$bestval == "s1cv10f"){sCV_idx <- 2}
  if(valfile$bestval == "s10cv10f"){sCV_idx <- 3}
  if(valfile$bestval == "s50cv10f"){sCV_idx <- 4}
  if(valfile$bestval == "s100cv10f"){sCV_idx <- 5}
  sCV <- cvfile[[sCV_idx]]
  ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function:
  if(valfile$Rsq != valfile$Rsq_bt) {sCV$pcvpred_bt <- (exp(sCV$pcvpred) - 1)}
  if(valfile$Rsq == valfile$Rsq_bt) {sCV$pcvpred_bt <- sCV$pcvpred}
  ## Spatial CV stats
  sRsq <- 1-var(sCV$prop_t - sCV$pcvpred, na.rm=TRUE)/var(sCV$prop_t, na.rm=TRUE)
  sRMSE <- sqrt(mean((sCV$prop - sCV$pcvpred_bt)^2, na.rm=TRUE))
  sMAE <- mean(abs(sCV$prop - sCV$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  sMedAE <- median(abs(sCV$prop - sCV$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  sBias <- mean(sCV$prop - sCV$pcvpred_bt, na.rm=TRUE)/mean(sCV$prop, na.rm=T)
  ## Calcs just with scd lab data
  sRsq_l <- 1-var(sCV[sCV$tid == "scd",]$prop_t - sCV[sCV$tid == "scd",]$pcvpred, na.rm=TRUE)/var(sCV[sCV$tid == "scd",]$prop_t, na.rm=TRUE)
  sRMSE_l <- sqrt(mean((sCV[sCV$tid == "scd",]$prop - sCV[sCV$tid == "scd",]$pcvpred_bt)^2, na.rm=TRUE))
  sMAE_l <- mean(abs(sCV[sCV$tid == "scd",]$prop - sCV[sCV$tid == "scd",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  sMedAE_l <- median(abs(sCV[sCV$tid == "scd",]$prop - sCV[sCV$tid == "scd",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  sBias_l <- mean(sCV[sCV$tid == "scd",]$prop - sCV[sCV$tid == "scd",]$pcvpred_bt, na.rm=TRUE)/mean(sCV[sCV$tid == "scd",]$prop, na.rm=T)
  ## Calcs just with nasis pedon data
  sRsq_n <- 1-var(sCV[sCV$tid == "nasis",]$prop_t - sCV[sCV$tid == "nasis",]$pcvpred, na.rm=TRUE)/var(sCV[sCV$tid == "nasis",]$prop_t, na.rm=TRUE)
  sRMSE_n <- sqrt(mean((sCV[sCV$tid == "nasis",]$prop - sCV[sCV$tid == "nasis",]$pcvpred_bt)^2, na.rm=TRUE))
  sMAE_n <- mean(abs(sCV[sCV$tid == "nasis",]$prop - sCV[sCV$tid == "nasis",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  sMedAE_n <- median(abs(sCV[sCV$tid == "nasis",]$prop - sCV[sCV$tid == "nasis",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  sBias_n <- mean(sCV[sCV$tid == "nasis",]$prop - sCV[sCV$tid == "nasis",]$pcvpred_bt, na.rm=TRUE)/mean(sCV[sCV$tid == "nasis",]$prop, na.rm=T)
  ## Calcs just with gNATSGO random point data
  sRsq_g <- 1-var(sCV[sCV$tid == "gNAT",]$prop_t - sCV[sCV$tid == "gNAT",]$pcvpred, na.rm=TRUE)/var(sCV[sCV$tid == "gNAT",]$prop_t, na.rm=TRUE)
  sRMSE_g <- sqrt(mean((sCV[sCV$tid == "gNAT",]$prop - sCV[sCV$tid == "gNAT",]$pcvpred_bt)^2, na.rm=TRUE))
  sMAE_g <- mean(abs(sCV[sCV$tid == "gNAT",]$prop - sCV[sCV$tid == "gNAT",]$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  sMedAE_g <- median(abs(sCV[sCV$tid == "gNAT",]$prop - sCV[sCV$tid == "gNAT",]$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  sBias_g <- mean(sCV[sCV$tid == "gNAT",]$prop - sCV[sCV$tid == "gNAT",]$pcvpred_bt, na.rm=TRUE)/mean(sCV[sCV$tid == "gNAT",]$prop, na.rm=T)
  ## Now grap gRPI values
  grpifile <- read.delim(dirfiles[grep("gRPIs_", dirfiles,ignore.case = TRUE)])
  grpifile <- grpifile[grpifile$rank_final==1,] # Just grab row of top ranked model
  gRPIave <- grpifile$gRPI.ave
  gRPImed <- grpifile$gRPI.med
  gRPIall <- grpifile$gRPIfinal
  newrow <- data.frame(folder = prop_fldr,model = valfile$model, prop=t, propsym = propsym, depth="NA", valtype = valfile$valtype, data=grpifile$datagrid,
                       n = n_a, n_t = valfile$n, n_l = n_l, nt_l = nt_l, n_n = n_n, nt_n = nt_n, n_g = n_g,
                       min = min, max = max, mean = mean, median = median,
                       cvRPIave=valfile$RPI.cvave,cvRPImed=valfile$RPI.cvmed,cvRPIall=valfile$RPIall, PICP=valfile$PICP,
                       Rsq = Rsq, RMSE=RMSE, MAE=MAE, MedAE=MedAE, Bias=Bias,
                       Rsq_l = Rsq_l, RMSE_l=RMSE_l, MAE_l=MAE_l, MedAE_l=MedAE_l, Bias_l=Bias_l,
                       Rsq_n=Rsq_n, RMSE_n=RMSE_n, MAE_n=MAE_n, MedAE_n=MedAE_n, Bias_n=Bias_n,
                       Rsq_g=Rsq_g, RMSE_g=RMSE_g, MAE_g=MAE_g, MedAE_g=MedAE_g, Bias_g=Bias_g,
                       sRsq=sRsq, sRMSE=sRMSE, sMAE=sMAE, sMedAE=sMedAE, sBias=sBias,
                       sRsq_l=sRsq_l, sRMSE_l=sRMSE_l, sMAE_l=sMAE_l, sMedAE_l=sMedAE_l, sBias_l=sBias_l,
                       sRsq_n=sRsq_n, sRMSE_n=sRMSE_n, sMAE_n=sMAE_n, sMedAE_n=sMedAE_n, sBias_n=sBias_n,
                       sRsq_g=sRsq_g, sRMSE_g=sRMSE_g, sMAE_g=sMAE_g, sMedAE_g=sMedAE_g, sBias_g=sBias_g,
                       bestval=valfile$bestval, gRPIave=gRPIave, gRPImed=gRPImed, gRPIall=gRPIall)
  cv_sum_tab <- rbind(cv_sum_tab, newrow)
}

## Now Save Table
write.csv(cv_sum_tab,paste0(resultfolder,"/model_cv_summary_",gsub("-","",Sys.Date()),".csv"))
cv_sum_tab <- read.csv(paste0(resultfolder,"/model_cv_summary_",gsub("-","",Sys.Date()),".csv"))

## Summary figure of Rsq
cv_sum_tab$depthn <- as.numeric(cv_sum_tab$depth)
cv_sum_tab$depth <- ifelse(is.na(cv_sum_tab$depth) ,"none",cv_sum_tab$depth)
cv_sum_tab$depth <- ordered(cv_sum_tab$depth,levels=c("0","5","15","30","60","100","150","none"))
cv_sum_tab$propnm <- gsub("_250k","",cv_sum_tab$prop)
cv_sum_tab <- cv_sum_tab[cv_sum_tab$propnm != "VCoarseSand",]
cv_sum_tab[cv_sum_tab$propnm == "VCoarseSand_wtrans",c("propnm")] <- "VCoarseSand"
viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
Rsq.plt <- ggplot(cv_sum_tab, aes(factor(propnm),Rsq,fill=depth)) + #,order = -as.numeric(depthn)
  geom_boxplot(fill="white")+ylim(0.44,0.9) + ylab(expression(R^{2})) + xlab("") +
  geom_dotplot(binaxis='y', stackdir='center',stackratio=1.5, dotsize=0.5) +ggtitle("Cross Validation Coefficients of Determination")+
  theme(axis.text=element_text(size=10), legend.text=element_text(size=10), axis.title=element_text(size=11),plot.title = element_text(size=12,hjust=0.5)) +
  coord_flip()+#scale_fill_manual(values=c(gray(0.1), gray(0.2), gray(0.3),gray(0.4),gray(0.5),gray(0.6),gray(0.8)))
  scale_fill_manual(values=c(viridis(7),gray(0.1)))#
Rsq.plt
ggsave(paste0(resultfolder,'/Rsq_model_compare_plots_',gsub("-","",Sys.Date()),'.tif'), plot = Rsq.plt, device = "tiff", dpi = 600, limitsize = TRUE, width = 9, height = 6, units = 'in',compression = c("lzw"))




