
# Workspace setup
# Install packages if not already installed
required.packages <- c("raster","rgdal", "rasterVis","maptools","RColorBrewer","ggplot2","gridExtra","classInt","RStoolbox","hexbin", "ranger", "parallel", "doParallel" ,"dplyr","Hmisc","viridisLite")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
# memory.limit(100000)
rasterOptions(maxmemory = 5e+09)

## Key Folder Locations
topfolder <- "/media/sped/Hyb100m_gdrv/gRPI_preds"
repofolder <- "/home/tnaum/data/repos/DSMprops"
resultfolder <- "/media/sped/Hyb100m_gdrv/results"

## Create lists of folders an files to cycle through
topdirs <-  list.dirs(topfolder,recursive = F, full.names = T)
topdirsnm <- list.dirs(topfolder,recursive = F, full.names = F)
topdirsnmdepth <- topdirsnm[grep("RestrDepth_gRPI",topdirsnm)]
topdirsdepth <- topdirs[grep("RestrDepth_gRPI",topdirs)]
topdirs <- subset(topdirs, !topdirs %in% topdirsdepth)
topdirsnm <- subset(topdirsnm, !topdirsnm %in% topdirsnmdepth)

flist <- 1:length(topdirs)

cvsum_fn <- function(f) {
#for(f in 1:length(topdirs)){
  dirfiles <- list.files(topdirs[f], full.names = T, recursive = T)
  prop_fldr <- topdirsnm[f]
  newfile <- data.frame(model = "", prop="", depth="", valtype = "", data="", n = "", Rsq="", cvRPIave="",
                        cvRPImed="",cvRPIall="",PICP="",RMSE="", MAE="", MedAE="", Bias="",bestval="", gRPIave="", gRPImed="", gRPIall="")
  newfile <- newfile[-1,] # empty df
  for(d in c("_0_","_5_","_15_","_30_","_60_","_100_","_150_")){
    depthfiles <- dirfiles[grep(d, dirfiles,ignore.case = TRUE)]
    valfile <- read.delim(depthfiles[grep("sCVstats", depthfiles,ignore.case = TRUE)])[1,] ## [1,] Just grabs the first row (CV)
    # newfile <- data.frame(model = valfile$model, valtype = valfile$valtype, n = valfile$n, Rsq=valfile$Rsqpre, cvRPIave=valfile$RPI.cvave,
    #                       cvRPImed=valfile$RPI.cvmed,cvRPIall=valfile$RPIall,PICP=valfile$PICP,bestval=valfile$bestval)
    cvfile <- readRDS(depthfiles[grep("CVlist", depthfiles,ignore.case = TRUE)])
    ## Pull out CV prediction matrix to recompute metrics without RF bias adjustments for the normal 10-fold cross-val
    cv10f <- cvfile[[1]]
    ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function:
    if(valfile$Rsq != valfile$Rsq_bt) {cv10f$pcvpred_bt <- (exp(cv10f$pcvpredpre) - 1)}
    if(valfile$Rsq == valfile$Rsq_bt) {cv10f$pcvpred_bt <- cv10f$pcvpredpre}
    ## Untransformed calcs
    RMSE <- sqrt(mean((cv10f$prop - cv10f$pcvpred_bt)^2, na.rm=TRUE))
    MAE <- mean(abs(cv10f$prop - cv10f$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
    MedAE <- median(abs(cv10f$prop - cv10f$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
    Bias <- mean(cv10f$prop - cv10f$pcvpred_bt, na.rm=TRUE)/mean(cv10f$prop, na.rm=T)
    ## Now grap gRPI values
    grpifile <- read.delim(depthfiles[grep("gRPIs_", depthfiles,ignore.case = TRUE)])
    grpifile <- grpifile[grpifile$rank_final==1,] # Just grab row of top ranked model
    gRPIave <- grpifile$gRPI.ave
    gRPImed <- grpifile$gRPI.med
    gRPIall <- grpifile$gRPIfinal
    newrow <- data.frame(model = valfile$model, prop=gsub("_gRPI","",prop_fldr), depth=gsub("_","",d), valtype = valfile$valtype, data=grpifile$datagrid, n = valfile$n, Rsq=valfile$Rsqpre, cvRPIave=valfile$RPI.cvave,cvRPImed=valfile$RPI.cvmed,cvRPIall=valfile$RPIall,
                         PICP=valfile$PICP, RMSE=RMSE, MAE=MAE, MedAE=MedAE, Bias=Bias, bestval=valfile$bestval, gRPIave=gRPIave, gRPImed=gRPImed, gRPIall=gRPIall)
    newfile <- rbind(newfile, newrow)
  }
  return(newfile)
}

## Linux parallel list apply
cpus <- detectCores() - 2
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
  valfile <- read.delim(dirfiles[grep("sCVstats", dirfiles,ignore.case = TRUE)])[1,] ## [1,] Just grabs the first row (CV)
  cvfile <- readRDS(dirfiles[grep("CVlist", dirfiles,ignore.case = TRUE)])
  ## Pull out CV prediction matrix to recompute metrics without RF bias adjustments for the normal 10-fold cross-val
  cv10f <- cvfile[[1]]
  ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function:
  if(valfile$Rsq != valfile$Rsq_bt) {cv10f$pcvpred_bt <- (exp(cv10f$pcvpredpre) - 1)}
  if(valfile$Rsq == valfile$Rsq_bt) {cv10f$pcvpred_bt <- cv10f$pcvpredpre}
  ## Untransformed calcs
  RMSE <- sqrt(mean((cv10f$prop - cv10f$pcvpred_bt)^2, na.rm=TRUE))
  MAE <- mean(abs(cv10f$prop - cv10f$pcvpred_bt), na.rm=TRUE) # Mean Absolute Accuracy
  MedAE <- median(abs(cv10f$prop - cv10f$pcvpred_bt), na.rm=TRUE) # Median Absolute Accuracy
  Bias <- mean(cv10f$prop - cv10f$pcvpred_bt, na.rm=TRUE)/mean(cv10f$prop, na.rm=T)
  ## Now grap gRPI values
  grpifile <- read.delim(dirfiles[grep("gRPIs_", dirfiles,ignore.case = TRUE)])
  grpifile <- grpifile[grpifile$rank_final==1,] # Just grab row of top ranked model
  gRPIave <- grpifile$gRPI.ave
  gRPImed <- grpifile$gRPI.med
  gRPIall <- grpifile$gRPIfinal
  newrow <- data.frame(model = valfile$model, prop=t, depth=NA, valtype = valfile$valtype, data=grpifile$datagrid, n = valfile$n, Rsq=valfile$Rsqpre, cvRPIave=valfile$RPI.cvave,cvRPImed=valfile$RPI.cvmed,cvRPIall=valfile$RPIall,
                       PICP=valfile$PICP, RMSE=RMSE, MAE=MAE, MedAE=MedAE, Bias=Bias, bestval=valfile$bestval, gRPIave=gRPIave, gRPImed=gRPImed, gRPIall=gRPIall)
  cv_sum_tab <- rbind(cv_sum_tab, newrow)
}

## Now Save Table
cv_sum_tab[cv_sum_tab$prop == "Gypsum",]$data <- "geo_gps_gps2_gps3_srce_scd_direct_home" ## Correcting the decision to use this tier for all gypsum models
write.csv(cv_sum_tab,paste0(resultfolder,"/model_cv10f_summary_",gsub("-","",Sys.Date()),".csv"))


