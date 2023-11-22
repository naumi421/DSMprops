
# Workspace setup
# Install packages if not already installed
required.packages <- c("raster","rgdal", "rasterVis","maptools","RColorBrewer","ggplot2","gridExtra","classInt","RStoolbox","hexbin", "ranger", "parallel", "doParallel" ,"dplyr","Hmisc",
                       "viridisLite","DSMprops","stringr")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
# memory.limit(100000)
rasterOptions(maxmemory = 5e+09)

## Set cpu limit
cpus <- min(124,detectCores() - 2)

## Key Folder Locations
topfolder <- "/mnt/covs/solus_preds/v2tst_gnat_pts"
resultfolder <- "/mnt/covs/solus_preds/v2tst_gnat_pts/analysis"
ptsfolder <- "/mnt/solus100/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final"

## Create lists of folders an files to cycle through
topdirs <-  list.dirs(topfolder,recursive = F, full.names = T)
topdirs
#topdirs <- c("/mnt/covs/solus_preds/v2tst_gnat_pts/Rock_tot_gRPI_250k")
## Subsetting
topdirs <- topdirs[grepl("_250k",topdirs)]
topdirs <- topdirs[!grepl("first_test",topdirs)]
topdirs <- topdirs[!grepl("analysis",topdirs)]
topdirs <- topdirs[!grepl("models_more",topdirs)]
topdirsnm <- basename(topdirs)
#topdirsnm <- list.dirs(topfolder,recursive = F, full.names = F)
## Pedon level predictions
topdirsnmdepth <- topdirsnm[grep("RestrDepth_gRPI",topdirsnm)]
topdirsdepth <- topdirs[grep("RestrDepth_gRPI",topdirs)]
topdirs <- subset(topdirs, !topdirs %in% topdirsdepth)
topdirsnm <- subset(topdirsnm, !topdirsnm %in% topdirsnmdepth)

## Modified list for iteratively creating CVlist objects
CVtopdirs <- topdirs[grepl("Gypsum",topdirs)]
# CVtopdirs <- CVtopdirs[!grepl("/CoarseSand",CVtopdirs)]
# CVtopdirs <- CVtopdirs[!grepl("/Sand",CVtopdirs)]
# CVtopdirs <- CVtopdirs[!grepl("Silt",CVtopdirs)]
# CVtopdirs <- CVtopdirs[c(1:7,9:10)]
# CVtopdirs <- topdirs

## Read in gRPI points
pts.gRPI <- readRDS(paste(ptsfolder,"/CONUS_random_gRPIsamp_covs.rds",sep=""))

## Recreating CV objects
trn.params <- list(ntrees = 100, min.node.size = 1)
img10kf <- readRDS("/mnt/solus100/2020_Pedons/img10kf.rds")
for(d in CVtopdirs){
  dirfiles <- list.files(d, full.names = T, recursive = T)
  dirfiles <- dirfiles[!grepl("first_test", dirfiles)]
  dirfiles <- dirfiles[!grepl("results_oldSCD", dirfiles)]
  for(dpt in c("_0_","_5_","_15_","_30_","_60_","_100_","_150_")){
    depthfiles <- dirfiles[grep(dpt, dirfiles,ignore.case = TRUE)]
    ptfile <- depthfiles[grepl("TrainPTS", depthfiles)]
    rpifile <- depthfiles[grepl("gRPIs",depthfiles)]
    rffile <- depthfiles[grepl("ranger",depthfiles)]
    rf <- readRDS(rffile)
    valfile <- read.delim(depthfiles[grep("sCVstats", depthfiles,ignore.case = TRUE)])[1,] ## [1,] Just grabs the first row (CV)
    ## Determine transformation (if used)
    if(valfile$Rsq != valfile$Rsq_bt) {trans <- "log"}
    if(valfile$Rsq == valfile$Rsq_bt) {trans <- "none"}
    formulaStringRF <- as.formula(paste('prop_t ~', paste(names(rf$variable.importance), collapse="+")))
    avptfile <- depthfiles[grepl("AvailPTS_",depthfiles)]
    ptfilenm <- basename(ptfile)
    pts.pcv <- readRDS(ptfile)
    pts.all <- readRDS(avptfile)
    data_grid_df <- read.delim(rpifile)
    datagrid <- data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),c("datagrid")]
    datagrid_spl <- str_split(datagrid, "_")
    srce_idx <- which(datagrid_spl[[1]] == "srce")
    geo_idx <- which(datagrid_spl[[1]] == "geo")
    geo_vec <- paste(datagrid_spl[[1]][(geo_idx + 1):(srce_idx - 1)], collapse = "_")
    geo_vec <- paste(geo_vec, "gNAT", sep = "_")
    srce_vec <- paste(datagrid_spl[[1]][(srce_idx + 1):length(datagrid_spl[[1]])], collapse = "_")
    srce_vec <- paste(srce_vec, "gNAT", sep = "_")
    grid_vec <- expand.grid(geo=geo_vec, srce=srce_vec, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
    data_grid_new <- DSMprops::gRPI_estim_ranger(x = pts.all@data, gsamp = pts.gRPI, fm = formulaStringRF, griddf = grid_vec, os="linux",
                                                    train.params = trn.params, nthreads = cpus)
    gRPIfinal <-  (data_grid_new$gRPI.ave + data_grid_new$gRPI.med) / 2
    quants_vec <- c(data_grid_new$quant_l,data_grid_new$quant_h)
    ## Normal 10-fold cross validation
    cv10f <- DSMprops::CVranger(x = pts.pcv@data, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                nfolds = 10, nthreads = cpus, os = "linux") # 5min
    ## Spatial 10-fold cross validation on 1km blocks
    s1cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                      rast = img10kf, nfolds = 10, nthreads = cpus, resol = 1, os = "linux") # 6 min
    ## Spatial 10-fold cross validation on 10km blocks
    s10cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                       rast = img10kf, nfolds = 10, nthreads = cpus, resol = 10, os = "linux") # 6min
    ## Spatial 10-fold cross validation on 50km blocks
    s50cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                       rast = img10kf, nfolds = 10, nthreads = cpus, resol = 50, os = "linux") # 6min
    ## Spatial 10-fold cross validation on 100km blocks
    s100cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                        rast = img10kf, nfolds = 10, nthreads = cpus, resol = 100, os = "linux")
    ## Create list of CV files
    cv.lst <- list(cv10f,s1cv10f,s10cv10f,s50cv10f,s100cv10f)
    ## Validation metrics for CVs at different spatial supports
    valmets_sCV <- DSMprops::valmetrics(xlst = cv.lst, trans = trans, quants = quants_vec, prop = gsub("_gRPI_250k", "",basename(d)), depth = gsub("_","",dpt))
    valmets_sCV$RPIall <- (valmets_sCV$RPI.cvave + valmets_sCV$RPI.cvmed) / 2
    idx_val <- which(abs(valmets_sCV$RPIall-gRPIfinal)==min(abs(valmets_sCV$RPIall-gRPIfinal)))
    bestval <- valmets_sCV[idx_val,c("valtype")]
    valmets_sCV$bestval <- bestval
    # ## Pick best CV resolution based on RPI match
    # resol_rpi <- as.numeric(gsub("s", "", str_split(bestval,"cv")[[1]][1]))
    ## Save results
    write.table(valmets_sCV, paste(d,"/sCVstatsV2_", gsub("_gRPI_250k", "",basename(d)),"_", gsub("_","",dpt), "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
    ## Save file
    saveRDS(cv.lst,paste0(d,"/CVlist", gsub("TrainPTS","",ptfilenm))) # takes forever...
    print(paste(ptfilenm, "was done at", Sys.time(), sep=" "))
  }
}


## Seperate loop for depth models
for(d in topdirsdepth){
  dirfiles <- list.files(d, full.names = T, recursive = T)
  dirfiles <- dirfiles[!grepl("first_test", dirfiles)]
  dirfiles <- dirfiles[!grepl("results_oldSCD", dirfiles)]
  for(dpt in c("anylithic","resdept")){
    depthfiles <- dirfiles[grep(dpt, dirfiles,ignore.case = TRUE)]
    ptfile <- depthfiles[grepl("TrainPTS", depthfiles)]
    rpifile <- depthfiles[grepl("gRPIs",depthfiles)]
    rffile <- depthfiles[grepl("ranger",depthfiles)]
    rf <- readRDS(rffile)
    valfile <- read.delim(depthfiles[grep("sCVstats", depthfiles,ignore.case = TRUE)])[1,] ## [1,] Just grabs the first row (CV)
    ## Determine transformation (if used)
    if(valfile$Rsq != valfile$Rsq_bt) {trans <- "log"}
    if(valfile$Rsq == valfile$Rsq_bt) {trans <- "none"}
    formulaStringRF <- as.formula(paste('prop_t ~', paste(names(rf$variable.importance), collapse="+")))
    avptfile <- depthfiles[grepl("AvailPTS_",depthfiles)]
    ptfilenm <- basename(ptfile)
    pts.pcv <- readRDS(ptfile)
    pts.all <- readRDS(avptfile)
    data_grid_df <- read.delim(rpifile)
    datagrid <- data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),c("datagrid")]
    datagrid_spl <- str_split(datagrid, "_")
    srce_idx <- which(datagrid_spl[[1]] == "srce")
    geo_idx <- which(datagrid_spl[[1]] == "geo")
    geo_vec <- paste(datagrid_spl[[1]][(geo_idx + 1):(srce_idx - 1)], collapse = "_")
    geo_vec <- paste(geo_vec, "gNAT", sep = "_")
    srce_vec <- paste(datagrid_spl[[1]][(srce_idx + 1):length(datagrid_spl[[1]])], collapse = "_")
    srce_vec <- paste(srce_vec, "gNAT", sep = "_")
    grid_vec <- expand.grid(geo=geo_vec, srce=srce_vec, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
    data_grid_new <- DSMprops::gRPI_estim_ranger(x = pts.all@data, gsamp = pts.gRPI, fm = formulaStringRF, griddf = grid_vec, os="linux",
                                                 train.params = trn.params, nthreads = cpus)
    gRPIfinal <-  (data_grid_new$gRPI.ave + data_grid_new$gRPI.med) / 2
    quants_vec <- c(data_grid_new$quant_l,data_grid_new$quant_h)
    ## Normal 10-fold cross validation
    cv10f <- DSMprops::CVranger(x = pts.pcv@data, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                nfolds = 10, nthreads = cpus, os = "linux") # 5min
    ## Spatial 10-fold cross validation on 1km blocks
    s1cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                      rast = img10kf, nfolds = 10, nthreads = cpus, resol = 1, os = "linux") # 6 min
    ## Spatial 10-fold cross validation on 10km blocks
    s10cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                       rast = img10kf, nfolds = 10, nthreads = cpus, resol = 10, os = "linux") # 6min
    ## Spatial 10-fold cross validation on 50km blocks
    s50cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                       rast = img10kf, nfolds = 10, nthreads = cpus, resol = 50, os = "linux") # 6min
    ## Spatial 10-fold cross validation on 100km blocks
    s100cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,quants=quants_vec,
                                        rast = img10kf, nfolds = 10, nthreads = cpus, resol = 100, os = "linux")
    ## Create list of CV files
    cv.lst <- list(cv10f,s1cv10f,s10cv10f,s50cv10f,s100cv10f)
    ## Validation metrics for CVs at different spatial supports
    valmets_sCV <- DSMprops::valmetrics(xlst = cv.lst, trans = trans, quants = quants_vec, prop = gsub("_gRPI_250k", "",basename(d)), depth = gsub("_","",dpt))
    valmets_sCV$RPIall <- (valmets_sCV$RPI.cvave + valmets_sCV$RPI.cvmed) / 2
    idx_val <- which(abs(valmets_sCV$RPIall-gRPIfinal)==min(abs(valmets_sCV$RPIall-gRPIfinal)))
    bestval <- valmets_sCV[idx_val,c("valtype")]
    valmets_sCV$bestval <- bestval
    # ## Pick best CV resolution based on RPI match
    # resol_rpi <- as.numeric(gsub("s", "", str_split(bestval,"cv")[[1]][1]))
    ## Save results
    write.table(valmets_sCV, paste(d,"/sCVstatsV2_", gsub("_gRPI_250k", "",basename(d)),"_", gsub("_","",dpt), "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
    ## Save file
    saveRDS(cv.lst,paste0(d,"/CVlist", gsub("TrainPTS","",ptfilenm))) # takes forever...
    print(paste(ptfilenm, "was done at", Sys.time(), sep=" "))
  }
}




