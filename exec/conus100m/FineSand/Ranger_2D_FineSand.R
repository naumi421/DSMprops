######################
## Random Forest script that includes:
## Extraction of covariates to points
## Prediction interval creation
## Cross Validation
## Most steps parallelized
######################


# Workspace setup
# Install packages if not already installed
required.packages <- c("raster", "sp", "rgdal", "ranger", "snow", "snowfall", "dplyr", "ggplot2",
                       "hexbin","doParallel","aqp","Hmisc","spatstat","maptools","DSMprops","RSQLite")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package: Windows only
#memory.limit(500000)
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)

## Key Folder Locations
predfolder <- "/push/HYBconus100m/FineSand"
repofolder <- "/push/repos/DSMprops/exec/conus100m/FineSand"
covfolder <- "/push/SG100_covars"
ptsfolder <- "/push/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final"

######## Load soil profile collection ##############
pts <- readRDS(paste(ptsfolder,"/NASIS_all_component_horizon_match_SPC_ssurgo20.rds",sep=""))
pts.proj <- proj4string(pts)
shp.pts <- pts@site
coordinates(shp.pts) <- ~ x_std + y_std
shp.pts <- as(pts,'SpatialPointsDataFrame')
projection(shp.pts) <- pts.proj
shp.pts@data <- shp.pts@data[,c("peiid","mtchtype","compname")]


######### Grid Prep #################
## Make list of grids
cov.grids <- list.files(path = covfolder,pattern=".tif$",full.names = T, recursive = F)
cov.grids.names <- basename(cov.grids)
## If points need to be matched up to grids ###
projgrid <- raster(paste(covfolder,"/T07PRI5.tif",sep="")) # This raster is fully clipped
cov.proj <- projection(projgrid)
shp.pts <- spTransform(shp.pts, CRS(cov.proj)) # project to match rasters

####### Polygon boundary if needed to clip down
polybound <- readOGR("/ped/GIS_Archive/US_boundaries/states_21basic", "conus_bound")
polybound <- spTransform(polybound, cov.proj)
shp.pts <- shp.pts[polybound,]

## Plot to ensure alignment bw points and rasters
# plot(projgrid)
# plot(shp.pts, add=TRUE)

## Parallelized extract: (larger datasets)
# rasterOptions(maxmemory = 4e+09)
# pts.ext <- DSMprops::parPTextr(sp = shp.pts, gridlist = cov.grids, os = "linux",nthreads = 50)
# ## Save points
# saveRDS(pts.ext, paste(predfolder,"/CONUS_nasis_SSURGO_SG100_covarsc.rds",sep=""))
## Updated extract for CONUS
pts.ext <- readRDS(paste(predfolder,"/CONUS_nasis_SSURGO_SG100_covarsc.rds",sep=""))

## Weed out duplicates
pts.ext$locid <- paste(pts.ext$x_std,pts.ext$y_std,sep="_")
pts.ext <- pts.ext[!duplicated(pts.ext$locid),]
pts.ext$tid <- "nasis"

## Join extracted pts to horizon data
pts.ext.hor <- left_join(pts@horizons[pts@horizons$peiid %in% pts.ext$peiid,],pts.ext, by="peiid")

## Prep nasis training data for Random Forest
pts.ext.hor$prop <- pts.ext.hor$sandfine_r ## UPDATE EVERY TIME
prop <- "sandfine_r" ## Dependent variable
hist(pts.ext.hor$prop)
summary(pts.ext.hor$prop)
## Set transformation and scaling: UPDATE EVERY TIME!!!!!!!!!!!!!!!!
trans <- "none" # none, log10, log, or sqrt
data_type <- "INT1U" # from raster::dataType - INT1U, INT1S, INT2S, INT2U, INT4S, INT4U, FLT4S, FLT8S
datastretch <- 1
datastretchlab <- paste(datastretch,"x",sep="")


##### Load and prep SCD data
## RSQlite workflow form https://github.com/ncss-tech/gsp-sas/blob/master/lab_data.Rmd
con <- dbConnect(RSQLite::SQLite(), "/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map/2020_Pedons/KSSL-snapshot-draft/KSSL-data.sqlite")
(ldm_names <- dbListTables(con))
ldm <- lapply(c("NCSS_Layer","NCSS_Site_Location","PSDA_and_Rock_Fragments"), function(x) dbReadTable(con , x))
names(ldm) <- c("NCSS_Layer","NCSS_Site_Location","PSDA_and_Rock_Fragments")
dbDisconnect(con)

## Now wrangle tables
scd.hor <- inner_join(ldm$NCSS_Layer,ldm$PSDA_and_Rock_Fragments[!duplicated(ldm$PSDA_and_Rock_Fragments$labsampnum),],by="labsampnum")
scd.pts <- ldm$NCSS_Site_Location

# ### SCD prep: Weed out points with imprecise coordinates ###
# scd.pts$latnchar <- nchar(abs(scd.pts$latitude_decimal_degrees))
# scd.pts$longnchar <- nchar(abs(scd.pts$longitude_decimal_degrees))
# scd.pts <- subset(scd.pts, scd.pts$latnchar > 5 & scd.pts$longnchar > 6)
## Location ID for later use and remove duplicates
scd.pts$locid <- paste(scd.pts$latitude_decimal,scd.pts$longitude_decima,sep="_")
scd.pts <- scd.pts[!duplicated(scd.pts$locid),] # Over 20k duplicates...

### Turn into spatial file
coordinates(scd.pts) <- ~ longitude_decima + latitude_decimal # Typos in sqlite KSSL snapshot...
temp.proj <- CRS("+proj=longlat +datum=WGS84") ## specify projection
projection(scd.pts) <- temp.proj
# ######## Clip with boundary if necessary ###########
scd.pts <- spTransform(scd.pts, cov.proj)
scd.pts <- scd.pts[polybound,]

## Further SCD prep
## Extract covariates for prediction onto SCD points
# scd.pts.ext <- DSMprops::parPTextr(scd.pts, cov.grids, os = "linux", nthreads=50)
# ## Save scd pts with extraction
# saveRDS(scd.pts.ext, paste(predfolder,"/SCD_PSD","_extracted_nasisSSURGO_SG100.rds",sep=""))
scd.pts.ext <- readRDS(paste(predfolder,"/SCD_PSD","_extracted_nasisSSURGO_SG100.rds",sep=""))

## Now join scd.hor to scd.pts and get rid of any duplicates
scd.pts.ext.hor <- inner_join(scd.hor,scd.pts.ext, by="site_key")
scd.pts.ext.hor$hzn_bot_locid <- paste(scd.pts.ext.hor$hzn_bot,scd.pts.ext.hor$locid,sep="_") # compound ID to remove horizon duplicates
scd.pts.ext.hor <- scd.pts.ext.hor[!duplicated(scd.pts.ext.hor$hzn_bot_locid),]

## SCD prep for RF
scd.pts.ext.hor$prop <- scd.pts.ext.hor$sand_f_psa ## UPDATE everytime!
scdprop <- "sand_f_psa"
scd.pts.ext.hor$tid <- "scd"

## Prep base raster for density calculations and spatial cross validation
# rasterOptions(maxmemory = 5e+08,chunksize = 5e+07)
# cpus <- 50 #detectCores() - 1
# beginCluster(cpus,type='SOCK')
# crs.img <- CRS(projection(projgrid))
# img10k <- projectRaster(projgrid,res=10000,method='ngb',crs=crs.img)
# endCluster()
# imgfun <- function(x){
#   ind <- ifelse(x>0,1,NA) ## need to make sure this makes sense for raster used
#   return(ind)
# }
# ## 10k grid with valus of 1 for all pixels except nodata areas
# img10kf <- calc(img10k,fun=imgfun,progress="text")
# ## Save file for future use
# saveRDS(img10kf,paste(predfolder,"/img10kf.rds",sep=""))
img10kf <- readRDS(paste(predfolder,"/img10kf.rds",sep=""))
## Create another grid with unique values for each pixel excepting nodata areas
## This is for use in the point density mapping
img10kfid <- img10kf
values(img10kfid) <- 1:ncell(img10kfid)
img10kfid <- overlay(stack(img10kfid,img10kf),fun=function(a,b){a*b})

##### Loop to train and predict properties for all depths
depths <- c(0,5,15,30,60,100,150)
for(d in depths){
  pts.extc <- subset(pts.ext.hor, as.numeric(pts.ext.hor$hzdept_r) <= d & as.numeric(pts.ext.hor$hzdepb_r) > d) # subset to chosen depth
  # pedonLocations <- unique(pts.extc$LocID) # if length differs from # of rows in pts, there are duplicates
  # pts.extc <- subset(pts.extc, !duplicated(pts.extc[c("LocID")])) #removes duplicates
  ptspred.list <- gsub(".tif","", cov.grids.names)# Take .tif off of the grid list to just get the variable names
  ptspred.list <- c(ptspred.list,"prop","mtchtype","peiid","tid","x_std","y_std") #Add dependent variable
  pts.extcc <- pts.extc[,c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use
  pts.extcc <- na.omit(pts.extcc)# Remove any record with NA's (in any column - be careful)
  ### Check for duplication with lab pedons
  scd.pts.d <- subset(scd.pts.ext.hor, as.numeric(scd.pts.ext.hor$hzn_top) <= d & as.numeric(scd.pts.ext.hor$hzn_bot) > d)
  ## Remove SCD dups
  scd.pts.d <- scd.pts.d[!duplicated(scd.pts.d$locid),]
  scd.pts.d$x_std <- scd.pts.d$longitude_decima
  scd.pts.d$y_std <- scd.pts.d$latitude_decimal
  scd.pts.d$peiid <- scd.pts.d$locid
  scd.pts.d$mtchtype <- "scd"
  scd.pts.d$tid <- "scd"
  ## Check for duplicated pedons betweeen NASIS and SCD using buffers
  coordinates(scd.pts.d) <- ~ longitude_decima + latitude_decimal
  projection(scd.pts.d) <- cov.proj
  pts.buf <- shp.pts[shp.pts$peiid %in% pts.extcc$peiid,]
  scd.pts.d.bufc <- raster::buffer(scd.pts.d,width=5,dissolve=F)
  pts.buf.exlude <- pts.buf[scd.pts.d.bufc,]
  pts.extcc <- pts.extcc[!pts.extcc$peiid %in% pts.buf.exlude@data$peiid,]
  pts.buf <- pts.buf[!pts.buf$peiid %in% pts.buf.exlude@data$peiid,c("peiid")]
  ## Now combine cleaned up datasets
  scd.pts.d <- scd.pts.d@data
  scd.pts.d <- scd.pts.d[,c(ptspred.list)]
  scd.pts.d <- na.omit(scd.pts.d)
  pts.extcc <- rbind(pts.extcc,scd.pts.d)
  coordinates(pts.extcc) <- ~ x_std + y_std
  projection(pts.extcc) <- cov.proj
  ## Characterize pt density for case weights
  polywin    <- as(polybound, "owin")
  pts.ppp    <- as(pts.buf, "ppp")
  ## Bound the data
  Window(pts.ppp) <- polywin
  ## Set the raster grid to align with
  pop  <- as.im.RasterLayer(img10kfid)
  E    <- tess(image=pop)  # Create a tesselated surface
  Q   <- quadratcount(pts.ppp, tess = E)  # Tally counts
  Q.d <- intensity(Q, image=TRUE) # create density raster
  ## Rescale to per km basis and convert to raster
  Qd.sgdf <- as.SpatialGridDataFrame.im(Q.d*((raster::res(img10kfid)[1])^2)) # convert to SGDF and scale to pts per cell
  Qd.rast <- raster(Qd.sgdf)
  writeRaster(Qd.rast, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",d,"_ptdensity_per_cell.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  names(Qd.rast) <- "ptspercell"
  pts.extcc <- extract(Qd.rast, pts.extcc, df=T, sp=T)
  ## Create spatial density based weights: better way to do this??
  pts.extcc@data$sp_wghts <- 1-(pts.extcc@data$ptspercell / (max(pts.extcc@data$ptspercell)*1.1))
  ## Data quality weights
  pts.extcc@data$qual_wts <- ifelse(pts.extcc@data$tid == "scd", 1.0,0.5)
  ## Combined sample weights: better way to do?
  pts.extcc@data$tot_wts <- pts.extcc@data$sp_wghts * pts.extcc@data$qual_wts
  ## Apply transformation
  if(trans=="log10") {pts.extcc$prop_t <- log10(pts.extcc$prop + 0.1)}
  if(trans=="log") {pts.extcc$prop_t <- log(pts.extcc$prop + 1)}
  if(trans=="sqrt") {pts.extcc$prop_t <- sqrt(pts.extcc$prop)}
  if(trans=="none") {pts.extcc$prop_t <- pts.extcc$prop}
  ## Save regression matrix file
  # saveRDS(pts.extcc, paste(predfolder,"/trainPTS_", prop, '_',d, "_cm_nasisSSURGO_SG100.rds",sep=""))
  pts.extcc <- readRDS(paste(predfolder,"/trainPTS_", prop, '_',d, "_cm_nasisSSURGO_SG100.rds",sep=""))
  ## Set up RF formula
  formulaStringRF <- as.formula(paste('prop_t ~', paste(gsub(".tif","", cov.grids.names), collapse="+")))
  ## Determine 95% interquartile range for relative prediction interval
  varrange <- as.numeric(quantile(pts.extcc@data$prop, probs=c(0.975), na.rm=T)-quantile(pts.extcc@data$prop, probs=c(0.025),na.rm=T)) ## TRANSFORM IF NEEDED!

  ############### Build quantile Random Forest
  ## Set training parameters
  trn.params <- list(ntrees = 100, min.node.size = 1)
  rf.qrf <- ranger(formulaStringRF, data=pts.extcc@data, num.trees = trn.params$ntrees, quantreg = T, num.threads = 50,
                   min.node.size = trn.params$min.node.size,
                   case.weights = pts.extcc@data$tot_wts)
  ## OOB error
  rf.qrf
  ## Linear Adjustment for bias in low and high predictions
  pts.extcc@data$trainpreds <- predict(rf.qrf, data=pts.extcc@data, num.threads = 50)$predictions
  attach(pts.extcc@data)
  rf_lm_adj <- lm(prop_t ~ trainpreds)
  detach(pts.extcc@data)
  pts.extcc@data$trainpredsadj <- predict(rf_lm_adj, newdata=pts.extcc@data)
  ## plot model performance stuff: How to get OOB preds??
  # plot(pts.extcc@data$prop_t~predict(rf.qrf,data=pts.extcc@data)$predictions) ## Fit plot
  x1 <-c(-100,0,100,10000,100000000)
  y1 <-c(-100,0,100,10000,100000000)
  # lines(x1,y1, col = 'red')#1:1 line
  # plot(pts.extcc@data$prop_t~pts.extcc@data$trainpreds) #Fit plot
  # lines(x1,y1, col = 'red')#1:1 line
  # plot(pts.extcc@data$prop_t~pts.extcc@data$trainpredsadj) #Adjusted Fit plot
  # lines(x1,y1, col = 'red')#1:1 line
  # TODO update for ranger variable importance
  #varImpPlot(rf.qrf)
  ## Save important documentation and intermediate files as R objects
  saveRDS(rf.qrf, paste(predfolder,"/rangerQRF_", prop, '_',d, "_cm_nasisSSURGO_SG100.rds",sep=""))
  saveRDS(rf_lm_adj, paste(predfolder,"/rflmadj_RFmodel_",prop,"_", d, "_cm_nasisSSURGO_SG100.rds",sep=""))
  ## Block to open prior files for updating work
  rf.qrf <- readRDS(paste(predfolder,"/rangerQRF_", prop, '_',d, "_cm_nasisSSURGO_SG100.rds",sep=""))
  rf_lm_adj <- readRDS(paste(predfolder,"/rflmadj_RFmodel_",prop,"_", d, "_cm_nasisSSURGO_SG100.rds",sep=""))



  ############ Cross Validate and examine metrics among Lab and Nasis pedons ##########
  ################### Manual Cross validation ################################
  pts.pcv <- DSMprops::CVranger(x = pts.extcc@data, fm = formulaStringRF, train.params = trn.params,
                                nfolds = 10, nthreads = 50, os = "linux") # 5min
  pts.pcv$valtype <- "pcv10f"

  pts.s1cv <- DSMprops::SpatCVranger(sp = pts.extcc, fm = formulaStringRF, train.params = trn.params,
                                     rast = img10kf, nfolds = 10, nthreads = 50, resol = 1, os = "linux") # 6 min
  pts.s1cv$valtype <- "s1cv10f"

  pts.s10cv <- DSMprops::SpatCVranger(sp = pts.extcc, fm = formulaStringRF, train.params = trn.params,
                                      rast = img10kf, nfolds = 10, nthreads = 50, resol = 10, os = "linux") # 6min
  pts.s10cv$valtype <- "s10cv10f"

  pts.s100cv <- DSMprops::SpatCVranger(sp = pts.extcc, fm = formulaStringRF, train.params = trn.params,
                                       rast = img10kf, nfolds = 10, nthreads = 50, resol = 100, os = "linux")
  pts.s100cv$valtype <- "s100cv10f"

  ## Combine CV tables and save in list as R object
  # cv.lst <- list(pts.pcv,pts.s1cv,pts.s10cv,pts.s100cv)
  # saveRDS(cv.lst,paste(predfolder,"/CVlist_", prop, '_',d, "_cm_nasisSSURGO_SG100.rds",sep="")) # takes forever...
  cv.lst <- readRDS(paste(predfolder,"/CVlist_", prop, '_',d, "_cm_nasisSSURGO_SG100.rds",sep=""))


  ## Validation metrics
  valmets <- DSMprops::valmetrics(xlst = cv.lst, trans = trans, varrange = varrange, prop = prop, depth = d)
  write.table(valmets, paste(predfolder,"/CVstats_", prop,"_", d, "_cm_nasisSSURGO_SG100.txt",sep=""), sep = "\t", row.names = FALSE)

  # ## PCV statistics: all data
  # cvp.RMSE = sqrt(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2, na.rm=TRUE))
  # cvp.Rsquared = 1-var(pts.extpcv$prop_t - pts.extpcv$pcvpred, na.rm=TRUE)/var(pts.extpcv$prop_t, na.rm=TRUE)
  # ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function: Using Duan's smearing estimator
  # if(trans=="log10") {pts.extpcv$pcvpred_bt <- (10^(pts.extpcv$pcvpred) - 0.1)*(mean(10^(pts.extpcv$prop_t - pts.extpcv$trainpredsadj)))}
  # if(trans=="log") {pts.extpcv$pcvpred_bt <- (exp(pts.extpcv$pcvpred) - 1)*(mean(exp(pts.extpcv$prop_t - pts.extpcv$trainpredsadj)))}
  # if(trans=="sqrt") {pts.extpcv$pcvpred_bt <- ((pts.extpcv$pcvpred)^2)*(mean((pts.extpcv$prop_t - pts.extpcv$trainpredsadj)^2))}
  # if(trans=="none") {pts.extpcv$pcvpred_bt <- pts.extpcv$pcvpred}
  # ## Untransformed calcs
  # cvp.RMSE_bt = sqrt(mean((pts.extpcv$prop - pts.extpcv$pcvpred_bt)^2, na.rm=TRUE))
  # cvp.Rsquared_bt = 1-var(pts.extpcv$prop - pts.extpcv$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv$prop, na.rm=TRUE)
  # ## PCV stats for scd points
  # pts.extpcv.scd <- subset(pts.extpcv, pts.extpcv$tid == "scd")
  # cvp.RMSE.scd <- sqrt(mean((pts.extpcv.scd$prop_t - pts.extpcv.scd$pcvpred)^2, na.rm=TRUE))
  # cvp.Rsquared.scd <- 1-var(pts.extpcv.scd$prop_t - pts.extpcv.scd$pcvpred, na.rm=TRUE)/var(pts.extpcv.scd$prop_t, na.rm=TRUE)
  # # cvp.RMSE.scdc <- sqrt(mean((pts.extpcv.scd$prop_tc - pts.extpcv.scd$pcvpred)^2, na.rm=TRUE))
  # # cvp.Rsquared.scdc <- 1-var(pts.extpcv.scd$prop_tc - pts.extpcv.scd$pcvpred, na.rm=TRUE)/var(pts.extpcv.scd$prop_t, na.rm=TRUE)
  # ## PCV stats for scd points: backtransformed
  # cvp.RMSE.scd_bt <- sqrt(mean((pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt)^2, na.rm=TRUE))
  # cvp.Rsquared.scd_bt <- 1-var(pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv.scd$prop, na.rm=TRUE)
  # ## Number of SCD samples
  # n_scd <- length(pts.extpcv.scd[,1])
  # ## RPI
  # ## Back transform low PI
  # # TODO Not sure if smearing estimator should be used on PIs? Probably not since the linear
  # # adjustment does not apply.
  # if(trans=="log10") {pts.extpcv$pcvpredpre.025_bt <- 10^(pts.extpcv$pcvpredpre.025) - 0.1}
  # if(trans=="log") {pts.extpcv$pcvpredpre.025_bt <- exp(pts.extpcv$pcvpredpre.025) - 1}
  # if(trans=="sqrt") {pts.extpcv$pcvpredpre.025_bt <- (pts.extpcv$pcvpredpre.025)^2}
  # if(trans=="none") {pts.extpcv$pcvpredpre.025_bt <- pts.extpcv$pcvpredpre.025}
  # ## Back transform Upper PI
  # if(trans=="log10") {pts.extpcv$pcvpredpre.975_bt <- 10^(pts.extpcv$pcvpredpre.975) - 0.1}
  # if(trans=="log") {pts.extpcv$pcvpredpre.975_bt <- exp(pts.extpcv$pcvpredpre.975) - 1}
  # if(trans=="sqrt") {pts.extpcv$pcvpredpre.975_bt <- (pts.extpcv$pcvpredpre.975)^2}
  # if(trans=="none") {pts.extpcv$pcvpredpre.975_bt <- pts.extpcv$pcvpredpre.975}
  # ## Calculate different metrics
  # pts.extpcv$abs.resid <- abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt)
  # pts.extpcv$RPI <- (pts.extpcv$pcvpredpre.975_bt - pts.extpcv$pcvpredpre.025_bt)/varrange
  # plot(pts.extpcv$abs.resid~pts.extpcv$RPI) # Quick look at relationship
  # ## Summarize RPI and residuals
  # ## Back transform original property to avoid bias in PICP
  # if(trans=="log10") {pts.extpcv$prop_bt <- 10^(pts.extpcv$prop_t) - 0.1}
  # if(trans=="log") {pts.extpcv$prop_bt <- exp(pts.extpcv$prop_t) - 1}
  # if(trans=="sqrt") {pts.extpcv$prop_bt <- (pts.extpcv$prop_t)^2}
  # if(trans=="none") {pts.extpcv$prop_bt <- pts.extpcv$prop_t}
  # pts.extpcv$rel.abs.resid <- pts.extpcv$abs.resid/varrange
  # RPI.cvave <- mean(pts.extpcv$RPI)
  # RPI.cvmed <- median(pts.extpcv$RPI)
  # rel.abs.res.ave <- mean(pts.extpcv$rel.abs.resid)
  # rel.abs.res.med <- median(pts.extpcv$rel.abs.resid)
  # pts.extpcv$BTbias <- pts.extpcv$prop_bt - pts.extpcv$prop
  # BTbias.abs.max <- max(abs(pts.extpcv$BTbias))
  # BTbias.ave <- mean(pts.extpcv$BTbias)
  # PICP <- sum(ifelse(pts.extpcv$prop_bt <= pts.extpcv$pcvpredpre.975_bt & pts.extpcv$prop_bt >= pts.extpcv$pcvpredpre.025_bt,1,0))/length(pts.extpcv[,1])
  # ## Create PCV table
  # CVdf <- data.frame(cvp.RMSE, cvp.Rsquared, cvp.RMSE_bt, cvp.Rsquared_bt, cvp.RMSE.scd, cvp.Rsquared.scd,  cvp.RMSE.scd_bt, cvp.Rsquared.scd_bt,n_scd,RPI.cvave,RPI.cvmed,PICP,rel.abs.res.ave,rel.abs.res.med,BTbias.abs.max,BTbias.ave)
  # names(CVdf) <- c("cvp.RMSE","cvp.Rsquared","cvp.RMSE_bt", "cvp.Rsquared_bt", "cvp.RMSE.scd", "cvp.Rsquared.scd", "cvp.RMSE.scd_bt", "cvp.Rsquared.scd_bt","n_scd","RPI.CVave","RPI.CVmed","PICP","rel.abs.res.ave","rel.abs.res.med","BTbias.abs.max","BTbias.ave")
  # write.table(CVdf, paste(predfolder,"/PCVstats_", prop,"_", d, "_cm_nasisSSURGO_SG100.txt",sep=""), sep = "\t", row.names = FALSE)

  ## Save Cross validation graph and data for future plotting
  saveRDS(pts.extpcv, paste(predfolder,"/cvlm_preds_2D_",prop, "_", d, "_cm_nasisSSURGO_SG100.rds", sep=""))

  ###### Cross validation plots
  ## All data
  viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
  gplt.dcm.2D.CV <- ggplot(data=pts.extpcv, aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1) + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri)) +
    ggtitle(paste("Cross val", prop, d, "cm",sep=" "))
  gplt.dcm.2D.CV
  ggsave(paste(predfolder,'/ValPlot_1to1_all_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  ## Now just SCD pedons
  gplt.dcm.2D.CV.SCD <- ggplot(data=pts.extpcv.scd, aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri)) +
    ggtitle(paste("Cross val", prop, d, "cm",sep=" "))
  gplt.dcm.2D.CV.SCD
  ggsave(paste(predfolder,'/ValPlot_1to1_scd_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV.SCD, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))


  ############################## Raster Preditions ######################################################
  ##### Reference covar rasters to use in prediction
  rasters <- stack(cov.grids)
  #names(rasters)

  ## Ranger Predict
  predfun <- function(model, ...) predict(model, ...)$predictions
  # TODO figure out why internal ranger parallelization did not seem to work.
  # predtst <- predict(rasters, rf.qrf, fun=predfun, type="response", progress="text", num.threads = 50)


  ## Predict onto covariate grid
  ## Parallelized predict
  rasterOptions(maxmemory = 5e+09,chunksize = 9e+08)# maxmemory = 1e+09,chunksize = 1e+08 for soilmonster
  beginCluster(50,type='SOCK')
  Sys.time()
  predl <- clusterR(rasters, predict, args=list(model=rf.qrf, fun=predfun,type = "quantiles", quantiles = c(0.025)),progress="text")
  Sys.time()
  predh <- clusterR(rasters, predict, args=list(model=rf.qrf, fun=predfun,type = "quantiles", quantiles = c(0.975)),progress="text")
  Sys.time()
  #pred <- predict(rasters, rf.qrf, fun=predfun, progress="text")
  pred <- clusterR(rasters, predict, args=list(model=rf.qrf, fun=predfun), progress="text") # works ok
  #pred <- clusterR(rasters, predict, args=list(model=soiclass),progress="text")
  Sys.time()
  ## End cluster to clear cache and memory
  # endCluster()
  # gc()
  # beginCluster(50,type='SOCK')
  ## Rename for lm adjustment
  names(pred) <- "trainpreds"
  # ## Linear Adjustment
  predlm <- clusterR(pred, predict, args=list(model=rf_lm_adj),progress="text")
  ## Scale function for saving rasters as integers
  rastscale.fn <- function(x) {
    ind <- x*datastretch
    ind <- ifelse(ind<0,0,ind)
    return(ind)
  }
  ## Logic statement to handle transformed versus non transformed
  if(trans=="none"){ # For models that have not been transformed
    # ## PI widths
    s <- stack(predh,predl)
    # PIwidth.fn <- function(a,b) {
    #   ind <- a-b
    #   return(ind)
    # }
    # PIwidth <- clusterR(s, overlay, args=list(fun=PIwidth.fn),progress = "text")
    # # Determine 95% interquantile range of original training data for horizons that include the depth being predicted
    PIrelwidth.fn <- function(a,b) {
      ind <- ((a-b)/varrange)*1000
      return(ind)
    }
    PIrelwidth <- clusterR(s, overlay, args=list(fun=PIrelwidth.fn),progress = "text",export='varrange')
    ## Rescale rasters for saving
    predsc <- clusterR(pred, calc, args=list(fun=rastscale.fn),progress='text',export='datastretch')
    predlmsc <- clusterR(predlm, calc, args=list(fun=rastscale.fn),progress='text',export='datastretch')
    predlsc <- clusterR(predl, calc, args=list(fun=rastscale.fn),progress='text',export='datastretch')
    predhsc <- clusterR(predh, calc, args=list(fun=rastscale.fn),progress='text',export='datastretch')
    ## Write rasters
    writeRaster(predsc, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRF.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
    writeRaster(predlmsc, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRFadj.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
    writeRaster(predlsc, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRF_95PI_l.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
    writeRaster(predhsc, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRF_95PI_h.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
    writeRaster(PIrelwidth, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_1000x_",d,"_cm_2D_QRF_95PI_relwidth.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT2U', progress="text")
    # # writeRaster(PIwidth, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",d,"_cm_2D_QRF_95PI_width.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
  } else {
    ## Back transformation Stuff
    if(trans=="log10"){
      smrest <- mean(10^(pts.extcc@data$prop_t - pts.extcc@data$trainpredsadj))
      bt.fnlm <- function(x) { # Duans smearing est
      ind <-  ((10^(x))-0.1)*smrest
      return(ind)
      }
      bt.fn <- function(x) {
        ind <-  (10^(x))-0.1
        return(ind)
      }
    }
    if(trans=="log"){
      smrest <- mean(exp(pts.extcc@data$prop_t - pts.extcc@data$trainpredsadj))
      bt.fnlm <- function(x) { # Duans smearing est
      ind <-  ((exp(x))-1)*smrest
      return(ind)
      }
      bt.fn <- function(x) {
        ind <-  (exp(x))-1
        return(ind)
      }
      }
    if(trans=="sqrt"){
      smrest <- mean((pts.extcc@data$prop_t - pts.extcc@data$trainpredsadj)^2)
      bt.fnlm <- function(x) { # Duans smearing est
      ind <-  (x^2)*smrest
      return(ind)
      }
      bt.fn <- function(x) {
        ind <-  x^2
        return(ind)
      }
      }
    predh_bt <- clusterR(predh, calc, args=list(fun=bt.fn),progress='text')
    predl_bt <- clusterR(predl, calc, args=list(fun=bt.fn),progress='text')
    pred_bt <- clusterR(pred, calc, args=list(fun=bt.fn),progress='text')
    predlm_bt <- clusterR(predlm, calc, args=list(fun=bt.fnlm),progress='text',export='smrest')
    ## Stack PIs for RPI calculation
    s_bt <- stack(predh_bt,predl_bt)
    # PIwidth_bt.fn <- function(a,b) {
    #   ind <- a-b
    #   return(ind)
    # }
    # PIwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIwidth_bt.fn),progress = "text")
    ## If transformed, use the following code for PI width prep steps
    PIrelwidth_bt.fn <- function(a,b) {
      ind <- ((a-b)/varrange)*100
      ind <- ifelse(ind>65533,65533,ind)
      return(ind)
    }
    PIrelwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIrelwidth_bt.fn),progress = "text", export='varrange')
    ## Rescale rasters for saving
    pred_btsc <- clusterR(pred_bt, calc, args=list(fun=rastscale.fn),progress='text',export='datastretch')
    predlm_btsc <- clusterR(predlm_bt, calc, args=list(fun=rastscale.fn),progress='text',export='datastretch')
    predl_btsc <- clusterR(predl_bt, calc, args=list(fun=rastscale.fn),progress='text',export='datastretch')
    predh_btsc <- clusterR(predh_bt, calc, args=list(fun=rastscale.fn),progress='text',export='datastretch')
    ## Write rasters
    writeRaster(pred_btsc, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRF_bt.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
    writeRaster(predlm_btsc, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRFadj_bt.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
    writeRaster(predl_btsc, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRF_95PI_l_bt.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
    writeRaster(predh_btsc, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRF_95PI_h_bt.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
    #writeRaster(PIwidth_bt, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",datastretchlab,"_",d,"_cm_2D_QRF_95PI_width_bt.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
    writeRaster(PIrelwidth_bt, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_100x_",d,"_cm_2D_QRF_95PI_relwidth_bt.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT2U', progress="text")
  }
  ## Close out raster cluster
  endCluster()
}


################### Manual Cross validation ################################
# pts.extcvm <- pts.extcc
# nfolds <- 10
# pts.extcvm$folds <- sample.int(nfolds,size =length(pts.extcvm[,1]),replace=T)
# pts.extcvm$prop_t <- pts.extcvm$prop ## UPDATE: tranform if needed else just create new version of prop
# formulaStringCVm <- as.formula(paste('prop_t ~', paste(gsub(".tif","", cov.grids), collapse="+")))
# #for (g in seq(nfolds)){
# CV_factorRF <- function(g,pts.extcvm, formulaStringCVm){
#   traindf <- subset(pts.extcvm, pts.extcvm$folds != g)
#   testdf <- subset(pts.extcvm, pts.extcvm$folds == g)
#   xtrain.t <- as.matrix(traindf[c(gsub(".tif","", cov.grids))])
#   ytrain.t <- c(as.matrix(traindf$prop_t))
#   rf.pcv <- quantregForest(x=xtrain.t, y=ytrain.t, importance=TRUE, ntree=100, keep.forest=TRUE)
#   rf.pcvc <- rf.pcv
#   class(rf.pcvc) <- "randomForest"
#   traindf$pcvpredpre <- predict(rf.pcvc, newdata=traindf)
#   testdf$pcvpredpre <- predict(rf.pcvc, newdata=testdf)
#   #traindf$pcvpredpre <- predict(rf.pcv, newdata=traindf, what=c(0.5)) ## If median is desired
#   #testdf$pcvpredpre <- predict(rf.pcv, newdata=testdf,, what=c(0.5)) ## If median is desired
#   testdf$pcvpredpre.025 <- predict(rf.pcv, newdata=testdf, what=c(0.025))
#   testdf$pcvpredpre.975 <- predict(rf.pcv, newdata=testdf, what=c(0.975))
#   attach(traindf)
#   lm.pcv <- lm(prop_t~pcvpredpre)
#   detach(traindf)
#   testdf$pcvpred <- predict(lm.pcv, newdata=testdf)
#   return(testdf)
# }
# snowfall::sfInit(parallel=TRUE, cpus=nfolds)
# snowfall::sfExport("pts.extcvm","formulaStringCVm","CV_factorRF","cov.grids")
# snowfall::sfLibrary(randomForest)
# snowfall::sfLibrary(quantregForest)
# pts.extpcv <- snowfall::sfLapply(1:nfolds, function(g){CV_factorRF(g, pts.extcvm=pts.extcvm,formulaStringCVm=formulaStringCVm)})
# snowfall::sfStop()
# pts.extpcv <- plyr::rbind.fill(pts.extpcv)
# pts.extpcv$pcvpred = as.numeric(pts.extpcv$pcvpred)

## Validate with lab pedons
# scd.pts.d <- subset(scd.pts, as.numeric(scd.pts$hzn_top) <= d & as.numeric(scd.pts$hzn_bot) > d)
# scd.pts.d <- spTransform(scd.pts.d, projection(pred))
# pts.extpcv <- extract(pred, scd.pts.d, df=TRUE, sp=T)
# #pts.extpcv@data$prop_t <- log10(pts.extpcv@data$ec_satp + 0.1)
# pts.extpcv@data$prop <- pts.extpcv@data$ec_satp
# pts.extpcv@data$pcvpred <- pts.extpcv@data$layer
# #pts.extpcv@data$pcvpred <- log10(pts.extpcv@data$layer + 0.1)
# ## PCV statistics
# cvp.RMSE = sqrt(mean((pts.extpcv@data$prop_t - pts.extpcv@data$pcvpred)^2, na.rm=TRUE))
# cvp.Rsquared = 1-var(pts.extpcv@data$prop_t - pts.extpcv@data$pcvpred, na.rm=TRUE)/var(pts.extpcv@data$prop_t, na.rm=TRUE)
# ## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function
# #pts.extpcv@data$pcvpred_bt <- pts.extpcv@data$pcvpred
# cvp.RMSE_bt = sqrt(mean((pts.extpcv@data$prop - pts.extpcv@data$pcvpred)^2, na.rm=TRUE))
# cvp.Rsquared_bt = 1-var(pts.extpcv@data$prop - pts.extpcv@data$pcvpred, na.rm=TRUE)/var(pts.extpcv@data$prop, na.rm=TRUE)
# ## PCV stats for scd points
# # pts.extpcv.scd <- subset(pts.extpcv, pts.extpcv$tid == "scd")
# # cvp.RMSE.scd <- sqrt(mean((pts.extpcv.scd$prop_t - pts.extpcv.scd$pcvpred)^2, na.rm=TRUE))
# # cvp.Rsquared.scd <- 1-var(pts.extpcv.scd$prop_t - pts.extpcv.scd$pcvpred, na.rm=TRUE)/var(pts.extpcv.scd$prop_t, na.rm=TRUE)
# # ## PCV stats for scd points: backtransformed
# # cvp.RMSE.scd_bt <- sqrt(mean((pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt)^2, na.rm=TRUE))
# # cvp.Rsquared.scd_bt <- 1-var(pts.extpcv.scd$prop - pts.extpcv.scd$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv.scd$prop, na.rm=TRUE)
# # ## Number of SCD samples
# # n_scd <- length(pts.extpcv.scd[,1])
# # ## RPI
# # pts.extpcv$prop_bt <- pts.extpcv$prop_t # UPDATE: backtransform if necessary. Used for PICP and to characterize backtransformation bias
# # pts.extpcv$pcvpredpre.025_bt <- pts.extpcv$pcvpredpre.025 # UPDATE: backtransform if necessary
# # pts.extpcv$pcvpredpre.975_bt <- pts.extpcv$pcvpredpre.975 # UPDATE: backtransform if necessary
# # pts.extpcv$abs.resid <- abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt)
# # pts.extpcv$RPI <- (pts.extpcv$pcvpredpre.975_bt - pts.extpcv$pcvpredpre.025_bt)/varrange
# # plot(pts.extpcv$abs.resid~pts.extpcv$RPI) # Quick look at relationship
# ## Summarize RPI and residuals
# # pts.extpcv$rel.abs.resid <- pts.extpcv$abs.resid/varrange
# # RPI.cvave <- mean(pts.extpcv$RPI)
# # RPI.cvmed <- median(pts.extpcv$RPI)
# # rel.abs.res.ave <- mean(pts.extpcv$rel.abs.resid)
# # rel.abs.res.med <- median(pts.extpcv$rel.abs.resid)
# # pts.extpcv$BTbias <- pts.extpcv$prop_bt - pts.extpcv$prop
# # BTbias.abs.max <- max(abs(pts.extpcv$BTbias))
# # BTbias.ave <- mean(pts.extpcv$BTbias)
# # PICP <- sum(ifelse(pts.extpcv$prop_bt <= pts.extpcv$pcvpredpre.975_bt & pts.extpcv$prop_bt >= pts.extpcv$pcvpredpre.025_bt,1,0))/length(pts.extpcv[,1])
# # ## Create PCV table
# # CVdf <- data.frame(cvp.RMSE, cvp.Rsquared, cvp.RMSE_bt, cvp.Rsquared_bt, cvp.RMSE.scd, cvp.Rsquared.scd, cvp.RMSE.scd_bt, cvp.Rsquared.scd_bt,n_scd,RPI.cvave,RPI.cvmed,PICP,rel.abs.res.ave,rel.abs.res.med,BTbias.abs.max,BTbias.ave)
# # names(CVdf) <- c("cvp.RMSE","cvp.Rsquared","cvp.RMSE_bt", "cvp.Rsquared_bt", "cvp.RMSE.scd", "cvp.Rsquared.scd", "cvp.RMSE.scd_bt", "cvp.Rsquared.scd_bt","n_scd","RPI.CVave","RPI.CVmed","PICP","rel.abs.res.ave","rel.abs.res.med","BTbias.abs.max","BTbias.ave")
# # setwd(predfolder)
# # write.table(CVdf, paste("PCVstats", prop, d, "cm_nasisSSURGO_ART_SG100.txt",sep="_"), sep = "\t", row.names = FALSE)
# CVdf <- data.frame(cvp.RMSE, cvp.Rsquared)
# names(CVdf) <- c("cvp.RMSE","cvp.Rsquared")
# #setwd(predfolder)
# write.table(CVdf, paste(predfolder,"/PCVstats_", prop, "_", d, "_cm_nasisSSURGO_ART_SG100.txt",sep=""), sep = "\t", row.names = FALSE)
# # plot(pts.extpcv$prop~pts.extpcv$pcvpred_bt)
# # plot(pts.extpcv$prop~pts.extpcv$pcvpred_bt, xlim=c(0,10),ylim=c(0,10))
# # plot(pts.extpcv$prop~pts.extpcv$pcvpred_bt, xlim=c(0,5),ylim=c(0,5))
# # plot(pts.extpcv$prop~pts.extpcv$pcvpred_bt, xlim=c(0,1),ylim=c(0,1))
# # plot(pts.extpcv$prop_t~pts.extpcv$pcvpred)
# #lines(x1,y1, col = 'red')#1:1 line
# # CV plots
# viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
# gplt.dcm.2D.CV <- ggplot(data=pts.extpcv@data, aes(prop, pcvpred)) +
#   stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + xlim(-5,105) + ylim(-5,105) +
#   theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
#   xlab("Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri)) +
#   ggtitle(paste("Cross val", prop, d, "cm",sep=" "))
# gplt.dcm.2D.CV
# ggsave(paste(predfolder,'/ValPlot_1to1_scd_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
# ## Save Cross validation graph and data for future plotting
# saveRDS(pts.extpcv, paste(prop, "cvlm_preds_2D", d, "cm_nasisSSURGO_ART_SG100.rds", sep="_"))



############# Masking water pixels out ############
# nlcd <- raster("/home/tnaum/data/UCRB_Covariates/NLCDcl.tif")
# beginCluster(30,type='SOCK')
## Make a mask raster
# mask_fn <- function(nlcd){ind <- ifelse(nlcd!=11,1,NA)
#   return(ind)
# }
# mask <- clusterR(nlcd, calc, args=list(fun=mask_fn),progress='text')
# endCluster()
# plot(mask)
# writeRaster(mask, overwrite=TRUE,filename="/home/tnaum/data/BLMsoils/nlcd_watermask.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text",datatype='INT1U')
# rm(mask)
## Now set up a list of rasters and function to mask out water
# rasterOptions(maxmemory = 1e+09,chunksize = 1e+08)
# setwd("/home/tnaum/data/BLMsoils/Sand_2D_NASIS_SSURGO_SCD")
# grids <- list.files(pattern=".tif$")
# mskfn <- function(rast,mask){
#   ind <- rast*mask
#   ind[ind<0]<-0 # to bring the slighly negative predictions back to zero
#   return(ind)
# }
# ## par list apply fn
# watermask_fn <- function(g){
#   setwd("/home/tnaum/data/BLMsoils/Sand_2D_NASIS_SSURGO_SCD")
#   rast <- raster(g)
#   names(rast) <- "rast"
#   setwd("/home/tnaum/data/BLMsoils")
#   mask <- raster("/home/tnaum/data/BLMsoils/nlcd_watermask.tif")
#   h2ostk <- stack(rast,mask)
#   setwd("/home/tnaum/data/BLMsoils/Sand_2D_NASIS_SSURGO_SCD/masked")
#   overlay(h2ostk,fun=mskfn,progress='text',filename=g, options=c("COMPRESS=DEFLATE", "TFW=YES"))
#   gc()
# }
# snowfall::sfInit(parallel=TRUE, cpus=30)
# snowfall::sfExport("watermask_fn","mskfn","grids")
# snowfall::sfLibrary(raster)
# Sys.time()
# snowfall::sfLapply(grids, function(g){watermask_fn(g)})
# Sys.time()
# snowfall::sfStop()



