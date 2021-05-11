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
                       "hexbin","doParallel","aqp","Hmisc","spatstat","maptools","DSMprops","RSQLite","lubridate","stringr")
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
## Geographic coordinate quality levels
pts_geocode <- readRDS("/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map/2020_Pedons/geocode_weighting.RDS") # From Dave White
# pts_geocode_wt <- readRDS("/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map/2020_Pedons/geocode_weighting_v2.RDS") # From Dave White
pts_geocode$geo_wt <- pts_geocode$wt
pts_geocode$wt <- NULL
pts_geocode$peiid <- as.character(pts_geocode$peiid)
# ## Pedon quality clases and weights
# pts_qual_cls <- readRDS("/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map/2020_Pedons/pedonquality_weighting.RDS") # Dave White
# pts_qual_cls <- pts_qual_cls[!duplicated(pts_qual_cls$peiid),]
## Pedons with extracted SSUGO component data
pts <- readRDS(paste(ptsfolder,"/NASIS_all_component_horizon_match_SPC_ssurgo20.rds",sep=""))
pts.proj <- proj4string(pts)
n.pts <- pts@site
## Bring in orig nasis points to eliminate those with partial coords
load("/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map/2020_Pedons/nasis_sites_20210325.RData")
s$latnchar <- nchar(abs(s$y_std))
s$longnchar <- nchar(abs(s$x_std))
s <- subset(s, s$latnchar > 5 & s$longnchar > 6)
okpeiids <- s$peiid
n.pts <- n.pts[n.pts$peiid %in% okpeiids,]
## Create sp object
coordinates(n.pts) <- ~ x_std + y_std
projection(n.pts) <- pts.proj
n.pts@data <- n.pts@data[,c("peiid","mtchtype","compname")]


######### Grid Prep #################
## Make list of grids
cov.grids <- list.files(path = covfolder,pattern=".tif$",full.names = T, recursive = F)
cov.grids.names <- basename(cov.grids)
## If points need to be matched up to grids ###
projgrid <- raster(paste(covfolder,"/T07PRI5.tif",sep="")) # This raster is fully clipped
cov.proj <- projection(projgrid)
n.pts <- spTransform(n.pts, CRS(cov.proj)) # project to match rasters

####### Polygon boundary if needed to clip down
polybound <- readOGR("/ped/GIS_Archive/US_boundaries/states_21basic", "conus_bound")
polybound <- spTransform(polybound, cov.proj)
n.pts <- n.pts[polybound,]

####### Now create a random sample of points in CONUS to use for gRPI estimation
pts.gRPI <- spsample(polybound[1,], 1000000, type = 'random')
pts.gRPI <- SpatialPointsDataFrame(pts.gRPI, data.frame(row.names=row.names(pts.gRPI), ID=1:length(pts.gRPI)))

## Plot to ensure alignment bw points and rasters
# plot(projgrid)
# plot(n.pts, add=TRUE)
# plot(pts.gRPI, add=TRUE)

## Parallelized extract for gRPI points
# rasterOptions(maxmemory = 4e+09)
# pts.gRPI <- DSMprops::parPTextr(sp = pts.gRPI, gridlist = cov.grids, os = "linux",nthreads = 50)
# pts.gRPI <- na.omit(pts.gRPI)
# ## Save points
# saveRDS(pts.gRPI, paste(ptsfolder,"/CONUS_random_gRPIsamp.rds",sep=""))
## Updated extract for CONUS
pts.gRPI <- readRDS(paste(ptsfolder,"/CONUS_random_gRPIsamp.rds",sep=""))


## Parallelized extract for nasis points: (larger datasets)
# rasterOptions(maxmemory = 4e+09)
# pts.ext <- DSMprops::parPTextr(sp = n.pts, gridlist = cov.grids, os = "linux",nthreads = 50)
# ## Save points
# saveRDS(pts.ext, paste(predfolder,"/CONUS_nasis_extracted.rds",sep=""))
## Updated extract for CONUS
pts.ext <- readRDS(paste(predfolder,"/CONUS_nasis_extracted.rds",sep=""))
pts.ext <- left_join(pts.ext, pts_geocode, by = "peiid")
pts.ext$geo_wt <- ifelse(is.na(pts.ext$geo_wt),8,pts.ext$geo_wt) # 18.8k NAs put in class 8

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
ldm <- lapply(c("NCSS_Layer","NCSS_Site_Location","PSDA_and_Rock_Fragments","NCSS_Pedon_Taxonomy"), function(x) dbReadTable(con , x))
names(ldm) <- c("NCSS_Layer","NCSS_Site_Location","PSDA_and_Rock_Fragments","NCSS_Pedon_Taxonomy")
dbDisconnect(con)

## Now wrangle tables
scd.hor <- inner_join(ldm$NCSS_Layer,ldm$PSDA_and_Rock_Fragments[!duplicated(ldm$PSDA_and_Rock_Fragments$labsampnum),],by="labsampnum")
ldm$NCSS_Pedon_Taxonomy$peiid <- ldm$NCSS_Pedon_Taxonomy$pedoniid
scd.pts <- ldm$NCSS_Site_Location
scd.pts <- left_join(scd.pts,ldm$NCSS_Pedon_Taxonomy[ldm$NCSS_Pedon_Taxonomy$site_key %in% scd.pts$site_key, c("site_key","peiid")], by="site_key")

# ### SCD prep: Weed out points with imprecise coordinates ###
scd.pts$latnchar <- nchar(abs(scd.pts$latitude_decimal))
scd.pts$longnchar <- nchar(abs(scd.pts$longitude_decima))
scd.pts <- subset(scd.pts, scd.pts$latnchar > 5 & scd.pts$longnchar > 6)
## Location ID for later use and remove duplicates
scd.pts$locid <- paste(scd.pts$latitude_decimal,scd.pts$longitude_decima,sep="_")
scd.pts <- scd.pts[!duplicated(scd.pts$locid),] # Over 20k duplicates...
scd.pts <- scd.pts[!duplicated(scd.pts$peiid),]

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
# saveRDS(scd.pts.ext, paste(predfolder,"/SCD","_extracted.rds",sep=""))
scd.pts.ext <- readRDS(paste(predfolder,"/SCD","_extracted.rds",sep=""))
## Create geo-coordinate quality weights
scd.pts.ext$date <- as.Date(scd.pts.ext$site_obsdate, format = "%m/%d/%Y")
scd.pts.ext$decdate <- lubridate::decimal_date(scd.pts.ext$date)
scd.pts.ext$peiid <- as.character(scd.pts.ext$peiid)
scd.pts.ext <- left_join(scd.pts.ext, pts_geocode, by="peiid")
scd.pts.ext$geo_wt <- ifelse(is.na(scd.pts.ext$geo_wt) & scd.pts.ext$date >= "2010-01-01", 6, scd.pts.ext$geo_wt)
scd.pts.ext$geo_wt <- ifelse(is.na(scd.pts.ext$geo_wt) & scd.pts.ext$date >= "2005-01-01", 7, scd.pts.ext$geo_wt)
scd.pts.ext$geo_wt <- ifelse(is.na(scd.pts.ext$geo_wt), 8, scd.pts.ext$geo_wt) ## Awful lot of class 8s...

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
## Bring in feature space weights reference distributions
ft_wts_ref <- readRDS("/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map/2020_Pedons/ref_df.RDS") # From Stephen Roecker, 5/3/2021

##### Loop to train and predict properties for all depths
depths <- c(0,5,15,30,60,100,150)
for(d in depths){
  pretime <- Sys.time()
  pts.extc <- subset(pts.ext.hor, as.numeric(pts.ext.hor$hzdept_r) <= d & as.numeric(pts.ext.hor$hzdepb_r) > d) # subset to chosen depth
  ptspred.list <- gsub(".tif","", cov.grids.names)# Take .tif off of the grid list to just get the variable names
  ptspred.list <- c(ptspred.list,"prop","mtchtype","peiid","tid","x_std","y_std","geo_wt","locid") #Add dependent variable
  pts.extcc <- pts.extc[,c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use
  pts.extcc <- na.omit(pts.extcc)# Remove any record with NA's (in any column - be careful)
  ### Check for duplication with lab pedons
  scd.pts.d <- subset(scd.pts.ext.hor, as.numeric(scd.pts.ext.hor$hzn_top) <= d & as.numeric(scd.pts.ext.hor$hzn_bot) > d)
  ## Remove SCD dups
  scd.pts.d <- scd.pts.d[!duplicated(scd.pts.d$locid),]
  scd.pts.d$x_std <- scd.pts.d$longitude_decima
  scd.pts.d$y_std <- scd.pts.d$latitude_decimal
  scd.pts.d$mtchtype <- "scd"
  scd.pts.d$tid <- "scd"
  ## Check for duplicated pedons betweeen NASIS and SCD using buffers
  coordinates(scd.pts.d) <- ~ longitude_decima + latitude_decimal
  projection(scd.pts.d) <- cov.proj
  pts.buf <- n.pts[n.pts$peiid %in% pts.extcc$peiid,]
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
  ## Geocode classes
  pts.extcc@data <- dplyr::mutate(pts.extcc@data, geo_cls = ifelse(geo_wt<6,"gps",NA),
                                  geo_cls = ifelse(geo_wt==6,"gps2",geo_cls),
                                  geo_cls = ifelse(geo_wt==7,"gps3",geo_cls),
                                  geo_cls = ifelse(geo_wt==8,"unk",geo_cls))
  #pts.extcc@data <- left_join(pts.extcc@data, pts_qual_cls[,c("peiid","source")], by="peiid")
  ## Apply transformation
  if(trans=="log10") {pts.extcc$prop_t <- log10(pts.extcc$prop + 0.1)}
  if(trans=="log") {pts.extcc$prop_t <- log(pts.extcc$prop + 1)}
  if(trans=="sqrt") {pts.extcc$prop_t <- sqrt(pts.extcc$prop)}
  if(trans=="none") {pts.extcc$prop_t <- pts.extcc$prop}
  ## Save pts available for modeling file
  saveRDS(pts.extcc, paste(predfolder,"/AvailPTS_", prop, '_',d, "_cm.rds",sep=""))
  #pts.extcc <- readRDS(paste(predfolder,"/AvailPTS_", prop, '_',d, "_cm.rds",sep=""))
  ## Set up RF formula
  formulaStringRF <- as.formula(paste('prop_t ~', paste(gsub(".tif","", cov.grids.names), collapse="+")))

  ############ Cross Validate and examine metrics among Lab and Nasis pedons ##########
  ################### Manual Cross validation ################################
  ## Set training parameters for trials
  trn.params <- list(ntrees = 100, min.node.size = 1)

  #### summarize RPI in full raster prediction using sample for speed (tested against full average)
  ## Determine 95% interquartile range for relative prediction interval
  varrange_gRPI <- as.numeric(quantile(pts.extcc@data$prop, probs=c(0.975), na.rm=T)-quantile(pts.extcc@data$prop, probs=c(0.025),na.rm=T)) ## TRANSFORM IF NEEDED!
  gRPI_rf <- ranger(formulaStringRF, data=pts.extcc@data, num.trees = trn.params$ntrees, quantreg = T, num.threads = 60,
                    min.node.size = trn.params$min.node.size)
  ## Predict onto random gRPI sample pts
  pts.gRPI$lowpred <- predict(gRPI_rf, data=pts.gRPI, num.threads = 60,type = "quantiles", quantiles = c(0.025))$predictions
  pts.gRPI$highpred <- predict(gRPI_rf, data=pts.gRPI, num.threads = 60,type = "quantiles", quantiles = c(0.975))$predictions
  pts.gRPI$RPI <- (pts.gRPI@data$highpred - pts.gRPI@data$lowpred) / varrange_gRPI
  gRPI.ave <- mean(pts.gRPI@data$RPI)
  gRPI.med <- median(pts.gRPI@data$RPI)
  gRPI.n <- length(rpi_samp)
  RPIg_df <- data.frame(gRPI.ave,gRPI.med,gRPI.n)
  write.table(RPIg_df, paste(predfolder,"/gRPI_", prop,"_", d, "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
  ## Match best CV scheme using gRPI
  gRPIall <- ave(gRPI.ave,gRPI.med)


  ## Normal 10-fold cross validation
  cv10f <- DSMprops::CVranger(x = pts.pcv@data, fm = formulaStringRF, train.params = trn.params,
                              nfolds = 10, nthreads = 60, os = "linux") # 5min
  ## Spatial 10-fold cross validation on 1km blocks
  s1cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,
                                    rast = img10kf, nfolds = 10, nthreads = 60, resol = 1, os = "linux") # 6 min
  ## Spatial 10-fold cross validation on 10km blocks
  s10cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,
                                     rast = img10kf, nfolds = 10, nthreads = 60, resol = 10, os = "linux") # 6min
  ## Spatial 10-fold cross validation on 30km blocks
  s50cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,
                                     rast = img10kf, nfolds = 10, nthreads = 60, resol = 50, os = "linux") # 6min
  ## Spatial 10-fold cross validation on 100km blocks
  s100cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF, train.params = trn.params,
                                      rast = img10kf, nfolds = 10, nthreads = 60, resol = 100, os = "linux")

  ## Combine CV tables and save in list as R object
  cv.lst <- list(cv10f,s1cv10f,s10cv10f,s50cv10f, s100cv10f)
  #saveRDS(cv.lst,paste(predfolder,"/CVlist_", prop, '_',d, "_cm.rds",sep="")) # takes forever...
  #cv.lst <- readRDS(paste(predfolder,"/CVlist_", prop, '_',d, "_cm.rds",sep=""))

  ## Validation metrics for CVs at different spatial supports
  valmets_sCV <- DSMprops::valmetrics(xlst = cv.lst, trans = trans, prop = prop, depth = d)
  valmets_sCV$RPIall <- (valmets_sCV$RPI.cvave + valmets_sCV$RPI.cvmed) / 2
  idx_val <- which(abs(valmets_sCV$RPIall-gRPIall)==min(abs(valmets_sCV$RPIall-gRPIall)))
  bestval <- valmets_sCV[idx_val,c("valtype")]
  bestvalpts <- get(bestval)
  valmets_sCV$bestval <- bestval
  ## Pick best CV resolution based on RPI match
  resol_rpi <- as.numeric(gsub("s", "", str_split(bestval,"cv")[[1]][1]))
  ## Save results
  write.table(valmets_sCV, paste(predfolder,"/sCVstats_", prop,"_", d, "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
  saveRDS(bestvalpts,paste(predfolder,"/CVfull_best_pts_", prop, '_',d, "_cm.rds",sep=""))


  ## Now prepare a data tier grid search using the chosen cross validation approach
  geo_vec <- c("gps","gps_gps2","gps_gps2_gps3", "gps_gps2_gps3_unk")
  srce_vec <- c("scd","scd_direct","scd_direct_home","scd_direct_home_adjacent")
  grid_vec <- expand.grid(geo=geo_vec, srce=srce_vec, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
  CV_grid_fn <- function(x){
    levs <- data.frame(grid_vec[x,])
    colnames(levs) <- colnames(grid_vec)
    geo_levs <- str_split(levs$geo, "_")[[1]]
    srce_levs <- str_split(levs$srce, "_")[[1]]
    ptseval <- pts.extcc
    ptseval <- ptseval[ptseval@data$geo_cls %in% geo_levs,]
    ptseval <- ptseval[ptseval@data$mtchtype %in% srce_levs,]
    nfolds <- 10
    resol <- resol_rpi # kilometers
    ptsevalcv <- DSMprops::SpatCVranger(sp = ptseval, fm = formulaStringRF, train.params = trn.params,
                                        rast = img10kf, nfolds = nfolds, nthreads = 60, resol = resol, os = "linux",casewts = "tot_wts")
    ptsevalcv$cvgrid <- paste(colnames(levs),levs[1,],collapse="_",sep="_")
    return(ptsevalcv)
  }
  ## List apply implementation
  Sys.time()
  CV_grid_lst <- lapply(1:nrow(grid_vec),CV_grid_fn)
  Sys.time()
  #saveRDS(CV_grid_lst,paste(predfolder,"/CV_grid_list_", prop, '_',d, "_cm.rds",sep="")) # takes forever
  ## Subset CVs to recent scd points
  CV_grid_lst_scd_gps <- lapply(1:length(CV_grid_lst),function(x){
    df <- CV_grid_lst[[x]]
    df <- df[df$mtchtype=="scd"&(df$geo_cls=="gps"|df$geo_cls=="gps2"),]
    return(df)
  })
  CV_grid_lst_scd_all <- lapply(1:length(CV_grid_lst),function(x){
    df <- CV_grid_lst[[x]]
    df <- df[df$mtchtype=="scd",]
    return(df)
  })
  ## Validation metrics for CVs with different case weighting schemes
  Grid_valmets <- DSMprops::valmetrics(xlst = CV_grid_lst, trans = trans, prop = prop, depth = d)
  Grid_valmets_scd_gps <- DSMprops::valmetrics(xlst = CV_grid_lst_scd_gps, trans = trans, prop = prop, depth = d)
  Grid_valmets_scd_all <- DSMprops::valmetrics(xlst = CV_grid_lst_scd_all, trans = trans, prop = prop, depth = d)
  valmets_rank <- data.frame( cvgrid = Grid_valmets$cvgrid, Rsq.all = Grid_valmets$Rsq, RPI.all = Grid_valmets$RPI.cvave,
                              Rsq.scd = Grid_valmets_scd_all$Rsq, Rsq.gscd = Grid_valmets_scd_gps$Rsq, QRMSE_bt.gscd = Grid_valmets$QRMSE_bt,
                             QMedAE_bt.gscd = Grid_valmets_scd_gps$QMedAE_bt, n = Grid_valmets$n)
  grid_valmets_lst <- list(Grid_valmets,Grid_valmets_scd_gps,Grid_valmets_scd_all,valmets_rank)
  names(grid_valmets_lst) <- c("Grid_valmets","Grid_valmets_scd_gps","Grid_valmets_scd_all","valmets_rank")
  saveRDS(grid_valmets_lst,paste(predfolder,"/CV_grid_valmets_", prop, '_',d, "_cm.rds",sep=""))
  ## Now pick model mode ranking on multiple criteria
  valmets_rank$RPI.cvave.rnk <- rank(valmets_rank$RPI.all,ties.method = "average")
  valmets_rank$QRMSE_bt.gscd.rnk <- (nrow(valmets_rank)+1) - rank(valmets_rank$QRMSE_bt.gscd,ties.method = "average")
  valmets_rank$QMedAE_bt.gscd.rnk <- (nrow(valmets_rank)+1) - rank(valmets_rank$QMedAE_bt.gscd,ties.method = "average")
  valmets_rank$Rsq.scd.rnk <- (nrow(valmets_rank)+1) - rank(valmets_rank$Rsq.scd,ties.method = "average")
  valmets_rank$Rsq.gscd.rnk <- (nrow(valmets_rank)+1) - rank(valmets_rank$Rsq.gscd,ties.method = "average")
  valmets_rank$Rsq.all.rnk <- (nrow(valmets_rank)+1) - rank(valmets_rank$Rsq.all,ties.method = "average")
  valmets_rank$rank_ave <- apply(valmets_rank[,c("RPI.cvave.rnk","QRMSE_bt.gscd.rnk","QMedAE_bt.gscd.rnk",
                                                 "Rsq.scd.rnk","Rsq.gscd.rnk","Rsq.all.rnk")], MARGIN = 1, FUN=mean)
  valmets_rank$rank_final <- rank(valmets_rank$rank_ave,ties.method = "random")
  write.table(valmets_rank, paste(predfolder,"/Grid_valmets_", prop,"_", d, "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
  ## Select points from overall dataset based on combined rankings
  model_params <- str_split(valmets_rank[valmets_rank$rank_final==min(valmets_rank$rank_final),c("cvgrid")][1],"_")[[1]]
  srce_params <- model_params[(match("srce",model_params)+1):length(model_params)]
  geo_params <- model_params[(match("geo",model_params)+1):(match("srce",model_params)-1)]
  ## Subset training points to new optimized dataset
  pts.pcv <- pts.extcc[pts.extcc$mtchtype %in% srce_params,]
  pts.pcv <- pts.pcv[pts.pcv$geo_cls %in% geo_params,]

  ## Characterize pt density for case weights
  polywin    <- as(polybound, "owin")
  pts.ppp    <- as(pts.pcv, "ppp")
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
  pts.pcv <- extract(Qd.rast, pts.pcv, df=T, sp=T)
  ## Create spatial density based weights: better way to do this??
  pts.pcv@data$sp_wts <- (1-(pts.pcv@data$ptspercell / (max(pts.pcv@data$ptspercell)*1.1)))^5
  ## Feature space weights
  vars <- c("SLPNED6", "POSNED6", "EVI", "PPT", "TEMP")
  p1 <- seq(0, 1, 0.1)
  p2 <- seq(0, 1, 0.5)
  probs <- list(p1, p2, p1, p1, p1, NULL)
  names(probs) <- vars
  pts.pcv@data$feat_wts <- DSMprops::feat_wts(ref_df = ft_wts_ref, obs_df = pts.pcv@data, vars, probs)$wts # ~700 NAs for d = 0
  pts.pcv@data$feat_wts <- ifelse(is.na(pts.pcv@data$feat_wts), median(pts.pcv@data$feat_wts, na.rm=T), pts.pcv@data$feat_wts) # Stopgap for NAs created
  pts.pcv@data$feat_wts <- pts.pcv@data$feat_wts / max(pts.pcv@data$feat_wts, na.rm = T) # Scale 0-1
  ## Combined sample weights: better way to do?
  pts.pcv@data$tot_wts <- pts.pcv@data$sp_wts * pts.pcv@data$feat_wts

  ######## Test for case weights schemes
  ## Now prepare a case.weights grid search using the chosen cross validation approach
  feat_vec <- c("full","mid","none")
  sp_vec <- c("full","mid","none")
  grid_vec <- expand.grid(feat=feat_vec, spat=sp_vec, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
  CV_wts_grid_fn <- function(x){
    levs <- data.frame(grid_vec[x,])
    colnames(levs) <- colnames(grid_vec)
    ptseval <- pts.pcv
    if(levs$feat == "mid"){ptseval@data$feat_wts <- ptseval@data$feat_wts^(1/7)}
    if(levs$feat == "none"){ptseval@data$feat_wts <- 1}
    if(levs$spat == "mid"){ptseval@data$sp_wts <- sqrt(ptseval@data$sp_wts)}
    if(levs$spat == "none"){ptseval@data$sp_wts <- 1}
    nfolds <- 10
    resol <- resol_rpi # kilometers
    ptseval@data$tot_wts <- ptseval@data$feat_wts * ptseval@data$sp_wts
    # ptsevalsamp <- sample(ptseval@data$peiid, size = floor(0.75 * nrow(ptseval@data)), prob = ptseval@data$tot_wts)
    # ptseval <- ptseval[ptseval@data$peiid %in% ptsevalsamp,]
    # ptseval@data$tot_wts <- NULL
    # ptsevalcv <- DSMprops::SpatCVranger(sp = ptseval, fm = formulaStringRF, train.params = trn.params,
    #                                     rast = img10kf, nfolds = nfolds, nthreads = 60, resol = resol, os = "linux")
    ptsevalcv <- DSMprops::SpatCVranger(sp = ptseval, fm = formulaStringRF, train.params = trn.params,
                                        rast = img10kf, nfolds = nfolds, nthreads = 60, resol = resol, os = "linux", casewts = "tot_wts")
    ptsevalcv$cvgrid <- paste(colnames(levs),levs[1,],collapse="_",sep="_")
    return(ptsevalcv)
  }
  ## List apply implementation
  Sys.time()
  CV_wts_grid_lst <- lapply(1:nrow(grid_vec),CV_wts_grid_fn)
  Sys.time()
  #saveRDS(CV_wts_grid_lst,paste(predfolder,"/CV_wts_grid_lst_", prop, '_',d, "_cm.rds",sep="")) # takes forever
  ## Validation metrics for CVs with different case weighting schemes
  ## Subset CVs to recent scd points
  CV_wts_grid_lst_scd_gps <- lapply(1:length(CV_wts_grid_lst),function(x){
    df <- CV_wts_grid_lst[[x]]
    df <- df[df$mtchtype=="scd"&(df$geo_cls=="gps"|df$geo_cls=="gps2"),]
    return(df)
  })
  CV_wts_grid_lst_scd_all <- lapply(1:length(CV_wts_grid_lst),function(x){
    df <- CV_wts_grid_lst[[x]]
    df <- df[df$mtchtype=="scd",]
    return(df)
  })
  ## Validation metrics for CVs with different case weighting schemes
  Grid_wt_valmets <- DSMprops::valmetrics(xlst = CV_wts_grid_lst, trans = trans, prop = prop, depth = d)
  Grid_wt_valmets_scd_gps <- DSMprops::valmetrics(xlst = CV_wts_grid_lst_scd_gps, trans = trans, prop = prop, depth = d)
  Grid_wt_valmets_scd_all <- DSMprops::valmetrics(xlst = CV_wts_grid_lst_scd_all, trans = trans, prop = prop, depth = d)
  wt_valmets_rank <- data.frame( cvgrid = Grid_wt_valmets$cvgrid, Rsq.all = Grid_wt_valmets$Rsq, RPI.all = Grid_wt_valmets$RPI.cvave,
                              Rsq.scd = Grid_wt_valmets_scd_all$Rsq, Rsq.gscd = Grid_wt_valmets_scd_gps$Rsq, QRMSE_bt.gscd = Grid_wt_valmets$QRMSE_bt,
                              QMedAE_bt.gscd = Grid_wt_valmets_scd_gps$QMedAE_bt, n = Grid_wt_valmets$n)
  Grid_wt_valmets_lst <- list(Grid_wt_valmets,Grid_wt_valmets_scd_gps,Grid_wt_valmets_scd_all,wt_valmets_rank)
  names(Grid_wt_valmets_lst) <- c("Grid_wt_valmets","Grid_wt_valmets_scd_gps","Grid_wt_valmets_scd_all","wt_valmets_rank")
  saveRDS(Grid_wt_valmets_lst,paste(predfolder,"/CV_wts_grid_wt_valmets_", prop, '_',d, "_cm.rds",sep=""))
  ## Now pick model mode ranking on multiple criteria
  wt_valmets_rank$RPI.cvave.rnk <- rank(wt_valmets_rank$RPI.all,ties.method = "average")
  wt_valmets_rank$QRMSE_bt.gscd.rnk <- (nrow(wt_valmets_rank)+1) - rank(wt_valmets_rank$QRMSE_bt.gscd,ties.method = "average")
  wt_valmets_rank$QMedAE_bt.gscd.rnk <- (nrow(wt_valmets_rank)+1) - rank(wt_valmets_rank$QMedAE_bt.gscd,ties.method = "average")
  wt_valmets_rank$Rsq.scd.rnk <- (nrow(wt_valmets_rank)+1) - rank(wt_valmets_rank$Rsq.scd,ties.method = "average")
  wt_valmets_rank$Rsq.gscd.rnk <- (nrow(wt_valmets_rank)+1) - rank(wt_valmets_rank$Rsq.gscd,ties.method = "average")
  wt_valmets_rank$Rsq.all.rnk <- (nrow(wt_valmets_rank)+1) - rank(wt_valmets_rank$Rsq.all,ties.method = "average")
  wt_valmets_rank$rank_ave <- apply(wt_valmets_rank[,c("RPI.cvave.rnk","QRMSE_bt.gscd.rnk","QMedAE_bt.gscd.rnk",
                                                 "Rsq.scd.rnk","Rsq.gscd.rnk","Rsq.all.rnk")], MARGIN = 1, FUN=mean)
  wt_valmets_rank$rank_final <- rank(wt_valmets_rank$rank_ave,ties.method = "random")
  write.table(wt_valmets_rank, paste(predfolder,"/Grid_wt_valmets_", prop,"_", d, "_cm.txt",sep=""), sep = "\t", row.names = FALSE)

  ## Pick out CV matrix to save and plot with
  idx_wts_val <- which(wt_valmets_rank$rank_final==1)
  #bestwtval <- wt_valmets_rank[idx_val,c("cvgrid")]
  bestvalwtpts <- CV_wts_grid_lst[[idx_wts_val]]
  saveRDS(bestvalwtpts,paste(predfolder,"/CV_finalwts_best_pts_", prop, '_',d, "_cm.rds",sep=""))

  ## Now set up final case weights based on grid search
  model_params <- str_split(wt_valmets_rank[wt_valmets_rank$rank_final==min(wt_valmets_rank$rank_final),c("cvgrid")][1],"_")[[1]]
  spat_param <- model_params[(match("spat",model_params)+1):length(model_params)]
  feat_param <- model_params[(match("feat",model_params)+1):(match("spat",model_params)-1)]
  ## Adjust weights according to chosen model
  pts.pcv@data <- dplyr::mutate(pts.pcv@data, sp_wts = ifelse(spat_param == "mid",sqrt(sp_wts),sp_wts),
                                sp_wts = ifelse(spat_param == "none",1,sp_wts)) # spatial wts
  pts.pcv@data <- dplyr::mutate(pts.pcv@data, feat_wts = ifelse(feat_param == "mid",feat_wts^(1/7),feat_wts),
                                feat_wts = ifelse(feat_param == "none",1,feat_wts)) # feat wts
  ## Save training points file
  saveRDS(pts.pcv, paste(predfolder,"/TrainPTS_", prop, '_',d, "_cm.rds",sep=""))
  #pts.pcv <- readRDS(paste(predfolder,"/TrainPTS_", prop, '_',d, "_cm.rds",sep=""))

  ###### Cross validation plots
  ## All data
  viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
  scaleFUN <- function(x) round(x,0)
  gplt.dcm.2D.CV <- ggplot(data=bestvalwtpts, aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1) + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste(bestval,"Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV
  ggsave(paste(predfolder,'/ValPlot_1to1_all_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  ## Now just SCD pedons
  gplt.dcm.2D.CV.SCD <- ggplot(data=bestvalwtpts[bestvalwtpts$mtchtype=="scd",], aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("SCD Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste(bestval,"Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV.SCD
  ggsave(paste(predfolder,'/ValPlot_1to1_scd_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV.SCD, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  ## Now just with SCD pedons with GPS or after 2010
  gplt.dcm.2D.CV.gSCD <- ggplot(data=bestvalwtpts[bestvalwtpts$mtchtype=="scd"&(bestvalwtpts$geo_cls=="gps"|bestvalwtpts$geo_cls=="gps2"),], aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("gps SCD Measured") + ylab("CV Prediction") + scale_fill_gradientn( name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste(bestval,"Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV.gSCD
  ggsave(paste(predfolder,'/ValPlot_1to1_gscd_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV.gSCD, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))


  ## TODO ? model tuning step for mtry, min node size, ntrees, other params???

  ############### Build quantile Random Forest
  ## TODO incorporate updated weighting or tuning parameters???????
  ## Determine 95% interquartile range for relative prediction interval
  varrange <- as.numeric(quantile(pts.pcv@data$prop, probs=c(0.975), na.rm=T)-quantile(pts.pcv@data$prop, probs=c(0.025),na.rm=T)) ## TRANSFORM IF NEEDED!
  ## Train global ranger model
  rf.qrf <- ranger(formulaStringRF, data=pts.pcv@data, num.trees = trn.params$ntrees, quantreg = T, num.threads = 60,
                   min.node.size = trn.params$min.node.size,
                   case.weights = pts.pcv$tot_wts)
  ## OOB error
  # rf.qrf
  ## Linear Adjustment for bias in low and high predictions
  pts.pcv@data$trainpreds <- predict(rf.qrf, data=pts.pcv@data, num.threads = 60)$predictions
  attach(pts.pcv@data)
  rf_lm_adj <- lm(prop_t ~ trainpreds)
  detach(pts.pcv@data)
  pts.pcv@data$trainpredsadj <- predict(rf_lm_adj, newdata=pts.pcv@data)
  ## Save important documentation and intermediate files as R objects
  saveRDS(rf.qrf, paste(predfolder,"/rangerQRF_", prop, '_',d, "_cm.rds",sep=""))
  saveRDS(rf_lm_adj, paste(predfolder,"/rflmadj_RFmodel_",prop,"_", d, "_cm.rds",sep=""))
  ## Block to open prior files for updating work
  # rf.qrf <- readRDS(paste(predfolder,"/rangerQRF_", prop, '_',d, "_cm.rds",sep=""))
  # rf_lm_adj <- readRDS(paste(predfolder,"/rflmadj_RFmodel_",prop,"_", d, "_cm.rds",sep=""))


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
  rasterOptions(maxmemory = 6e+09,chunksize = 9e+08)# maxmemory = 1e+09,chunksize = 1e+08 for soilmonster
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
    PIrelwidth <- clusterR(s_bt, overlay, args=list(fun=PIrelwidth_bt.fn),progress = "text", export='varrange')
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
    writeRaster(PIrelwidth, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_100x_",d,"_cm_2D_QRF_95PI_relwidth_bt.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT2U', progress="text")
  }
  ## Close out raster cluster
  endCluster()


  ## Wrap up loop
  posttime <- Sys.time()
  runtime <- posttime - pretime
  print(paste(d, " cm was done at", posttime,"in",runtime, sep=" "))
}


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



