######################
## Random Forest script that includes:
## Extraction of covariates to points
## Prediction interval creation
## Cross Validation
## Most steps parallelized
######################


# Workspace setup
# Install packages if not already installed
required.packages <- c("raster", "sp", "rgdal", "ranger","caret", "snow", "snowfall", "dplyr", "ggplot2",
                       "hexbin","doParallel","aqp","Hmisc","spatstat","maptools","DSMprops","RSQLite","lubridate","stringr")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package: Windows only
#memory.limit(500000)
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)
cpus <- min(detectCores()-2, 124)

## Key Folder Locations
predfolder <- "/mnt/disks/sped/solus100preds/v2tst_gnat_pts/Rock_tot_gRPI_250k"
repofolder <- "/mnt/disks/sped/repos/DSMprops"
covfolder <- "/mnt/disks/sped/covs100m"
propcovfldr <- "/mnt/disks/sped/covs100m_by_prop"
ptsfolder <- "/mnt/solus100/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final"
pedonfldr <- "/mnt/solus100/2020_Pedons"
ssurgo_fgdb <- "/home/tnaum/data/gSSURGO20/gSSURGO_CONUS.gdb"

######## Load soil profile collection ##############
## Geographic coordinate quality levels
pts_geocode <- readRDS(paste0(pedonfldr,"/geocode_weighting.RDS")) # From Dave White
pts_geocode$geo_wt <- pts_geocode$wt
pts_geocode$wt <- NULL
pts_geocode$peiid <- as.character(pts_geocode$peiid)
## Pedons with extracted SSUGO component data
pts <- readRDS(paste(ptsfolder,"/NASIS_all_component_horizon_match_SPC_ssurgo20.rds",sep=""))
pts.proj <- proj4string(pts)
n.pts <- pts@site
## Bring in orig nasis points to eliminate those with partial coords
load(paste0(pedonfldr,"/nasis_sites_20210325.RData")) # object s
s$latnchar <- nchar(abs(s$y_std))
s$longnchar <- nchar(abs(s$x_std))
s <- subset(s, s$latnchar > 5 & s$longnchar > 6)
okpeiids <- s$peiid
n.pts <- n.pts[n.pts$peiid %in% okpeiids,]
## Create sp object
coordinates(n.pts) <- ~ x_std + y_std
projection(n.pts) <- pts.proj
n.pts@data <- n.pts@data[,c("peiid","mtchtype","compname")]


### Define SSURGO property for modeling and covariate layer selection
prop <- "fragvol_r" ## Dependent variable: UPDATE EVERY TIME


######### Grid Prep #################
## Make list of grids
cov.grids <- list.files(path = covfolder,pattern=".tif$",full.names = T, recursive = F)
cov.grids.names <- basename(cov.grids)

## Grab gNATSGO property maps to use in creating point sample
prop.cov.grids <- list.files(path = paste0(propcovfldr,"/",prop),pattern=prop,full.names = T, recursive = F)
prop.cov.grids <- prop.cov.grids[grepl(".tif",prop.cov.grids)]
prop.cov.grids <- prop.cov.grids[!grepl("_NAs_",prop.cov.grids)]
prop.cov.grids.names <- basename(prop.cov.grids)

## Extraction list for gNATSGO training points
gNATpts_cov.grids <- append(cov.grids,prop.cov.grids)

## If points need to be matched up to grids ###
projgrid <- raster(paste(covfolder,"/T07PRI5.tif",sep="")) # This raster is fully clipped
cov.proj <- projection(projgrid)
n.pts <- spTransform(n.pts, CRS(cov.proj)) # project to match rasters

####### Polygon boundary if needed to clip down
polybound <- readOGR("/mnt/solus100/gis", "conus_bound_500k")
polybound <- spTransform(polybound, cov.proj)
n.pts <- n.pts[polybound,]

####### Now create a random sample of points in CONUS to use for gRPI estimation
# pts.gRPI <- spsample(polybound[1,], 1000000, type = 'random')
# pts.gRPI <- SpatialPointsDataFrame(pts.gRPI, data.frame(row.names=row.names(pts.gRPI), ID=1:length(pts.gRPI)))
# saveRDS(pts.gRPI, paste0(ptsfolder,"/CONUS_random_gRPIsamp.rds"))
# pts.gRPI <- readRDS(paste0(ptsfolder,"/CONUS_random_gRPIsamp.rds"))

## Plot to ensure alignment bw points and rasters
# plot(projgrid)
# plot(n.pts, add=TRUE)
# plot(pts.gRPI, add=TRUE)

## Parallelized extract for gRPI points
# rasterOptions(maxmemory = 1.5e+10)
# pts.gRPI <- DSMprops::parPTextr(sp = pts.gRPI, gridlist = cov.grids, os = "windows",nthreads = cpus)
# pts.gRPI <- na.omit(pts.gRPI)
# ## Save points
# saveRDS(pts.gRPI, paste(ptsfolder,"/CONUS_random_gRPIsamp_covs.rds",sep=""))
## Updated extract for CONUS
pts.gRPI <- readRDS(paste(ptsfolder,"/CONUS_random_gRPIsamp_covs.rds",sep=""))

## Set up gNATSGO property estimate training points
# pts.gNAT <- spsample(polybound[1,], 250000, type = 'random')
# pts.gNAT <- SpatialPointsDataFrame(pts.gNAT, data.frame(row.names=row.names(pts.gNAT), ID=1:length(pts.gNAT)))
# saveRDS(pts.gNAT, paste0(ptsfolder,"/CONUS_random_gNATsamp_250k.rds"))
# pts.gNAT <- readRDS(paste0(ptsfolder,"/CONUS_random_gNATsamp_250k.rds"))
# rasterOptions(maxmemory = 1.5e+10)
# pts.gNAT <- DSMprops::parPTextr(sp = pts.gNAT, gridlist = gNATpts_cov.grids, os = "windows",nthreads = cpus)
## Save points
# saveRDS(pts.gNAT, paste(ptsfolder,"/CONUS_random_gNATsamp_250k_covs_",prop,"_",".rds",sep=""))
pts.gNAT <- readRDS(paste(ptsfolder,"/CONUS_random_gNATsamp_250k_covs_",prop,"_",".rds",sep=""))

## Parallelized extract for nasis points: (larger datasets)
# rasterOptions(maxmemory = 2e+10)
# pts.ext <- DSMprops::parPTextr(sp = n.pts, gridlist = cov.grids, os = "windows",nthreads = cpus)
# ## Save points
# saveRDS(pts.ext, paste(ptsfolder,"/CONUS_nasis_extracted_covs.rds",sep=""))
## Updated extract for CONUS
pts.ext <- readRDS(paste(ptsfolder,"/CONUS_nasis_extracted_covs.rds",sep=""))
pts.ext <- left_join(pts.ext, pts_geocode, by = "peiid")
s$peiid <- as.character(s$peiid)
pts.ext <- left_join(pts.ext, s[,c("peiid","obsdate","obsdatekind")], by = "peiid")
pts.ext$obsdate <- as.Date(pts.ext$obsdate)
pts.ext$geo_wt <- ifelse(pts.ext$geo_wt==8 & pts.ext$obsdate > "2000-01-01" & pts.ext$obsdatekind == "actual site observation date",
                         7,pts.ext$geo_wt)
pts.ext$geo_wt <- ifelse(is.na(pts.ext$geo_wt),8,pts.ext$geo_wt) # 18.8k NAs put in class 8

## Weed out duplicates
pts.ext$locid <- paste(pts.ext$x_std,pts.ext$y_std,sep="_")
pts.ext <- pts.ext[!duplicated(pts.ext$locid),]
pts.ext$tid <- "nasis"

## Join extracted pts to horizon data
pts.ext.hor <- left_join(pts@horizons[pts@horizons$peiid %in% pts.ext$peiid,],pts.ext, by="peiid")

#### Rock fragment prep for SSURGO extracted NASIS pts
## Read in gdb with rgdal to get sp object
# chfrags <- sf::st_read(dsn = ssurgo_fgdb, layer = "chfrags")
# chfrags[] <- lapply(chfrags, function(x) if (is.factor(x)) as.character(x) else {x})
# chkeys <- unique(pts.ext.hor$chkey)
# chfrags <- chfrags[chfrags$chkey %in% chkeys,]
# saveRDS(chfrags,paste0(predfolder,"/chfrags_conus.rds"))
chfrags <- readRDS(paste0(predfolder,"/chfrags_conus.rds"))
## Sum together rock frag vols for all recorded rock sizes for each horizon
frags_chkey <- plyr::ddply(chfrags,~chkey,summarise,fragvol_r = sum(fragvol_r))
pts.ext.hor <- left_join(pts.ext.hor,frags_chkey,by="chkey")
pts.ext.hor$fragvol_r <- ifelse(is.na(pts.ext.hor$fragvol_r), 0, pts.ext.hor$fragvol_r) # NAs set to zero as these are not recorded in chfrags table

##### Bring in original NASIS data to use those estimates as primary source of rock fragments
## NOTE: NOT USING SCD DATA FOR ROCK DUE TO LACK OF >3" DATA
load(paste0(pedonfldr,"/nasis_pedons_20210325.RData")) # object named spc
# nasis_ped_pts <- spc@site[spc@site$peiid %in% okpeiids,c("peiid","x_std","y_std")]
# nasis_ped_pts <- subset(nasis_ped_pts,!is.na(nasis_ped_pts$x_std))
# ## Create sp object
# coordinates(nasis_ped_pts) <- ~ x_std + y_std
# proj4string(nasis_ped_pts) <- CRS(SRS_string = "EPSG:4326")
# nasis_ped_pts <- spTransform(nasis_ped_pts, cov.proj) # project to match rasters
# nasis_ped_pts <- nasis_ped_pts[polybound,]
# ## Extract covariates
# ## Parallelized extract for nasis points: (larger datasets)
# rasterOptions(maxmemory = 2e+09)
# nas.pedpts.ext <- DSMprops::parPTextr(sp = nasis_ped_pts, gridlist = cov.grids, os = "linux",nthreads = 62)
# ## Save points
# saveRDS(nas.pedpts.ext, paste(predfolder,"/CONUS_nasis_pedpts_extracted_spatcovs.rds",sep=""))
## Updated extract for CONUS
nas.pedpts.ext <- readRDS(paste(predfolder,"/CONUS_nasis_pedpts_extracted_spatcovs.rds",sep=""))
## Now attach geo weights
nas.pedpts.ext <- left_join(nas.pedpts.ext, pts_geocode, by = "peiid")
nas.pedpts.ext <- left_join(nas.pedpts.ext, s[,c("peiid","obsdate","obsdatekind")], by = "peiid")
nas.pedpts.ext$obsdate <- as.Date(nas.pedpts.ext$obsdate)
nas.pedpts.ext$geo_wt <- ifelse(nas.pedpts.ext$geo_wt==8 & nas.pedpts.ext$obsdate > "2000-01-01" & nas.pedpts.ext$obsdatekind == "actual site observation date",
                                7,nas.pedpts.ext$geo_wt)
nas.pedpts.ext$geo_wt <- ifelse(is.na(nas.pedpts.ext$geo_wt),8,nas.pedpts.ext$geo_wt) # 18.8k NAs put in class 8


## List of pedons to exclude from SCD data with errors in texture mod versus rock volume estimates
badpeiids <- read.csv(paste0(pedonfldr,"/nomatch_texmod_peiid.csv"),stringsAsFactors = F)
## remove bad peiids from NASIS field data
nas.pedpts.ext <- nas.pedpts.ext[!nas.pedpts.ext$peiid %in% badpeiids$nomatch_peiid,]
## Now remove pts from nasis-ssurgo data that already have a field estimate
pts.ext.hor <- pts.ext.hor[!pts.ext.hor$peiid %in% nas.pedpts.ext$peiid,]

## Prep nasis training data for Random Forest
pts.ext.hor$prop <- pts.ext.hor[,prop]  ## UPDATE EVERY TIME
hist(pts.ext.hor$prop)
summary(pts.ext.hor$prop)
## Set transformation and scaling: UPDATE EVERY TIME!!!!!!!!!!!!!!!!
trans <- "log" # none, log10, log, or sqrt
data_type <- "INT1U" # from raster::dataType - INT1U, INT1S, INT2S, INT2U, INT4S, INT4U, FLT4S, FLT8S
datastretch <- 1
datastretchlab <- paste(datastretch,"x",sep="")
pts.ext.hor$prop <- ifelse(pts.ext.hor$prop > 100, 100, pts.ext.hor$prop) # Can't have higher than 100% rock

##### Load and prep SCD data: Note, for rock, field observations (nas.pedpts.ext) will be the equivalent of 'SCD' data as they are our best estimates
scd.pts.ext <- nas.pedpts.ext
scd.pts.ext$locid <- paste(scd.pts.ext$y_std,scd.pts.ext$x_std,sep="_")
scd.pts.ext <- scd.pts.ext[!duplicated(scd.pts.ext$locid),] # Over 20k duplicates...
scd.pts.ext <- scd.pts.ext[!duplicated(scd.pts.ext$peiid),]

## Now join scd.hor to scd.pts and get rid of any duplicates
scd.pts.ext.hor <- inner_join(spc@horizons,scd.pts.ext, by="peiid")
scd.pts.ext.hor$hzn_bot_locid <- paste(scd.pts.ext.hor$hzdepb,scd.pts.ext.hor$locid,sep="_") # compound ID to remove horizon duplicates
scd.pts.ext.hor <- scd.pts.ext.hor[!duplicated(scd.pts.ext.hor$hzn_bot_locid),]

## SCD prep for RF
scd.pts.ext.hor$prop <- scd.pts.ext.hor$total_frags_pct_nopf  ## UPDATE everytime!
scdprop <- "total_frags_pct_nopf"
scd.pts.ext.hor$tid <- "scd"
hist(scd.pts.ext.hor$prop)
summary(scd.pts.ext.hor$prop)
scd.pts.ext.hor$prop <- ifelse(scd.pts.ext.hor$prop > 100, 100, scd.pts.ext.hor$prop) # Can't have > 100% rock

## Change NASIS terms to SCD terms for loop function
scd.pts.ext.hor$hzn_top <- scd.pts.ext.hor$hzdept
scd.pts.ext.hor$hzn_bot <- scd.pts.ext.hor$hzdepb

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
# saveRDS(img10kf,paste(pedonfldr,"/img10kf.rds",sep=""))
img10kf <- readRDS(paste(pedonfldr,"/img10kf.rds",sep=""))
## Create another grid with unique values for each pixel excepting nodata areas
## This is for use in the point density mapping
img10kfid <- img10kf
values(img10kfid) <- 1:ncell(img10kfid)
img10kfid <- overlay(stack(img10kfid,img10kf),fun=function(a,b){a*b})

##### Loop to train and predict properties for all depths
depths <- c(0,5,15,30,60,100,150)
for(d in depths){
  pretime <- Sys.time()
  ## Set up gNATSGO training points
  pts.gNAT.d <- pts.gNAT[,c("x","y",gsub(".tif","",cov.grids.names),colnames(pts.gNAT)[grepl(paste0("_",d,"_"),colnames(pts.gNAT))])]
  pts.gNAT.d$prop <- as.numeric(pts.gNAT.d[,paste0(prop,"_mean_",datastretchlab,"_",d,"_cm_gNATSGO")])/datastretch
  pts.gNAT.d[,paste0(prop,"_mean_",datastretchlab,"_",d,"_cm_gNATSGO")] <- NULL
  pts.gNAT.d <- na.omit(pts.gNAT.d)
  pts.gNAT.d$mtchtype <- "gNAT"
  pts.gNAT.d$peiid <- paste0("gNAT_",1:nrow(pts.gNAT.d))
  pts.gNAT.d$tid <- "gNAT"
  pts.gNAT.d$geo_wt <- 9
  pts.gNAT.d$locid <- paste(pts.gNAT.d$x,pts.gNAT.d$y,sep="_")
  pts.gNAT.d$x_std <- pts.gNAT.d$x
  pts.gNAT.d$y_std <- pts.gNAT.d$y
  pts.gNAT.d[,c("x","y")] <- NULL
  ## Now NASIS pts
  pts.extc <- subset(pts.ext.hor, as.numeric(pts.ext.hor$hzdept_r) <= d & as.numeric(pts.ext.hor$hzdepb_r) > d) # subset to chosen depth
  ptspred.list <- gsub(".tif","", cov.grids.names)# Take .tif off of the grid list to just get the variable names
  ptspred.list <- c(ptspred.list,"prop","mtchtype","peiid","tid","x_std","y_std","geo_wt","locid") #Add dependent variable
  pts.extcc <- pts.extc[,c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use
  pts.extcc <- na.omit(pts.extcc)# Remove any record with NA's (in any column - be careful)
  ### Check for duplication with lab pedons
  scd.pts.d <- subset(scd.pts.ext.hor, as.numeric(scd.pts.ext.hor$hzn_top) <= d & as.numeric(scd.pts.ext.hor$hzn_bot) > d)
  scd.pts.d$mtchtype <- "scd"
  scd.pts.d$tid <- "scd"
  scd.pts.d$lat <- scd.pts.d$y_std
  scd.pts.d$long <- scd.pts.d$x_std
  ## Check for duplicated pedons between NASIS and SCD using buffers
  coordinates(scd.pts.d) <- ~ long + lat
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
  pts.extcc <- rbind(pts.extcc,scd.pts.d,pts.gNAT.d)
  coordinates(pts.extcc) <- ~ x_std + y_std
  projection(pts.extcc) <- cov.proj
  ## Geocode classes
  pts.extcc@data <- dplyr::mutate(pts.extcc@data, geo_cls = ifelse(geo_wt<6,"gps",NA),
                                  geo_cls = ifelse(geo_wt==6,"gps2",geo_cls),
                                  geo_cls = ifelse(geo_wt==7,"gps3",geo_cls),
                                  geo_cls = ifelse(geo_wt==8,"unk",geo_cls),
                                  geo_cls = ifelse(geo_wt==9,"gNAT",geo_cls))
  #pts.extcc@data <- left_join(pts.extcc@data, pts_qual_cls[,c("peiid","source")], by="peiid")
  ## Apply transformation
  if(trans=="log10") {pts.extcc$prop_t <- log10(pts.extcc$prop + 0.1)}
  if(trans=="log") {pts.extcc$prop_t <- log(pts.extcc$prop + 1)}
  if(trans=="sqrt") {pts.extcc$prop_t <- sqrt(pts.extcc$prop)}
  if(trans=="none") {pts.extcc$prop_t <- pts.extcc$prop}
  ## Save pts available for modeling file
  saveRDS(pts.extcc, paste(predfolder,"/AvailPTS_", prop, '_',d, "_cm.rds",sep=""))
  #pts.extcc <- readRDS(paste(predfolder,"/AvailPTS_", prop, '_',d, "_cm.rds",sep=""))
  ## Set up QRF formula
  formulaStringRF <- as.formula(paste('prop_t ~', paste(gsub(".tif","", cov.grids.names), collapse="+")))

  ############ Cross Validate and examine metrics among Lab and Nasis pedons ##########
  ## Set training parameters for trials
  trn.params <- list(ntrees = 100, min.node.size = 1)


  ## Now prepare a data tier grid search using gRPI to rank
  geo_vec <- c("gps","gps_gps2","gps_gps2_gps3", "gps_gps2_gps3_unk")
  srce_vec <- c("scd","scd_direct","scd_direct_home","scd_direct_home_adjacent")
  grid_vec <- expand.grid(geo=geo_vec, srce=srce_vec, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
  #TODO: pull out train sets with only scd gps and gps2
  grid_vec <- grid_vec[!(grid_vec$srce=="scd"&grid_vec$geo=="gps"),]
  grid_vec <- grid_vec[!(grid_vec$srce=="scd"&grid_vec$geo=="gps_gps2"),]

  ## Call function to estimate gRPI using each data subset as training data
  Sys.time()
  ptsforgRPI.d <- pts.extcc@data[pts.extcc@data$mtchtype != "gNAT",]
  data_grid_df <- DSMprops::gRPI_estim_ranger(x = ptsforgRPI.d, gsamp = pts.gRPI, fm = formulaStringRF, griddf = grid_vec, os="linux",
                                              train.params = trn.params, nthreads = cpus)
  Sys.time()

  ## Now pick model mode ranking on gRPI and Rsq
  data_grid_df$gRPIfinal <-  (data_grid_df$gRPI.ave + data_grid_df$gRPI.med) / 2
  data_grid_df$gRPI.rank <- rank(data_grid_df$gRPIfinal,ties.method = "random")
  data_grid_df$Rsq.rank <- (nrow(data_grid_df)+1) - rank(data_grid_df$Rsq,ties.method = "average")
  data_grid_df$rank_ave <- apply(data_grid_df[,c("gRPI.rank","Rsq.rank")], MARGIN = 1, FUN=mean)
  data_grid_df$rank_final <- rank(data_grid_df$rank_ave,ties.method = "random")
  ## Determine lowest gRPI for spatial cross validation matching later
  gRPI_best <- data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),]$gRPIfinal
  ## Save RPI table
  write.table(data_grid_df, paste(predfolder,"/gRPIs_", prop,"_", d, "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
  # data_grid_df <- read.delim(paste(predfolder,"/gRPIs_", prop,"_", d, "_cm.txt",sep=""),stringsAsFactors = F)
  ## Select points from overall dataset based on combined rankings
  model_params <- str_split(data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),c("datagrid")][1],"_")[[1]]
  quant_l <- data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),c("quant_l")]
  quant_h <- data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),c("quant_h")]
  quants_vec <- c(quant_l,quant_h)
  srce_params <- model_params[(match("srce",model_params)+1):length(model_params)]
  geo_params <- model_params[(match("geo",model_params)+1):(match("srce",model_params)-1)]
  ## Subset training points to new optimized dataset
  pts.pcv <- pts.extcc[pts.extcc$mtchtype %in% srce_params | pts.extcc$mtchtype == "gNAT",]
  pts.pcv <- pts.pcv[(pts.pcv$mtchtype=="scd"&pts.pcv$geo_wt<8)|(pts.pcv$geo_cls %in% geo_params) | pts.pcv$mtchtype == "gNAT",] # subset, but leave scd > 2000

  ######### Match best CV scheme using gRPI.ave
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

  ## Combine CV tables and save in list as R object
  cv.lst <- list(cv10f,s1cv10f,s10cv10f,s50cv10f, s100cv10f)
  # saveRDS(cv.lst,paste(predfolder,"/CVlist_", prop, '_',d, "_cm.rds",sep="")) # takes forever...
  # #cv.lst <- readRDS(paste(predfolder,"/CVlist_", prop, '_',d, "_cm.rds",sep=""))

  ## Validation metrics for CVs at different spatial supports
  valmets_sCV <- DSMprops::valmetrics(xlst = cv.lst, trans = trans, quants = quants_vec, prop = prop, depth = d)
  valmets_sCV$RPIall <- (valmets_sCV$RPI.cvave + valmets_sCV$RPI.cvmed) / 2
  idx_val <- which(abs(valmets_sCV$RPIall-gRPI_best)==min(abs(valmets_sCV$RPIall-gRPI_best)))
  bestval <- valmets_sCV[idx_val,c("valtype")]
  bestvalpts <- get(bestval)
  valmets_sCV$bestval <- bestval
  ## Pick best CV resolution based on RPI match
  resol_rpi <- as.numeric(gsub("s", "", str_split(bestval,"cv")[[1]][1]))
  ## Save results
  write.table(valmets_sCV, paste(predfolder,"/sCVstats_", prop,"_", d, "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
  #saveRDS(bestvalpts,paste(predfolder,"/sCVfull_best_pts_", prop, '_',d, "_cm.rds",sep=""))
  #valmets_sCV <- read.delim(paste(predfolder,"/sCVstats_", prop,"_", d, "_cm.txt",sep=""),stringsAsFactors = F)

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

  ## Save training points file
  saveRDS(pts.pcv, paste(predfolder,"/TrainPTS_", prop, '_',d, "_cm.rds",sep=""))
  #pts.pcv <- readRDS(paste(predfolder,"/TrainPTS_", prop, '_',d, "_cm.rds",sep=""))

  ###### Best Spatial Cross validation plots
  ## All data
  viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
  scaleFUN <- function(x) round(x,0)
  gplt.dcm.2D.sCV <- ggplot(data=bestvalpts, aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1) + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("Measured") + ylab("Spat. CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste(bestval,"Spat. Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV
  ggsave(paste(predfolder,'/sValPlot_1to1_all_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.sCV, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  ## Now just SCD pedons
  gplt.dcm.2D.sCV.SCD <- ggplot(data=bestvalpts[bestvalpts$mtchtype=="scd",], aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("SCD Measured") + ylab("Spat. CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste(bestval,"Spat. Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV.SCD
  ggsave(paste(predfolder,'/sValPlot_1to1_scd_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.sCV.SCD, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  ## Now just with SCD pedons with GPS or after 2010
  gplt.dcm.2D.sCV.gSCD <- ggplot(data=bestvalpts[bestvalpts$mtchtype=="scd"&(bestvalpts$geo_cls=="gps"|bestvalpts$geo_cls=="gps2"|bestvalpts$geo_cls=="gps3"),], aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("gps SCD Measured") + ylab("Spat. CV Prediction") + scale_fill_gradientn( name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste(bestval,"Spat. Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV.gSCD
  ggsave(paste(predfolder,'/sValPlot_1to1_gscd_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.sCV.gSCD, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))

  ###### Cross validation plots
  ## All data
  viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
  scaleFUN <- function(x) round(x,0)
  gplt.dcm.2D.CV <- ggplot(data=cv.lst[[1]], aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1) + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste("Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV
  ggsave(paste(predfolder,'/ValPlot_1to1_all_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  ## Now just SCD pedons
  gplt.dcm.2D.CV.SCD <- ggplot(data=cv.lst[[1]][cv.lst[[1]]$mtchtype=="scd",], aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("SCD Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste("Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV.SCD
  ggsave(paste(predfolder,'/ValPlot_1to1_scd_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV.SCD, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  ## Now just with SCD pedons with GPS or after 2010
  gplt.dcm.2D.CV.gSCD <- ggplot(data=cv.lst[[1]][cv.lst[[1]]$mtchtype=="scd"&(cv.lst[[1]]$geo_cls=="gps"|cv.lst[[1]]$geo_cls=="gps2"|cv.lst[[1]]$geo_cls=="gps3"),], aes(prop_t, pcvpred)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("gps SCD Measured") + ylab("CV Prediction") + scale_fill_gradientn( name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste("Cross val", prop, d, "cm",sep=" "))
  #gplt.dcm.2D.CV.gSCD
  ggsave(paste(predfolder,'/ValPlot_1to1_gscd_',prop,'_',d,'_cm.tif',sep=""), plot = gplt.dcm.2D.CV.gSCD, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))

  ############### Build quantile Random Forest
  ## Determine 95% interquartile range for relative prediction interval
  varrange <- as.numeric(quantile(pts.pcv@data$prop, probs=c(0.975), na.rm=T)-quantile(pts.pcv@data$prop, probs=c(0.025),na.rm=T)) ## TRANSFORM IF NEEDED!
  if(varrange == 0) {varrange <- as.numeric(quantile(pts.pcv@data$prop, probs=c(0.995), na.rm=T)-quantile(pts.pcv@data$prop, probs=c(0.005),na.rm=T))}
  if(varrange == 0) {varrange <- as.numeric(quantile(pts.pcv@data$prop, probs=c(0.999), na.rm=T)-quantile(pts.pcv@data$prop, probs=c(0.001),na.rm=T))}
  if(varrange == 0) {varrange <- as.numeric(quantile(pts.pcv@data$prop, probs=c(0.9995), na.rm=T)-quantile(pts.pcv@data$prop, probs=c(0.0005),na.rm=T))}
  if(varrange == 0) {varrange <- as.numeric(quantile(pts.pcv@data$prop, probs=c(0.9999), na.rm=T)-quantile(pts.pcv@data$prop, probs=c(0.0001),na.rm=T))}
  if(varrange == 0) {varrange <- as.numeric(max(pts.pcv@data$prop, na.rm=T)-min(pts.pcv@data$prop,na.rm=T))}
  ## Train global ranger model
  rf.qrf <- ranger(formulaStringRF, data=pts.pcv@data, num.trees = trn.params$ntrees, quantreg = T, num.threads = 124,
                   min.node.size = trn.params$min.node.size,
                   importance = "permutation") #case.weights = pts.pcv$tot_wts,
  ## OOB error
  # rf.qrf
  ## Peak at importance
  # imp <- data.frame(var=names(importance(rf.qrf)),imp_meas = unname(importance(rf.qrf)))
  # imp[order(imp$imp_meas, decreasing = T),][1:30,]
  ## Linear Adjustment for bias in low and high predictions
  pts.pcv@data$trainpreds <- predict(rf.qrf, data=pts.pcv@data, num.threads = cpus)$predictions
  attach(pts.pcv@data)
  rf_lm_adj <- lm(prop_t ~ trainpreds)
  detach(pts.pcv@data)
  pts.pcv@data$trainpredsadj <- predict(rf_lm_adj, newdata=pts.pcv@data)
  ## Save important documentation and intermediate files as R objects
  saveRDS(rf.qrf, paste(predfolder,"/rangerQRF_", prop, '_',d, "_cm.rds",sep=""))
  saveRDS(rf_lm_adj, paste(predfolder,"/rflmadj_RFmodel_",prop,"_", d, "_cm.rds",sep=""))
  ## Block to open prior files for updating work
  #rf.qrf <- readRDS(paste(predfolder,"/rangerQRF_", prop, '_',d, "_cm.rds",sep=""))
  #rf_lm_adj <- readRDS(paste(predfolder,"/rflmadj_RFmodel_",prop,"_", d, "_cm.rds",sep=""))


  ############################## Raster Preditions ######################################################
  ##### Reference covar rasters to use in prediction
  rasters <- stack(cov.grids)
  # b <- brick(rasters)
  #rasters <- readAll(rasters) # Brings into active memory
  #names(rasters)

  ## Ranger Predict
  predfun <- function(model, ...) predict(model, ...)$predictions

  ## Predict onto covariate grid
  ## Parallelized predict
  rasterOptions(maxmemory = 6.5e+09,chunksize = 6e+08)# maxmemory = 6e+09,chunksize = 6e+08 40-50% cpu, 530GB
  beginCluster(cpus,type='SOCK')
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
      smrest <- mean(10^(pts.pcv@data$prop_t - pts.pcv@data$trainpredsadj))
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
      smrest <- mean(exp(pts.pcv@data$prop_t - pts.pcv@data$trainpredsadj))
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
      smrest <- mean((pts.pcv@data$prop_t - pts.pcv@data$trainpredsadj)^2)
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

