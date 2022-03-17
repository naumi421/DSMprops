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

## Key Folder Locations
predfolder <- "/home/tnaum/data/Hyb100m/Histosols"
repofolder <- "/home/tnaum/data/repos/DSMprops"
covfolder <- "/media/nped/SG100_covars"
ptsfolder <- "/media/sped/Hyb100m_gdrv/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final"

# ## Sourced functions
# source(paste0(repofolder,"/exec/covar_dev/rfe_rangerFuncs.R"))

######## Load soil profile collection ##############
## Geographic coordinate quality levels
pts_geocode <- readRDS("/media/sped/Hyb100m_gdrv/2020_Pedons/geocode_weighting.RDS") # From Dave White
# pts_geocode_wt <- readRDS("/home/tnaum/data/Hyb100m_gdrv/2020_Pedons/geocode_weighting_v2.RDS") # From Dave White
pts_geocode$geo_wt <- pts_geocode$wt
pts_geocode$wt <- NULL
pts_geocode$peiid <- as.character(pts_geocode$peiid)
# ## Pedon quality clases and weights
# pts_qual_cls <- readRDS("/home/tnaum/data/Hyb100m_gdrv/2020_Pedons/pedonquality_weighting.RDS") # Dave White
# pts_qual_cls <- pts_qual_cls[!duplicated(pts_qual_cls$peiid),]
## Pedons with extracted SSUGO component data
pts <- readRDS(paste(ptsfolder,"/NASIS_all_component_horizon_match_SPC_ssurgo20.rds",sep=""))
pts.proj <- proj4string(pts)
n.pts <- pts@site
## Bring in orig nasis points to eliminate those with partial coords
load("/media/sped/Hyb100m_gdrv/2020_Pedons/nasis_sites_20210325.RData") # object s
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
polybound <- readOGR("/media/sped/GIS_data/US_Census_500k/cb_2020_us_state_500k", "conus_bound")
polybound <- spTransform(polybound, cov.proj)
n.pts <- n.pts[polybound,]

####### Now create a random sample of points in CONUS to use for gRPI estimation
# pts.gRPI <- spsample(polybound[1,], 1000000, type = 'random')
# pts.gRPI <- SpatialPointsDataFrame(pts.gRPI, data.frame(row.names=row.names(pts.gRPI), ID=1:length(pts.gRPI)))

## Plot to ensure alignment bw points and rasters
# plot(projgrid)
# plot(n.pts, add=TRUE)
# plot(pts.gRPI, add=TRUE)

## Parallelized extract for gRPI points
# rasterOptions(maxmemory = 4e+09)
# pts.gRPI <- DSMprops::parPTextr(sp = pts.gRPI, gridlist = cov.grids, os = "linux",nthreads = 50)
# pts.gRPI <- na.omit(pts.gRPI)
# ## Save points
# saveRDS(pts.gRPI, paste(ptsfolder,"/CONUS_random_gRPIsamp_spatcovs.rds",sep=""))
## Updated extract for CONUS
pts.gRPI <- readRDS(paste(ptsfolder,"/CONUS_random_gRPIsamp_spatcovs.rds",sep=""))


## Parallelized extract for nasis points: (larger datasets)
# rasterOptions(maxmemory = 4e+09)
# pts.ext <- DSMprops::parPTextr(sp = n.pts, gridlist = cov.grids, os = "linux",nthreads = 50)
# ## Save points
# saveRDS(pts.ext, paste(predfolder,"/CONUS_nasis_extracted_spatcovs.rds",sep=""))
## Updated extract for CONUS
pts.ext <- readRDS(paste(predfolder,"/CONUS_nasis_extracted_spatcovs.rds",sep=""))
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

## Trim down dataset to more trusted data tiers and add taxonomy fields
pts.ext.hor <- pts.ext[pts.ext$mtchtype != "adjacent" & pts.ext$geo_wt < 8,]
pts.ext.hor <- left_join(pts.ext.hor, s[,c("peiid","taxorder","taxsuborder","taxgrtgroup","taxsubgrp","taxpartsize")], by = "peiid")
pts.ext.hor <- pts.ext.hor[!is.na(pts.ext.hor$taxsuborder),]

## Prep nasis training data for Random Forest
pts.ext.hor$prop <- ifelse(pts.ext.hor$taxorder == "histosols" & pts.ext.hor$taxsuborder != "folists", "wet_histsols", "none")  ## UPDATE EVERY TIME
prop <- "wet_histsols" ## Dependent variable
summary.factor(pts.ext.hor$prop)


##### Load and prep SCD data
## RSQlite workflow form https://github.com/ncss-tech/gsp-sas/blob/master/lab_data.Rmd
con <- dbConnect(RSQLite::SQLite(), "/media/sped/Hyb100m_gdrv/2020_Pedons/KSSL-snapshot-draft/KSSL-data.sqlite")
(ldm_names <- dbListTables(con))
ldm <- lapply(c("NCSS_Layer","NCSS_Site_Location","pH_and_Carbonates","NCSS_Pedon_Taxonomy"), function(x) dbReadTable(con , x))
names(ldm) <- c("NCSS_Layer","NCSS_Site_Location","pH_and_Carbonates","NCSS_Pedon_Taxonomy")
dbDisconnect(con)

## Now wrangle tables
ldm$NCSS_Pedon_Taxonomy$peiid <- ldm$NCSS_Pedon_Taxonomy$pedoniid
scd.pts <- ldm$NCSS_Site_Location
scd.pts <- left_join(scd.pts,ldm$NCSS_Pedon_Taxonomy[, c("site_key","peiid","SSL_taxorder","SSL_taxsuborder",
                                                         "corr_taxsuborder","corr_taxorder","samp_taxorder","samp_taxsuborder")], by="site_key")
scd.pts <- scd.pts[!is.na(scd.pts$corr_taxorder),]

# # ### SCD prep: Weed out points with imprecise coordinates ###
# scd.pts$latnchar <- nchar(abs(scd.pts$latitude_decimal))
# scd.pts$longnchar <- nchar(abs(scd.pts$longitude_decima))
# scd.pts <- subset(scd.pts, scd.pts$latnchar > 5 & scd.pts$longnchar > 6)
# ## Location ID for later use and remove duplicates
# scd.pts$locid <- paste(scd.pts$latitude_decimal,scd.pts$longitude_decima,sep="_")
# scd.pts <- scd.pts[!duplicated(scd.pts$locid),] # Over 20k duplicates...
# scd.pts <- scd.pts[!duplicated(scd.pts$peiid),]
#
# ### Turn into spatial file
# coordinates(scd.pts) <- ~ longitude_decima + latitude_decimal # Typos in sqlite KSSL snapshot...
# temp.proj <- CRS("+proj=longlat +datum=WGS84") ## specify projection
# projection(scd.pts) <- temp.proj
# # ######## Clip with boundary if necessary ###########
# scd.pts <- spTransform(scd.pts, cov.proj)
# scd.pts <- scd.pts[polybound,]

## Further SCD prep
## Extract covariates for prediction onto SCD points
# scd.pts.ext <- DSMprops::parPTextr(scd.pts, cov.grids, os = "linux", nthreads=50)
# ## Save scd pts with extraction
# saveRDS(scd.pts.ext, paste(predfolder,"/SCD","_extracted_spatcovs.rds",sep=""))
scd.pts.ext <- readRDS(paste(predfolder,"/SCD","_extracted_spatcovs.rds",sep=""))
## Create geo-coordinate quality weights
scd.pts.ext$date <- as.Date(scd.pts.ext$site_obsdate, format = "%m/%d/%Y")
scd.pts.ext$decdate <- lubridate::decimal_date(scd.pts.ext$date)
scd.pts.ext$peiid <- as.character(scd.pts.ext$peiid)
scd.pts.ext <- left_join(scd.pts.ext, pts_geocode, by="peiid")
scd.pts.ext$geo_wt <- ifelse(is.na(scd.pts.ext$geo_wt) & scd.pts.ext$date >= "2010-01-01", 6, scd.pts.ext$geo_wt)
scd.pts.ext$geo_wt <- ifelse(is.na(scd.pts.ext$geo_wt) & scd.pts.ext$date >= "2000-01-01", 7, scd.pts.ext$geo_wt)
scd.pts.ext$geo_wt <- ifelse(is.na(scd.pts.ext$geo_wt), 8, scd.pts.ext$geo_wt) ## Awful lot of class 8s...
scd.pts.ext$geo_wt <- ifelse(scd.pts.ext$geo_wt == 8 & scd.pts.ext$date >= "2000-01-01", 7, scd.pts.ext$geo_wt) #updated 7/8/21

## Now join scd.hor to scd.pts and get rid of any duplicates
scd.pts.ext.hor <- inner_join(scd.pts.ext, scd.pts[,c("site_key","peiid","corr_taxsuborder","corr_taxorder")], by="site_key")
scd.pts.ext.hor <- scd.pts.ext.hor[!duplicated(scd.pts.ext.hor$locid),]

## SCD prep for RF
scd.pts.ext.hor$prop <- ifelse(scd.pts.ext.hor$corr_taxorder == "histosols" & scd.pts.ext.hor$corr_taxsuborder != "folists", "wet_histsols", "none")  ## UPDATE EVERY TIME  ## UPDATE everytime!
scdprop <- "wet_histsols"
scd.pts.ext.hor$tid <- "scd"
summary.factor(scd.pts.ext.hor$prop)

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
ft_wts_ref <- readRDS("/media/sped/Hyb100m_gdrv/2020_Pedons/ref_df.RDS") # From Stephen Roecker, 5/3/2021

##### Blocks to train and predict depth (adapted from depthwise predictions for other properties...)
d <- "all"

  pretime <- Sys.time()
  pts.extc <- pts.ext.hor
  ptspred.list <- gsub(".tif","", cov.grids.names)# Take .tif off of the grid list to just get the variable names
  ptspred.list <- c(ptspred.list,"prop","mtchtype","peiid","tid","x_std","y_std","geo_wt","locid") #Add dependent variable
  pts.extcc <- pts.extc[,c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use
  pts.extcc <- na.omit(pts.extcc)# Remove any record with NA's (in any column - be careful)
  ### Check for duplication with lab pedons
  scd.pts.d <- scd.pts.ext.hor
  ## Remove SCD dups
  scd.pts.d <- scd.pts.d[!duplicated(scd.pts.d$locid),]
  scd.pts.d$x_std <- scd.pts.d$longitude_decima
  scd.pts.d$y_std <- scd.pts.d$latitude_decimal
  scd.pts.d$mtchtype <- "scd"
  scd.pts.d$tid <- "scd"
  ## Check for duplicated pedons between NASIS and SCD using buffers
  coordinates(scd.pts.d) <- ~ longitude_decima + latitude_decimal
  projection(scd.pts.d) <- cov.proj
  pts.buf <- n.pts[n.pts$peiid %in% pts.extcc$peiid,]
  scd.pts.d.bufc <- raster::buffer(scd.pts.d,width=5,dissolve=F)
  pts.buf.exlude <- pts.buf[scd.pts.d.bufc,]
  pts.extcc <- pts.extcc[!pts.extcc$peiid %in% pts.buf.exlude@data$peiid,]
  pts.buf <- pts.buf[!pts.buf$peiid %in% pts.buf.exlude@data$peiid,c("peiid")]
  ## Now combine cleaned up datasets
  scd.pts.d <- scd.pts.d@data
  scd.pts.d$peiid <- scd.pts.d$peiid.x
  scd.pts.d <- scd.pts.d[,c(ptspred.list)]
  scd.pts.d <- na.omit(scd.pts.d)
  pts.extcc <- rbind(pts.extcc,scd.pts.d)
  coordinates(pts.extcc) <- ~ x_std + y_std
  projection(pts.extcc) <- cov.proj
  ## Exclude pre-2000 points due to uncertain coordinates
  pts.extcc <- pts.extcc[pts.extcc$geo_wt < 8,]
  ## Save pts available for modeling file
  saveRDS(pts.extcc, paste(predfolder,"/AvailPTS_", prop, '_',d, "_cm.rds",sep=""))
  #pts.extcc <- readRDS(paste(predfolder,"/AvailPTS_", prop, '_',d, "_cm.rds",sep=""))
  ## Set up QRF formula and matrices
  formulaStringRF <- as.formula(paste('prop ~', paste(gsub(".tif","", cov.grids.names), collapse="+")))

  ############ Cross Validate and examine metrics among Lab and Nasis pedons ##########
  ## Set training parameters for trials
  trn.params <- list(ntrees = 100, min.node.size = 1)


  # ## Now prepare a data tier grid search using gRPI to rank
  # geo_vec <- c("gps","gps_gps2","gps_gps2_gps3", "gps_gps2_gps3_unk")
  # srce_vec <- c("scd","scd_direct","scd_direct_home","scd_direct_home_psuedo","scd_direct_home_psuedo_adjacent")
  # grid_vec <- expand.grid(geo=geo_vec, srce=srce_vec, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
  # ## pull out train sets with only scd gps and gps2
  # grid_vec <- grid_vec[!(grid_vec$srce=="scd"&grid_vec$geo=="gps"),]
  # grid_vec <- grid_vec[!(grid_vec$srce=="scd"&grid_vec$geo=="gps_gps2"),]
  #
  # ## Call function to estimate gRPI using each data subset as training data
  # Sys.time()
  # data_grid_df <- DSMprops::gRPI_estim_ranger(x = pts.extcc@data, gsamp = pts.gRPI, fm = formulaStringRF, griddf = grid_vec, os="linux",
  #                                              train.params = trn.params, nthreads = 30)
  # Sys.time()
  #
  # ## Now pick model mode ranking on gRPI and Rsq
  # data_grid_df$gRPIfinal <-  (data_grid_df$gRPI.ave + data_grid_df$gRPI.med) / 2
  # data_grid_df$gRPI.rank <- rank(data_grid_df$gRPIfinal,ties.method = "random")
  # data_grid_df$Rsq.rank <- (nrow(data_grid_df)+1) - rank(data_grid_df$Rsq,ties.method = "average")
  # data_grid_df$rank_ave <- apply(data_grid_df[,c("gRPI.rank","Rsq.rank")], MARGIN = 1, FUN=mean)
  # data_grid_df$rank_final <- rank(data_grid_df$rank_ave,ties.method = "random")
  # ## Determine lowest gRPI for spatial cross validation matching later
  # gRPI_best <- data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),]$gRPIfinal
  # ## Save RPI table
  # write.table(data_grid_df, paste(predfolder,"/gRPIs_", prop,"_", d, "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
  # data_grid_df <- read.delim(paste(predfolder,"/gRPIs_", prop,"_", d, "_cm.txt",sep=""),stringsAsFactors = F)
  # ## Select points from overall dataset based on combined rankings
  # model_params <- str_split(data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),c("datagrid")][1],"_")[[1]]
  # quant_l <- data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),c("quant_l")]
  # quant_h <- data_grid_df[data_grid_df$rank_final==min(data_grid_df$rank_final),c("quant_h")]
  # quants_vec <- c(quant_l,quant_h)
  # srce_params <- model_params[(match("srce",model_params)+1):length(model_params)]
  # geo_params <- model_params[(match("geo",model_params)+1):(match("srce",model_params)-1)]
  # ## Subset training points to new optimized dataset
  # pts.pcv <- pts.extcc[pts.extcc$mtchtype %in% srce_params,]
  # pts.pcv <- pts.pcv[(pts.pcv$mtchtype=="scd"&pts.pcv$geo_wt<8)|(pts.pcv$geo_cls %in% geo_params),] # subset, but leave scd > 2000
  #
  #
  # ########### Recursive feature elimination
  # #The simulation will fit models a set number of covariate subsets based on an initial full model variable importance
  # rf.rfe.init <- ranger(formulaStringRF, data=pts.pcv@data, num.trees = trn.params$ntrees, num.threads = 30,
  #                  min.node.size = trn.params$min.node.size, importance = "permutation")
  # vars.imp <- data.frame(import = unname(rf.rfe.init$variable.importance), var = names(rf.rfe.init$variable.importance))
  # vars.imp <- vars.imp[order(vars.imp$import),]
  # num.subsets <- 16
  # subset.increm <- floor(nrow(vars.imp)/num.subsets)
  # subsets <- as.vector(subset.increm)
  # rfe_tab <- data.frame(var_idx = nrow(vars.imp), rsq = rf.rfe.init$r.squared)
  # for(s in 1:(num.subsets-2)){subsets <- append(subsets,subsets[length(subsets)]+subset.increm)}
  # subsets <- subset(subsets, subsets > 20) # Decision to include at least 30 vars so one odd var doesn't dominate and create weird patterns
  # set.seed(12)
  # ranger_rfe_fn <- function(ss){
  #   subsvars <- vars.imp$var[(nrow(vars.imp)-ss):nrow(vars.imp)]
  #   formulaStringRF_subset <- as.formula(paste('prop_t ~', paste(subsvars, collapse="+")))
  #   rf.rfe.subs <- ranger(formulaStringRF_subset, data=pts.pcv@data, num.trees = trn.params$ntrees, num.threads = detectCores() - 1,
  #                         min.node.size = trn.params$min.node.size, importance = "none")
  #   rfe_tab.sub <- data.frame(var_idx = ss, rsq = rf.rfe.subs$r.squared)
  #   return(rfe_tab.sub)
  # }
  # ## List apply of rfe function
  # Sys.time()
  # rfe_tab.subs <- lapply(subsets,try(ranger_rfe_fn))
  # Sys.time()
  # rfe_tab.subsdf <- plyr::rbind.fill(rfe_tab.subs)
  # rfe_tab <- rbind(rfe_tab,rfe_tab.subsdf)
  # rfe_tab <- rfe_tab[order(rfe_tab$var_idx),]
  # ## Save rfe results
  # saveRDS(rfe_tab, paste(predfolder,"/rf.RFE_", prop,"_", d, "_cm",".rds",sep=""))
  # rfe_tab <- readRDS(paste(predfolder,"/rf.RFE_", prop,"_", d, "_cm",".rds",sep=""))
  # #plot(rfe_tab$rsq ~ rfe_tab$var_idx) ##
  # # 4.4 See list of predictors
  # ## Create idex to grab simpler models close to max Rsq
  # idx_rfe <- which(rfe_tab$rsq > 0.98*max(rfe_tab$rsq))
  # idx_rfe <- rfe_tab[min(idx_rfe),]$var_idx
  # ## Prep for Random Forest
  # rfe.list <- vars.imp$var[(nrow(vars.imp)-idx_rfe):nrow(vars.imp)] # Select variables from rfe
  # # varlist_pruned <- names(sort(rf.RFE$fit$importance[,1], decreasing=T)[0:40]) # chose top n variables
  # formulaStringRF_RFE <- as.formula(paste('prop_t ~ ',paste(rfe.list, collapse="+")))# put in dep variable name
  #
  # ## TODO!!!!!!!!!!!!!!!! Update formula handling in CV and sCV functions
  #
  # ######### Match best CV scheme using gRPI.ave
  # ## Normal 10-fold cross validation
  # cv10f <- DSMprops::CVranger(x = pts.pcv@data, fm = formulaStringRF_RFE, train.params = trn.params,quants=quants_vec,
  #                             nfolds = 10, nthreads = 30, os = "linux") # 5min
  # ## Spatial 10-fold cross validation on 1km blocks
  # s1cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF_RFE, train.params = trn.params,quants=quants_vec,
  #                                   rast = img10kf, nfolds = 10, nthreads = 30, resol = 1, os = "linux") # 6 min
  # ## Spatial 10-fold cross validation on 10km blocks
  # s10cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF_RFE, train.params = trn.params,quants=quants_vec,
  #                                    rast = img10kf, nfolds = 10, nthreads = 30, resol = 10, os = "linux") # 6min
  # ## Spatial 10-fold cross validation on 50km blocks
  # s50cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF_RFE, train.params = trn.params,quants=quants_vec,
  #                                    rast = img10kf, nfolds = 10, nthreads = 30, resol = 50, os = "linux") # 6min
  # ## Spatial 10-fold cross validation on 100km blocks
  # s100cv10f <- DSMprops::SpatCVranger(sp = pts.pcv, fm = formulaStringRF_RFE, train.params = trn.params,quants=quants_vec,
  #                                     rast = img10kf, nfolds = 10, nthreads = 30, resol = 100, os = "linux")
  #
  # ## Combine CV tables and save in list as R object
  # cv.lst <- list(cv10f,s1cv10f,s10cv10f,s50cv10f, s100cv10f)
  # saveRDS(cv.lst,paste(predfolder,"/CVlist_", prop, '_',d, "_cm.rds",sep="")) # takes forever...
  # #cv.lst <- readRDS(paste(predfolder,"/CVlist_", prop, '_',d, "_cm.rds",sep=""))
  #
  # ## Validation metrics for CVs at different spatial supports
  # valmets_sCV <- DSMprops::valmetrics(xlst = cv.lst, trans = trans, quants = quants_vec, prop = prop, depth = d)
  # valmets_sCV$RPIall <- (valmets_sCV$RPI.cvave + valmets_sCV$RPI.cvmed) / 2
  # idx_val <- which(abs(valmets_sCV$RPIall-gRPI_best)==min(abs(valmets_sCV$RPIall-gRPI_best)))
  # bestval <- valmets_sCV[idx_val,c("valtype")]
  # bestvalpts <- get(bestval)
  # valmets_sCV$bestval <- bestval
  # ## Pick best CV resolution based on RPI match
  # resol_rpi <- as.numeric(gsub("s", "", str_split(bestval,"cv")[[1]][1]))
  # ## Save results
  # write.table(valmets_sCV, paste(predfolder,"/sCVstats_", prop,"_", d, "_cm.txt",sep=""), sep = "\t", row.names = FALSE)
  # #saveRDS(bestvalpts,paste(predfolder,"/sCVfull_best_pts_", prop, '_',d, "_cm.rds",sep=""))
  # #valmets_sCV <- read.delim(paste(predfolder,"/sCVstats_", prop,"_", d, "_cm.txt",sep=""),stringsAsFactors = F)

  ## Characterize pt density for case weights
  polywin    <- as(polybound, "owin")
  pts.ppp    <- as(pts.extcc, "ppp")
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

  ## Save training points file
  saveRDS(pts.extcc, paste(predfolder,"/TrainPTS_", prop, '_',d, "_cm.rds",sep=""))
  #pts.pcv <- readRDS(paste(predfolder,"/TrainPTS_", prop, '_',d, "_cm.rds",sep=""))


  ## Train global ranger model
  pts.extcc@data$prop <- as.factor(pts.extcc@data$prop)
  rf.qrf <- ranger(formulaStringRF, data=pts.extcc@data, num.trees = trn.params$ntrees, num.threads = 30,
                   min.node.size = trn.params$min.node.size,
                   importance = "permutation")
  rf.qrf.prb <- ranger(formulaStringRF, data=pts.extcc@data, num.trees = trn.params$ntrees, num.threads = 30,
                   min.node.size = trn.params$min.node.size,
                   importance = "permutation",probability = TRUE)
  ## OOB error
  # rf.qrf
  ## Peak at importance
  # imp <- data.frame(var=names(importance(rf.qrf)),imp_meas = unname(importance(rf.qrf)))
  # imp[order(imp$imp_meas, decreasing = T),][1:30,]
  ## Save important documentation and intermediate files as R objects
  saveRDS(rf.qrf.prb, paste(predfolder,"/rangerQRFprob_", prop, '_',d, ".rds",sep=""))
  saveRDS(rf.qrf, paste(predfolder,"/rangerQRF_", prop, '_',d, ".rds",sep=""))
  ## Block to open prior files for updating work
  #rf.qrf <- readRDS(paste(predfolder,"/rangerQRF_", prop, '_',d, "_cm.rds",sep=""))


  ############################## Raster Preditions ######################################################
  ##### Reference covar rasters to use in prediction
  rasters <- stack(cov.grids)
  # b <- brick(rasters)
  #rasters <- readAll(rasters) # Brings into active memory
  #names(rasters)

  ## Ranger Predict
  predfun <- function(model, ...) predict(model, ...)$predictions
  # TODO figure out why internal ranger parallelization did not seem to work.
  # predtst <- predict(rasters, rf.qrf, fun=predfun, type="response", progress="text", num.threads = 50)


  ## Predict onto covariate grid
  ## Parallelized predict
  rasterOptions(maxmemory = 6e+09,chunksize = 8e+08)# maxmemory = 1e+09,chunksize = 1e+08 for soilmonster
  beginCluster(30,type='SOCK')
  Sys.time()
  predprob <- clusterR(rasters, predict, args=list(model=rf.qrf.prb, fun=predfun, index = 2),progress="text")
  Sys.time()
  #pred <- predict(rasters, rf.qrf, fun=predfun, progress="text")
  pred <- clusterR(rasters, predict, args=list(model=rf.qrf, fun=predfun), progress="text") # works ok
  #pred <- clusterR(rasters, predict, args=list(model=soiclass),progress="text")
  intfn <- function(x){
    ind <- x*100
    return(ind)}
  predprobint <- clusterR(predprob, calc, args=list(fun=intfn),progress='text')
  Sys.time()
  ## End cluster to clear cache and memory
  endCluster()

  ## Write rasters
  writeRaster(predprobint, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_prob_",d,".tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype="INT1U", progress="text")
  writeRaster(pred, overwrite=TRUE,filename=paste(predfolder,"/",prop,"_",d,".tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype="INT1U", progress="text")

