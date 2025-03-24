###############################################
## Script to create weighted area map unit soil property
## covariate layers from gNATSGO for use as DSM covariates
## T Nauman, 11/29/2022

### Packages
required.packages <- c("raster", "sp", "rgdal","doParallel", "plyr", "dplyr","sf","cluster","RSQLite","lubridate","rasterVis","maptools","classInt","mapview")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
##Raster settings
rasterOptions(maxmemory = 2e+09, chunksize = 2e+08)

### Key folders
projfldr <- "/mnt/covs/gNATSGO_CONUS_tif/covs100m_by_prop"
gnatsgofldr <- "/mnt/covs/gNATSGO_CONUS_tif"
gnatsgotabfldr <- "/mnt/covs/gNATSGO_CONUS_tif/gNATSGO_Tabular_CSV"
sg100fldr <- "/mnt/covs/covs100m"

### gNATSGO tif
sstif <- raster(paste0(gnatsgofldr,"/gNATSGO-mukey.tif"))
ss.crs <- raster::crs(sstif) ## Maintains both proj4 and proj6 properties

### SoilGrids 100m reference grids
sgtif <- raster(paste(sg100fldr,"/T07PRI5.tif",sep=""))
sg.crs <- raster::crs(sgtif)

## Snap gNATSGO grid to sq100
# rasterOptions(maxmemory = 5e+08,chunksize = 5e+07)
# cpus <- 116
# beginCluster(cpus,type='SOCK')
# sstif <- projectRaster(sstif,sgtif,method='ngb',datatype='INT4U')
# endCluster()
# writeRaster(sstif, overwrite=TRUE,filename=paste0(gnatsgofldr,"/mukeys_sg100grid_gNATSGO.tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT4U', progress="text")
rasterOptions(maxmemory = 6e+09, chunksize = 4e+08) # reset raster params
sstif <- raster(paste0(gnatsgofldr,"/mukeys_sg100grid_gNATSGO.tif"))

#### Get unique mukeys
mukeys_tif <- raster::unique(sstif)

## Mapunit table
mapunit.df <- read.csv(paste0(gnatsgotabfldr,"/mapunit.csv"))
mukeys <- mukeys_tif[mukeys_tif %in% mapunit.df$mukey]

## Check extent in mapview
# mapview(sstif)

##### Also load necessary data tables and subset to study area
## components
comp.df <- read.csv(paste0(gnatsgotabfldr,"/component.csv"))
comp.df$compname <- toupper(comp.df$compname)
horizon.df <- read.csv(paste0(gnatsgotabfldr,"/chorizon.csv"))
## Ecosites table
eco.df <- read.csv(paste0(gnatsgotabfldr,"/coecoclass.csv"))
## Restriction table
restr.df <- read.csv(paste0(gnatsgotabfldr,"/corestrictions.csv"))
## Surface Frags table
sfrag.df <- read.csv(paste0(gnatsgotabfldr,"/cosurffrags.csv"))
## Legend table
legend.df <- read.csv(paste0(gnatsgotabfldr,"/legend.csv"))
## subset tables down to area of interest
comps <- comp.df[comp.df$mukey %in% mukeys,]
cokeys <- comps$cokey
horizs <- horizon.df[horizon.df$cokey %in% cokeys,]
ecosites <- eco.df[eco.df$cokey %in% cokeys,]
sfrags <- sfrag.df[sfrag.df$cokey %in% cokeys,]
restrs <- restr.df[restr.df$cokey %in% cokeys,]
horiz_comps <- left_join(horizs,comps, by="cokey")
horiz_comps$chkey <- as.character(horiz_comps$chkey)
horiz_comps$cokey <- as.character(horiz_comps$cokey)
horiz_comps[] <- lapply(horiz_comps, function(x) if (is.factor(x)) as.character(x) else {x})
## Horizon Frag table
chfrags <- read.csv(paste0(gnatsgotabfldr,"/chfrags.csv"))
chfrags_bld <- chfrags[chfrags$fragsize_r > 599,]
boulders_chkey <- plyr::ddply(chfrags_bld,~chkey,summarise,fragvol_r = sum(fragvol_r))
boulders_chkey$boulders <- boulders_chkey$fragvol_r
boulders_chkey$fragvol_r <- NULL
saveRDS(boulders_chkey,paste0(projfldr,"boulders_chkey.rds"))
boulders_chkey <- readRDS(paste0(projfldr,"boulders_chkey.rds"))
horiz_comps$boulders <- ifelse(is.na(horiz_comps$boulders), 0, horiz_comps$boulders) # NAs set to zero as these are not recorded in chfrags table
horiz_comps <- subset(horiz_comps, horiz_comps$boulders < 101) ## Check for common errors, none for boulders

## Set gNATSGO property
prop <- 'boulders'
## Vector of depths to render
depths <- c(0,5,15,30,60,100,150)
## Scaling and datatype
data_type <- "INT1U" # from raster::dataType - INT1U, INT1S, INT2S, INT2U, INT4S, INT4U, FLT4S, FLT8S
datastretch <- 1
datastretchlab <- paste(datastretch,"x",sep="")

## Create folder for property subfiles
propdir <- paste0(projfldr,"/",prop)
dir.create(propdir)

## Prep for list apply parallellization
cpus <- min(detectCores()-2,124,length(mukeys))
mukeys.df <- data.frame(mukey=mukeys,grp = sample.int(cpus,size=length(mukeys),replace=T))
horiz_comps <- left_join(horiz_comps, mukeys.df, by = "mukey")
grps <- unique(horiz_comps$grp)
rasterOptions(maxmemory = 6e+09, chunksize = 8e+08)

## Set up depth loop for gNATSGO rasters
# d <- 60
for(d in depths){
  ## Prep property map to a depth interval
  propfn <- function(g){
    ud <- d - 0.5 # upper depth
    ld <- d + 0.5 # lower depth
    horiz_c <- horiz_comps[horiz_comps$grp == g,]
    horiz_c <- subset(horiz_c, (as.numeric(horiz_c$hzdept_r) <= ud & as.numeric(horiz_c$hzdepb_r) > ud) |
                        (as.numeric(horiz_c$hzdepb_r) >= ld & as.numeric(horiz_c$hzdept_r) < ld) ) # subset to horizon overlapping range
    horiz_c$indepthupper <- ifelse(as.numeric(horiz_c$hzdept_r) <= ud & as.numeric(horiz_c$hzdepb_r) < ld ,
                                   as.numeric(horiz_c$hzdepb_r) - ud, 0)
    horiz_c$indepthlower <- ifelse(as.numeric(horiz_c$hzdept_r) >= ud & as.numeric(horiz_c$hzdepb_r) >= ld,
                                   ld - as.numeric(horiz_c$hzdept_r),0)
    horiz_c$indepthoverlap <- ifelse(as.numeric(horiz_c$hzdept_r) < ud & as.numeric(horiz_c$hzdepb_r) > ld,
                                     ld - ud,0)
    horiz_c$indepth <- horiz_c$indepthupper + horiz_c$indepthlower + horiz_c$indepthoverlap
    horiz_c$depthwt <- horiz_c$indepth / (ld-ud)
    horiz_c <- horiz_c[,c(prop, "cokey","comppct_r","depthwt","mukey")]
    horiz_c <- na.omit(horiz_c)
    horiz_c_depthwt <- horiz_c[,c("cokey","depthwt")]
    horiz_c_depthwt <- ddply(horiz_c_depthwt,~cokey,summarise,depthwtsum=sum(depthwt)) # determine total depth weight to allow for post NA normalization
    horiz_c <- merge(horiz_c, horiz_c_depthwt,by="cokey")
    horiz_c$depthwtnorm <- horiz_c$depthwt / horiz_c$depthwtsum
    horiz_c$prop_depthwt <- horiz_c$depthwtnorm * horiz_c[,prop]
    comp_c <- ddply(horiz_c,~cokey,summarise,propc=sum(prop_depthwt),mukey=mean(as.numeric(mukey)),comppct_r=mean(comppct_r))

    ### Now summarize comps at mu level
    comp_c_cmpct <- comp_c[,c("mukey","comppct_r")]
    comp_c_cmpct <- ddply(comp_c_cmpct,~mukey,summarise,compsum=sum(comppct_r)) # determine total % named components in MU
    comp_c <- merge(comp_c, comp_c_cmpct,by="mukey")
    comp_cc <- subset(comp_c, comp_c$compsum > 10)
    comp_cc$propcwtd <- comp_cc$propc * (comp_cc$comppct_r / comp_cc$compsum)
    mus <- ddply(comp_cc,~mukey,summarise,propave=sum(propcwtd),compsum=mean(compsum))
    mus_min <- ddply(comp_cc,~mukey,summarise,propmin=min(propc))
    mus_max <- ddply(comp_cc,~mukey,summarise,propmax=max(propc))
    mus <- left_join(mus,mus_min,by="mukey")
    mus <- left_join(mus,mus_max,by="mukey")
    return(mus)
  }
  ## Run tabular steps in parallel list apply
  cl <- makeCluster(cpus, type="FORK")
  registerDoParallel(cl)
  mus_lst <- parLapply(cl,grps,try(propfn))
  stopCluster(cl)
  mus <- plyr::rbind.fill(mus_lst)

  ### Reclassify MU raster with new property values
  mu_tab <- mus[,c("mukey","propave","propmin","propmax")] # table to reclass by mukey
  mu_tab$propave <- mu_tab$propave * datastretch
  mu_tab$propmin <- mu_tab$propmin * datastretch
  mu_tab$propmax <- mu_tab$propmax * datastretch
  mukeys_na <- mukeys_tif[!mukeys_tif %in% mu_tab$mukey]
  mukeys_nadf <- data.frame(mukey=mukeys_na,propave=rep(NA,length(mukeys_na)),propmin=rep(NA,length(mukeys_na)),propmax=rep(NA,length(mukeys_na)))
  mu_tab <- rbind(mu_tab,mukeys_nadf)
  muave_tab <- mu_tab[,c("mukey","propave")]
  # mumin_tab <- mu_tab[,c("mukey","propmin")]
  # mumax_tab <- mu_tab[,c("mukey","propmax")]
  beginCluster(cpus,type='SOCK')
  proprast <- clusterR(sstif, reclassify, args=list(rcl=muave_tab))
  # minrast <- clusterR(sstif, reclassify, args=list(rcl=mumin_tab))
  # maxrast <- clusterR(sstif, reclassify, args=list(rcl=mumax_tab))
  endCluster()
  ## Save raster and reset environment for next depth
  writeRaster(proprast, overwrite=TRUE,filename=paste0(propdir,"/",prop,"_mean_",datastretchlab,"_",d,"_cm_gNATSGO.tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
  # writeRaster(minrast, overwrite=TRUE,filename=paste0(propdir,"/",prop,"_min_",datastretchlab,"_",d,"_cm_gNATSGO.tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
  # writeRaster(maxrast, overwrite=TRUE,filename=paste0(propdir,"/",prop,"_max_",datastretchlab,"_",d,"_cm_gNATSGO.tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
  rm(mus_lst,proprast,mus,mu_tab,mukeys_na,mukeys_nadf,d,minrast,maxrast,muave_tab,mumin_tab,mumax_tab) #horiz_c,horiz_c_depthwt,comp_c,comp_cc,comp_c_cmpct,
  gc()
}



################################ Create min, max and mean of map unit weighted mean property values across all depths for final covariates
#
# ### Bring in preped gNATSGO layers to summarize
# ss1 <- raster(paste0(propdir,"/",prop,"_mean_",datastretchlab,"_",depths[1],"_cm_gNATSGO.tif"))
# ss2 <- raster(paste0(propdir,"/",prop,"_mean_",datastretchlab,"_",depths[2],"_cm_gNATSGO.tif"))
# ss3 <- raster(paste0(propdir,"/",prop,"_mean_",datastretchlab,"_",depths[3],"_cm_gNATSGO.tif"))
# ss4 <- raster(paste0(propdir,"/",prop,"_mean_",datastretchlab,"_",depths[4],"_cm_gNATSGO.tif"))
# ss5 <- raster(paste0(propdir,"/",prop,"_mean_",datastretchlab,"_",depths[5],"_cm_gNATSGO.tif"))
# ss6 <- raster(paste0(propdir,"/",prop,"_mean_",datastretchlab,"_",depths[6],"_cm_gNATSGO.tif"))
# ss7 <- raster(paste0(propdir,"/",prop,"_mean_",datastretchlab,"_",depths[7],"_cm_gNATSGO.tif"))
#
# depthstk <- stack(ss1,ss2,ss3,ss4,ss5,ss6,ss7)
#
# ## Functions for raster summaries
# minfn <- function(x) min(x, na.rm =T)
# meanfn <- function(x) mean(x, na.rm =T)
# maxfn <- function(x) max(x, na.rm =T)
# NA_removefn <- function(x,y) {
#   ind <- ifelse(is.na(x)&!is.na(y), meanval,x)
#   return(ind)
# }
# NA_maskfn <- function(x,y){
#   ind <- ifelse(is.na(x), 0,1)
#   ind2 <- ifelse(is.na(y), NA,ind)
#   return(ind2)
# }
#
# ## Set up cluster
# beginCluster(cpus,type='SOCK')
#
# ## Calculate minimum value over all depths: fill in NAs with the mean value for use as SOLUS covarates
# minrast <- clusterR(depthstk,overlay, args=list(fun=minfn))
# nastack <- stack(minrast,sgtif)
# narast <- clusterR(nastack,overlay, args=list(fun=NA_maskfn),filename=paste0(propdir,"/",prop,"_NAs_gNATSGO.tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type)
# meanval <- cellStats(minrast, "mean")
# minstk <- stack(minrast,sgtif)
# minrast_na <- clusterR(minstk,overlay, args=list(fun=NA_removefn), export="meanval",filename=paste0(projfldr,"/",prop,"_min_",datastretchlab,"_gNATSGO.tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type)
# ## Calculate mean value over all depths
# meanrast <- clusterR(depthstk,overlay, args=list(fun=meanfn))
# meanval <- cellStats(meanrast, "mean")
# meanstk <- stack(meanrast,sgtif)
# meanrast_na <- clusterR(meanstk,overlay, args=list(fun=NA_removefn), export="meanval",filename=paste0(projfldr,"/",prop,"_mean_",datastretchlab,"_gNATSGO.tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type)
# ## Calculate max value over all depths
# maxrast <- clusterR(depthstk,overlay, args=list(fun=maxfn))
# meanval <- cellStats(maxrast, "mean")
# maxstk <- stack(maxrast,sgtif)
# maxrast_na <- clusterR(maxstk,overlay, args=list(fun=NA_removefn), export="meanval",filename=paste0(projfldr,"/",prop,"_max_",datastretchlab,"_gNATSGO.tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type)
#
#
# endCluster()
#



