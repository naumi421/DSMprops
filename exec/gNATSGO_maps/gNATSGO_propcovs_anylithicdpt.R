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
rasterOptions(maxmemory = 6.5e+09, chunksize = 6e+08) # reset raster params
sstif <- raster(paste0(gnatsgofldr,"/mukeys_sg100grid_gNATSGO.tif"))

#### Get unique mukeys from map
mukeys_tif <- raster::unique(sstif)

## Legend table
legend.df <- read.csv(paste0(gnatsgotabfldr,"/legend.csv"))
## Mapunit table
mapunit.df <- read.csv(paste0(gnatsgotabfldr,"/mapunit.csv"))
mukeys <- mukeys_tif[mukeys_tif %in% mapunit.df$mukey] # Eliminates mukeys that have no map unit table entries and thus no component data when adding zero values later

## Check extent in mapview
# mapview(sstif)

##### Also load necessary data tables and subset to study area
## components
comp.df <- read.csv(paste0(gnatsgotabfldr,"/component.csv"))
comp.df$compname <- toupper(comp.df$compname)
## subset tables down to area of interest
comps <- comp.df[comp.df$mukey %in% mukeys,]
cokeys <- comps$cokey

## Attribute common Misc area components
## Pick out rock outcrop components (from Colby Brungard and Shawn Salley)
# subset out bedrock type components for attribution
CMU.5 <- comps[grepl(comps$compname, pattern = "rock outcrop", ignore.case = TRUE),]
CMU.7 <- comps[grepl(comps$compname, pattern = "outcrop", ignore.case = TRUE),]
CMU.9 <- comps[grepl(comps$compname, pattern = "bedrock", ignore.case = TRUE),]
CMU.10 <- do.call("rbind", list(CMU.5, CMU.7, CMU.9))
CMU.11 <- unique(CMU.10)
CMU.12 <- CMU.11[!(CMU.11$compkind %in% c("Series", "Taxon above family", "Family", "Taxadjunct", "Variant")), ] 

## Restriction table
restr.df <- read.csv(paste0(gnatsgotabfldr,"/corestrictions.csv"))
restrs <- restr.df[restr.df$cokey %in% cokeys,]
## Order restrs to get first restriction for each category
restrs.all <- restrs[with(restrs,order(cokey,resdept_r)),]
## Depth paralithic or bedrock lithics contact
restrs.anyLithic <- subset(restrs.all, restrs.all$reskind == "Paralithic bedrock"|restrs.all$reskind=="Lithic bedrock")
comp_depth <- left_join(comps, restrs.anyLithic, by = "cokey")
comp_depth <- comp_depth[with(comp_depth,order(cokey,resdept_r)),] # Weed out some duplicates
comp_depth <- comp_depth[!duplicated(comp_depth$cokey),] # Weed out some duplicates
## Adds max obs depth values of 201 where no restriction exists (as per sdvattribute nullratingreplacementvalue)
comp_depth$resdept_r <- ifelse(comp_depth$cokey %in% CMU.12$cokey, 0, comp_depth$resdept_r)
comp_depth$resdept_r <- ifelse((is.na(comp_depth$resdept_r)) & (comp_depth$compkind != "Miscellaneous area"), 201, comp_depth$resdept_r)
comp_depth$resdept_r <- ifelse(comp_depth$resdept_r > 201, 201, comp_depth$resdept_r) # REclass a couple >201 to 201 for consistent right censorship
length(unique(comp_depth$cokey)) ## Check for duplicate cokeys, should be none (same number as the obs in df)


## Set gNATSGO property
prop <- 'anylithicdpt'
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
comp_depth <- left_join(comp_depth, mukeys.df, by = "mukey")
grps <- unique(comp_depth$grp)
rasterOptions(maxmemory = 6e+09, chunksize = 8e+08)

d <- "all"

## Set up function for gNATSGO MU weighted average calcs
propfn <- function(g){
  ## Create a cleaned up df
  comp_depthc <- comp_depth[comp_depth$grp == g,]
  comp_c <- ddply(comp_depthc,~cokey,summarise,propc=min(resdept_r),mukey=mean(as.numeric(mukey)),comppct_r=mean(comppct_r)) 

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





