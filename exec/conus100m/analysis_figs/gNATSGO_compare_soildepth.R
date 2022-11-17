###############################################
## Script to compare SOLUS property predictions
## to current gNATSGO weighted average mapunit
## property estimates
## T Nauman, 11/16/2022

### Packages
required.packages <- c("raster", "sp", "rgdal","doParallel", "plyr", "dplyr","sf","cluster","RSQLite","lubridate","rasterVis","maptools","classInt")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
##Raster settings
rasterOptions(maxmemory = 2e+09, chunksize = 2e+08)

### Key folders
projfldr <- "/mnt/disks/sped/solus100_analysis/comparisons"
gnatsgofldr <- "/mnt/covs/gNATSGO_CONUS_tif"
predfldr <- "/mnt/solus100/Predictionsv2/Predictionsv2tst/AWC_gRPI"
gnatsgotabfldr <- "/mnt/covs/gNATSGO_CONUS_tif/gNATSGO_Tabular_CSV"

### gNATSGO tif
sstif <- raster(paste0(gnatsgofldr,"/gNATSGO-mukey.tif"))
ss.proj <- projection(sstif)
ss.wkt <- sf::st_crs(ss.proj)
ss.crs <- raster::crs(sstif) ## Maintains both proj4 and proj6 properties
ss_rat <- #load dbf


### Area of Interest
## Using bbox
# aoi_ext <- extent(xmin,xmax,ymin,ymax) ## enter xmin,xmax,ymin,ymax in coordinate system
# proj4string(aoi_ext) <- sp::CRS("EPSG:4326") # PROJ6 compatible
## Using sp object
layernm <- "chico_solus_poly"
polybound <- readOGR(projfldr, layernm)
polybound_p <- spTransform(polybound, ss.crs)
aoi_ext <- extent(polybound_p)

#### Clip raster down and get unique mukeys
sstif_clp <- crop(sstif, aoi_ext)
mukeys <- raster::unique(sstif_clp)

##### Also load necessary data tables and subset to study area
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


#### Component level properties: e.g, Depth, frag cover
## Housekeeping for depth restrictions
restrs.all <- restrs[with(restrs,order(cokey,resdept_r)),] ## Need to order to pick first
restrs.anyLithic <- subset(restrs.all, restrs.all$reskind == "Paralithic bedrock"|restrs.all$reskind=="Lithic bedrock")
comp_restrdepth <- left_join(comps,restrs.all,by="cokey") # First restrictive layer of anykind
comp_depth <- left_join(comps,restrs.anyLithic,by="cokey") # First bedrock related layer
## Adds max obs depth values of 201 where no restriction exists (as per sdvattribute nullratingreplacementvalue)
comp_restrdepth$resdept_r <- ifelse(is.na(comp_restrdepth$resdept_r), 201, comp_restrdepth$resdept_r)
comp_restrdepth$resdept_r <- ifelse(comp_restrdepth$resdept_r > 201, 201, comp_restrdepth$resdept_r) # REclass a couple >201 to 201 for consistent right censorship
comp_depth$resdept_r <- ifelse(is.na(comp_depth$resdept_r), 201, comp_depth$resdept_r)
comp_depth$resdept_r <- ifelse(comp_depth$resdept_r > 201, 201, comp_depth$resdept_r) # REclass a couple >201 to 201 for consistent right censorship
comp_depthc <- ddply(comp_depth,~cokey,summarise,propc=sum(resdept_r),mukey=mean(as.numeric(mukey)),comppct_r=mean(comppct_r))
comp_restrdepthc <- ddply(comp_restrdepth,~cokey,summarise,propc=sum(resdept_r),mukey=mean(as.numeric(mukey)),comppct_r=mean(comppct_r))

##### Reassign property of interest to compc object for map unit summarization: Can only do one at a time thus far
comp_c <- comp_depthc # bedrock depth
comp_c2 <- comp_restrdepthc # restriction depth

### Now summarize comps at mu level: bedrock depth
comp_c_cmpct <- comp_c[,c("mukey","comppct_r")]
comp_c_cmpct <- ddply(comp_c_cmpct,~mukey,summarise,compsum=sum(comppct_r)) # determine total % named components in MU
comp_c <- merge(comp_c, comp_c_cmpct,by="mukey")
comp_cc <- subset(comp_c, comp_c$compsum > 10)
comp_cc$propcwtd <- comp_cc$propc * (comp_cc$comppct_r / comp_cc$compsum)
mus <- ddply(comp_cc,~mukey,summarise,propave=sum(propcwtd),compsum=mean(compsum))

### Now summarize comps at mu level: restriction depth
comp_c2_cmpct <- comp_c2[,c("mukey","comppct_r")]
comp_c2_cmpct <- ddply(comp_c2_cmpct,~mukey,summarise,compsum=sum(comppct_r)) # determine total % named components in MU
comp_c2 <- merge(comp_c2, comp_c2_cmpct,by="mukey")
comp_c2c <- subset(comp_c2, comp_c2$compsum > 10)
comp_c2c$propcwtd <- comp_c2c$propc * (comp_c2c$comppct_r / comp_c2c$compsum)
mus_r <- ddply(comp_c2c,~mukey,summarise,propave=sum(propcwtd),compsum=mean(compsum))

### Reclassify MU raster with new property values: Bedrock depth
mu_tab <- mus[,1:2] # table to reclass by mukey
mukeys_na <- mukeys[!mukeys %in% mu_tab$mukey]
mukeys_nadf <- data.frame(mukey=mukeys_na,propave=rep(NA,length(mukeys_na)))
mu_tab <- rbind(mu_tab,mukeys_nadf)
proprast <- reclassify(sstif_clp, mu_tab)

### Reclassify MU raster with new property values: restriction depth
mu_r_tab <- mus_r[,1:2] # table to reclass by mukey
mukeys_r_na <- mukeys[!mukeys %in% mu_r_tab$mukey]
mukeys_r_nadf <- data.frame(mukey=mukeys_r_na,propave=rep(NA,length(mukeys_r_na)))
mu_r_tab <- rbind(mu_r_tab,mukeys_r_nadf)
proprast_r <- reclassify(sstif_clp, mu_r_tab)





