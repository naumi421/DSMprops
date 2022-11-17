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
projfldr <- "/mnt/solus100/Predictionsv2/Predictionsv2tst/AWC_gRPI/compare"
gnatsgofldr <- "/mnt/covs/gNATSGO_CONUS_tif"
predfldr <- "/mnt/solus100/Predictionsv2/Predictionsv2tst/AWC_gRPI"
gnatsgotabfldr <- "/mnt/covs/gNATSGO_CONUS_tif/gNATSGO_Tabular_CSV"


### gNATSGO tif
sstif <- raster(paste0(gnatsgofldr,"/gNATSGO-mukey.tif"))
ss.proj <- projection(sstif)
ss.wkt <- sf::st_crs(ss.proj)
ss.crs <- raster::crs(sstif) ## Maintains both proj4 and proj6 properties


### Area of Interest
areaname <- "chico"
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

## Set gNATSGO property
prop <- 'awc_r'
## Vector of depths to render
depths <- c(0,5,15,30,60,100,150)
## Scaling and datatype
data_type <- "INT1U"
datastretch <- 100
datastretchlab <- paste(datastretch,"x",sep="")

## Set up depth loop for gNATSGO rasters
# d <- 60
for(d in depths){
  ## Prep property map to a depth interval
  ud <- d - 0.5 # upper depth
  ld <- d + 0.5 # lower depth
  horiz_c <- subset(horiz_comps, (as.numeric(horiz_comps$hzdept_r) <= ud & as.numeric(horiz_comps$hzdepb_r) > ud) |
                      (as.numeric(horiz_comps$hzdepb_r) >= ld & as.numeric(horiz_comps$hzdept_r) < ld) ) # subset to horizon overlapping range
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

  ### Reclassify MU raster with new property values
  mu_tab <- mus[,1:2] # table to reclass by mukey
  mu_tab$propave <- mu_tab$propave * datastretch
  mukeys_na <- mukeys[!mukeys %in% mu_tab$mukey]
  mukeys_nadf <- data.frame(mukey=mukeys_na,propave=rep(NA,length(mukeys_na)))
  mu_tab <- rbind(mu_tab,mukeys_nadf)
  proprast <- reclassify(sstif_clp, mu_tab)
  ## Save raster and reset environment for next depth
  writeRaster(proprast, overwrite=TRUE,filename=paste0(projfldr,"/",prop,"_",datastretchlab,"_",d,"_cm_gNATSGO_",areaname,".tif"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype=data_type, progress="text")
  rm(proprast,horiz_c,horiz_c_depthwt,comp_c,comp_cc,comp_c_cmpct,mus,mu_tab,mukeys_na,mukeys_nadf)
  gc()
}

### Bring in other raster predictions for comparison
lithic <- crop(raster("/mnt/solus100/Predictions/RestrDepth_gRPI/anylithicdpt_1x_all_cm_2D_QRF.tif"),aoi_ext) ## Depth to bedrock raster to clip
rast.files <- list.files(path = predfldr,pattern="QRFadj.tif$",full.names = T, recursive = F)
pred1 <- crop(raster(rast.files[grep(paste0('_',depths[1],'_cm'),rast.files)]),aoi_ext)
pred2 <- crop(raster(rast.files[grep(paste0('_',depths[2],'_cm'),rast.files)]),aoi_ext)
pred3 <- crop(raster(rast.files[grep(paste0('_',depths[3],'_cm'),rast.files)]),aoi_ext)
pred4 <- crop(raster(rast.files[grep(paste0('_',depths[4],'_cm'),rast.files)]),aoi_ext)
pred5 <- crop(raster(rast.files[grep(paste0('_',depths[5],'_cm'),rast.files)]),aoi_ext)
pred6 <- crop(raster(rast.files[grep(paste0('_',depths[6],'_cm'),rast.files)]),aoi_ext)
pred7 <- crop(raster(rast.files[grep(paste0('_',depths[7],'_cm'),rast.files)]),aoi_ext)

### Use mask function to assign pixels below bedrock depth to NA
pred1d <- overlay(lithic,pred1, fun=function(l,r) {return(ifelse(l>depths[1],r,NA))},progress='text',filename = paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[1],"_cm_SOLUS_",areaname,".tif"),datatype=data_type, overwrite=TRUE)
pred2d <- overlay(lithic,pred2, fun=function(l,r) {return(ifelse(l>depths[2],r,NA))},progress='text',filename = paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[2],"_cm_SOLUS_",areaname,".tif"),datatype=data_type, overwrite=TRUE)
pred3d <- overlay(lithic,pred3, fun=function(l,r) {return(ifelse(l>depths[3],r,NA))},progress='text',filename = paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[3],"_cm_SOLUS_",areaname,".tif"),datatype=data_type, overwrite=TRUE)
pred4d <- overlay(lithic,pred4, fun=function(l,r) {return(ifelse(l>depths[4],r,NA))},progress='text',filename = paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[4],"_cm_SOLUS_",areaname,".tif"),datatype=data_type, overwrite=TRUE)
pred5d <- overlay(lithic,pred5, fun=function(l,r) {return(ifelse(l>depths[5],r,NA))},progress='text',filename = paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[5],"_cm_SOLUS_",areaname,".tif"),datatype=data_type, overwrite=TRUE)
pred6d <- overlay(lithic,pred6, fun=function(l,r) {return(ifelse(l>depths[6],r,NA))},progress='text',filename = paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[6],"_cm_SOLUS_",areaname,".tif"),datatype=data_type, overwrite=TRUE)
pred7d <- overlay(lithic,pred7, fun=function(l,r) {return(ifelse(l>depths[7],r,NA))},progress='text',filename = paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[7],"_cm_SOLUS_",areaname,".tif"),datatype=data_type, overwrite=TRUE)

### Bring in preped gNATSGO layers to plot
ss1 <- raster(paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[1],"_cm_gNATSGO_",areaname,".tif"))
ss2 <- raster(paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[2],"_cm_gNATSGO_",areaname,".tif"))
ss3 <- raster(paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[3],"_cm_gNATSGO_",areaname,".tif"))
ss4 <- raster(paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[4],"_cm_gNATSGO_",areaname,".tif"))
ss5 <- raster(paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[5],"_cm_gNATSGO_",areaname,".tif"))
ss6 <- raster(paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[6],"_cm_gNATSGO_",areaname,".tif"))
ss7 <- raster(paste0(projfldr,"/",prop,"_",datastretchlab,"_",depths[7],"_cm_gNATSGO_",areaname,".tif"))

### Now create trellis graph to compare the datasets
model_maps<-stack(pred1d,pred2d,pred3d,pred4d,pred5d,pred6d,pred7d)
names(model_maps)<-list( paste0('Modeled ',prop,' ',depths[1],' cm'), paste0('Modeled ',prop,' ',depths[2],' cm'), paste0('Modeled ',prop,' ',depths[3],' cm'), paste0('Modeled ',prop,' ',depths[4],' cm'),
                   paste0('Modeled ',prop,' ',depths[5],' cm'),  paste0('Modeled ',prop,' ',depths[6],' cm'), paste0('Modeled ',prop,' ',depths[7],' cm'))

ss_maps<-stack(ss1,ss2,ss3,ss4,ss5,ss6,ss7)
names(ss_maps)<-list( paste0('gNATSGO ',prop,' ',depths[1],' cm'), paste0('gNATSGO ',prop,' ',depths[2],' cm'),paste0('gNATSGO ',prop,' ',depths[3],' cm'), paste0('gNATSGO ',prop,' ',depths[4],' cm'),
                         paste0('gNATSGO ',prop,' ',depths[5],' cm'), paste0('gNATSGO ',prop,' ',depths[6],' cm'), paste0('gNATSGO ',prop,' ',depths[7],' cm'))


## Ramp and Label parameters
my.at<-c(0,5,10,15,20,25,30,35) # color ramp breaks used for labels (will be excluded if colorkey truncates)
my.atc <- c(0,5,10,15,20,25,30,35) # color ramp break label positions
my.colorkey<-list(at=my.atc,
                  labels=list(
                    at=my.at
                  ))

### Actual graphing block
## Modeled plot
tiff(paste0(projfldr,"/",prop,"_SOLUS_",areaname,".tif"),width = 3.5, height = 14, units = 'in', res = 600 ) # begins rendering high res figure
p <- levelplot(model_maps, par.settings=viridisTheme,  scales=list(draw=F),colorkey=my.colorkey,at=my.at, main="",layout=c(1,7) ) #margin=list(FUN='median'),
p
dev.off()

## gNATSGO plot
tiff(paste0(projfldr,"/",prop,"_gNATSGO_",areaname,".tif"),width = 3.5, height = 14, units = 'in', res = 600 ) # begins rendering high res figure
p <- levelplot(ss_maps, par.settings=viridisTheme,  scales=list(draw=F),colorkey=my.colorkey,at=my.at, main="",layout=c(1,7) ) #margin=list(FUN='median'),
p
dev.off()

