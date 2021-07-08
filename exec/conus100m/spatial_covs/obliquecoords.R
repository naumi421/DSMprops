##############################
## Creation of oblique coordinate
## covariate rasters for
## the 100m CONUS dataset
## T. Nauman, 7/7/2021

# Install packages if not already installed
required.packages <- c("raster", "sp", "rgdal", "OGC", "doParallel")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)

## Key folders
covfolder <- "/push/SG100_covars"
ogcfolder <- "/push/ogc_covs"

## Base rater
projgrid <- raster(paste(covfolder,"/T07PRI5.tif",sep="")) # This raster is fully clipped

## Make oblique coordinate rasters: going with 30 based on results from Moller et al 2020,. https://doi.org/10.5194/soil-6-269-2020
n <- 5
ogcs <- makeOGC(projgrid, n) # Creates hard break in middle of USA???

## Function to split up raster brick and compress each OGC raster and save in covariate folder
ogc_splt_fn <- function(x){
  r <- ogcs[[x]]
  filenm <- names(r)
  writeRaster(r, overwrite=TRUE,filename=paste(ogcfolder,"/",filenm,".tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype="INT4S")
}

## Run parallel list apply to save rasters
xlst <- 1:n
cpus <- n
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
pretime <- Sys.time()
parLapply(cl,xlst,ogc_splt_fn)
posttime <- Sys.time()
runtime <- posttime - pretime

