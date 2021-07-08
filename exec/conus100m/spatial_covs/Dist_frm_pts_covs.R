##############################
## Creation of point distance
## covariate rasters for
## the 100m CONUS dataset
## Uses n random points to create n distance
## rasters to use as spatial covariates in model
## T. Nauman, 7/7/2021

# Install packages if not already installed
required.packages <- c("raster", "sp", "rgdal",  "doParallel")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)

## Key folders
covfolder <- "/push/SG100_covars"
ptdistfolder <- "/push/pt_dist_covs"

## Base rater
projgrid <- raster(paste(covfolder,"/T07PRI5.tif",sep="")) # This raster is fully clipped
cov.proj <- projection(projgrid)

## AOI
bnd <- readOGR("/ped/GIS_Archive/US_boundaries/states_21basic", "conus_bound")
bnd <- spTransform(bnd, cov.proj)

## Create random points within AOI
n <- 30
pts <- spsample(bnd[1,], n, type = 'nonaligned')
pts <- SpatialPointsDataFrame(pts, data.frame(row.names=row.names(pts), ID=1:length(pts)))
saveRDS(pts, paste(ptdistfolder,"/pts_fc.rds",sep=""))

## Function to split up points and create distance rasters for each point
pt_dist_fn <- function(x){
  pt <- pts[x,]
  r <- distanceFromPoints(projgrid, pt)
  r <- r / 100 # Scale for compression
  filenm <- paste("pt",x,"dist.tif",sep="")
  writeRaster(r, overwrite=TRUE,filename=paste(ptdistfolder,"/",filenm,sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype="INT2U")
}

## Run parallel list apply to save rasters
xlst <- 1:nrow(pts)
cpus <- 32
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
pretime <- Sys.time()
parLapply(cl,xlst,pt_dist_fn)
posttime <- Sys.time()
runtime <- posttime - pretime

