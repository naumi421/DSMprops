## Test out Landmap and MLR




# Workspace setup
# Install packages if not already installed
required.packages <- c("rgdal", "raster", "plotKML", "geoR", "ranger", "mlr", "forestError",
                       "xgboost", "glmnet", "matrixStats", "kernlab", "deepnet","landmap","mlbench","parallelMap")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package: Windows only
#memory.limit(500000)
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08)

## Key Folder Locations
predfolder <- "/push/Hyb100m_gdrv/Clay"
repofolder <- "/push/repos/DSMprops/exec/conus100m/Clay"
covfolder <- "/push/SG100_covars"
ptsfolder <- "/push/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final"

## Points
#demo(meuse, echo=FALSE)
#data(BostonHousing, package = "mlbench")
pts <- readRDS(paste0(predfolder,"/TrainPTS_claytotal_r_0_cm.rds"))

## Make list of grids
cov.grids <- list.files(path = covfolder,pattern=".tif$",full.names = T, recursive = F)
cov.grids.names <- basename(cov.grids)

## Trim training df
ptspred.list <- gsub(".tif","", cov.grids.names)# Take .tif off of the grid list to just get the variable names
ptspred.list <- c(ptspred.list,"prop")
pts_trn <- pts[,c(ptspred.list)]

## Create spatial blocking grid
img10kf <- readRDS(paste(predfolder,"/img10kf.rds",sep=""))
crs.rast <- CRS(projection(img10kf))
img10kf <- projectRaster(img10kf,res=50000,method='ngb',crs=crs.rast) # resample to 50km
## Create another grid with unique values for each pixel excepting nodata areas
## This is for use in the point density mapping
img10kfid <- img10kf
values(img10kfid) <- 1:ncell(img10kfid)
values(img10kfid) <- sample.int(10, size = ncell(img10kfid), replace=T) # replicates what I've done.
img10kfid <- overlay(stack(img10kfid,img10kf),fun=function(a,b){a*b})
pts_trn_blk <- raster::extract(img10kfid, pts_trn, df=T, sp=T)
pts_trn_blk <- pts_trn_blk@data
pts_trn_blk <- pts_trn_blk[!is.na(pts_trn_blk$layer),]

### Missing peice with Hengl: the grid id from grid topology: See splearner function: start at line 87


## Build MLR model an CV with blocking on grid ids
parallelMap::parallelStart(mode="multicore",cpus=10)
pts_trn_blk$layer <- as.factor(pts_trn_blk$layer)
tsk <- makeRegrTask(data = pts_trn_blk, target = "prop", blocking = pts_trn_blk$layer)
rdesc = makeResampleDesc("CV", iters = 10, blocking.cv = TRUE)
p = resample("regr.ranger", tsk, rdesc, par.vals = list(num.trees = 100, min.node.size = 1,num.threads = 5))
rsq <- measureRSQ(p$pred$data$truth,p$pred$data$response)
plot(p$pred$data$truth~p$pred$data$response)
parallelStop()

#base = c("regr.ranger", "regr.xgboost", "regr.nnet", "regr.ksvm", "regr.cvglmnet")
# lrns = lapply(base, makeLearner)
# m <- makeStackedLearner(base.learners = lrns, super.learner = "regr.lm", method = "stack.cv")
# tmp = train(m, tsk) # Rsq = 0.7266, took about 11 hrs with n = 99k
# res = predict(tmp, tsk)
#
# ## Save files
# saveRDS(tmp, paste(predfolder,"/mlr_stk_learner_tst_clay0cm.rds",sep=""))
#
#
# ## Prep grids for landmap
# rasters <- stack(cov.grids)
# sys.time()
# r_spdf <-  as(rasters, "SpatialPixelsDataFrame")
# sys.time()

