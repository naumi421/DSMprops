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
predfolder <- "/push/HYBconus100m/Clay"
repofolder <- "/push/repos/DSMprops/exec/conus100m/Clay"
covfolder <- "/push/SG100_covars"
ptsfolder <- "/push/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final"

## Points
#demo(meuse, echo=FALSE)
#data(BostonHousing, package = "mlbench")
pts <- readRDS("/push/HYBconus100m/Clay/TrainPTS_claytotal_r_0_cm.rds")

## Make list of grids
cov.grids <- list.files(path = covfolder,pattern=".tif$",full.names = T, recursive = F)
cov.grids.names <- basename(cov.grids)

## Trim training df
ptspred.list <- gsub(".tif","", cov.grids.names)# Take .tif off of the grid list to just get the variable names
ptspred.list <- c(ptspred.list,"prop")
pts_trn <- pts[,c(ptspred.list)]

## Build MLR ensemble model
parallelMap::parallelStart(mode="multicore",cpus=55)
tsk = makeRegrTask(data = pts_trn@data, target = "prop")
base = c("regr.ranger", "regr.xgboost", "regr.nnet", "regr.ksvm", "regr.cvglmnet")
lrns = lapply(base, makeLearner)
m <- makeStackedLearner(base.learners = lrns, super.learner = "regr.lm", method = "stack.cv")
tmp = train(m, tsk) # Rsq = 0.7266, took about 11 hrs with n = 99k
res = predict(tmp, tsk)
parallelStop()
## Save files
saveRDS(tmp, paste(predfolder,"/mlr_stk_learner_tst_clay0cm.rds",sep=""))


## Prep grids for landmap
rasters <- stack(cov.grids)
sys.time()
r_spdf <-  as(rasters, "SpatialPixelsDataFrame")
sys.time()

