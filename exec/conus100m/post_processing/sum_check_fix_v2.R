## Function to normalize sand silt and clay summed predictions to 100%


# Install packages if not already installed
required.packages <- c("raster","sp", "doParallel")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)

## Folders
oldfldr <- "/mnt/disks/sped/ssc_sum/orig_preds"
newfldr <- "/mnt/disks/sped/ssc_sum"

# Clay
cl0   <- raster("/mnt/disks/sped/ssc_sum/orig_preds/claytotal_r_1x_0_cm_2D_QRFadj.tif")
cl5   <- raster("/mnt/disks/sped/ssc_sum/orig_preds/claytotal_r_1x_5_cm_2D_QRFadj.tif")
cl15  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/claytotal_r_1x_15_cm_2D_QRFadj.tif")
cl30  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/claytotal_r_1x_30_cm_2D_QRFadj.tif")
cl60  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/claytotal_r_1x_60_cm_2D_QRFadj.tif")
cl100 <- raster("/mnt/disks/sped/ssc_sum/orig_preds/claytotal_r_1x_100_cm_2D_QRFadj.tif")
cl150 <- raster("/mnt/disks/sped/ssc_sum/orig_preds/claytotal_r_1x_150_cm_2D_QRFadj.tif")

# Sand
s0   <- raster("/mnt/disks/sped/ssc_sum/orig_preds/sandtotal_r_1x_0_cm_2D_QRFadj.tif")
s5   <- raster("/mnt/disks/sped/ssc_sum/orig_preds/sandtotal_r_1x_5_cm_2D_QRFadj.tif")
s15  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/sandtotal_r_1x_15_cm_2D_QRFadj.tif")
s30  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/sandtotal_r_1x_30_cm_2D_QRFadj.tif")
s60  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/sandtotal_r_1x_60_cm_2D_QRFadj.tif")
s100 <- raster("/mnt/disks/sped/ssc_sum/orig_preds/sandtotal_r_1x_100_cm_2D_QRFadj.tif")
s150 <- raster("/mnt/disks/sped/ssc_sum/orig_preds/sandtotal_r_1x_150_cm_2D_QRFadj.tif")

# Silt
si0   <- raster("/mnt/disks/sped/ssc_sum/orig_preds/silttotal_r_1x_0_cm_2D_QRFadj.tif")
si5   <- raster("/mnt/disks/sped/ssc_sum/orig_preds/silttotal_r_1x_5_cm_2D_QRFadj.tif")
si15  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/silttotal_r_1x_15_cm_2D_QRFadj.tif")
si30  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/silttotal_r_1x_30_cm_2D_QRFadj.tif")
si60  <- raster("/mnt/disks/sped/ssc_sum/orig_preds/silttotal_r_1x_60_cm_2D_QRFadj.tif")
si100 <- raster("/mnt/disks/sped/ssc_sum/orig_preds/silttotal_r_1x_100_cm_2D_QRFadj.tif")
si150 <- raster("/mnt/disks/sped/ssc_sum/orig_preds/silttotal_r_1x_150_cm_2D_QRFadj.tif")



ssc0 <- stack(s0, si0, cl0)
names(ssc0) <- c('sand', 'silt', 'clay')
sum0 <- sum(ssc0)
writeRaster(sum0, "/mnt/disks/sped/ssc_sum/sum0.tif")

ssc5 <- stack(s5, si5, cl5)
names(ssc5) <- c('sand', 'silt', 'clay')
sum5 <- sum(ssc5)
writeRaster(sum5, "/mnt/disks/sped/ssc_sum/sum5.tif")

ssc15 <- stack(s15, si15, cl15)
names(ssc15) <- c('sand', 'silt', 'clay')
sum15 <- sum(ssc15)
writeRaster(sum15, "/mnt/disks/sped/ssc_sum/sum15.tif")

ssc30 <- stack(s30, si30, cl30)
names(ssc30) <- c('sand', 'silt', 'clay')
sum30 <- sum(ssc30)
writeRaster(sum30, "/mnt/disks/sped/ssc_sum/sum30.tif")

ssc60 <- stack(s60, si60, cl60)
names(ssc60) <- c('sand', 'silt', 'clay')
sum60 <- sum(ssc60)
writeRaster(sum60, "/mnt/disks/sped/ssc_sum/sum60.tif")

ssc100 <- stack(s100, si100, cl100)
names(ssc100) <- c('sand', 'silt', 'clay')
sum100 <- sum(ssc100)
writeRaster(sum100, "/mnt/disks/sped/ssc_sum/sum100.tif")

ssc150 <- stack(s150, si150, cl150)
names(ssc150) <- c('sand', 'silt', 'clay')
sum150 <- sum(ssc150)
writeRaster(sum150, "/mnt/disks/sped/ssc_sum/sum150.tif")

## Files to include in function
oldfiles <- list.files(path = oldfldr,pattern=".tif$",full.names = T, recursive = T)

## Depths for list apply
depths <- c("_0_cm","_5_cm","_15_cm","_30_cm","_60_cm","_100_cm","_150_cm")

## Function for normalizing fractions
divide_by_sum <- function(d) {
  a <- raster(paste0(oldfldr,"/sandtotal_r_1x",d,"_2D_QRFadj.tif"))
  b <- raster(paste0(oldfldr,"/silttotal_r_1x",d,"_2D_QRFadj.tif"))
  c <- raster(paste0(oldfldr,"/claytotal_r_1x",d,"_2D_QRFadj.tif"))
  depthno <- as.numeric(unlist(strsplit(d, "_"))[2])
  # no values sum to 0 (and few that sum to exactly 100), but this is how a check could be incorporated
  #if (sum(vector) == 0) {
    #stop("Sum of the vector is zero. Division by zero is not allowed.")
  #}
  sand <- a / sum(a,b,c)*100
  silt <- b / sum(a,b,c)*100
  clay <- c / sum(a,b,c)*100
  #v <- vector / sum(vector)*100
  sand <- round(sand,0)
  silt <- round(silt,0)
  clay <- round(clay,0)
  #v2 <- round(v, 0)
  v2 <- brick(sand,silt,clay)
  ## rounding creates occasional sums greater than 100
  ## Errant values are either 99 or 101, so we are correcting by subtracting 1 from
  ## the texture class with the max value for 101 pixels, or adding 1 to the minimum
  ## texture fraction for sums of 99
  # Create rasters showing where errors are
  sumv2 <- sum(sand,silt,clay)
  sumv2_101 <- sumv2
  sumv2_101[sumv2_101 < 101] <- 0
  sumv2_101[sumv2_101 > 100] <- 1
  sumv2_99 <- sumv2
  sumv2_99[sumv2_99 < 99 | sumv2_99 > 99] <- 0
  sumv2_99[sumv2_99 == 99] <- 1
  ## Reference which fraction is max and min in stacks
  maxv2 <- which.max(v2)
  minv2 <- which.min(v2)
  # maxv2 <- terra::app(v2, which.max)
  # minv2 <- terra::app(v2, which.min)
  ## Create rasters for where clay values should be modified
  claymax <- maxv2
  claymax[claymax < 3] <- 0
  claymax[claymax == 3] <- 1
  claymax <- claymax * sumv2_101
  claymin <- minv2
  claymin[claymin < 3] <- 0
  claymin[claymin == 3] <- 1
  claymin <- claymin * sumv2_99
  ## Create rasters for where silt values should be modified
  siltmax <- maxv2
  siltmax[siltmax < 2 | siltmax > 2] <- 0
  siltmax[siltmax == 2] <- 1
  siltmax <- siltmax * sumv2_101
  siltmin <- minv2
  siltmin[siltmin < 2 | siltmin > 2] <- 0
  siltmin[siltmin == 2] <- 1
  siltmin <- siltmin * sumv2_99
  ## Create rasters for where sand values should be modified
  sandmax <- maxv2
  sandmax[sandmax > 1] <- 0
  sandmax <- sandmax * sumv2_101
  sandmin <- minv2
  sandmin[sandmin > 1] <- 0
  sandmin <- sandmin * sumv2_99
  gc()
  clay <- v2[["clay"]]
  clay <- clay + claymin
  clay <- clay - claymax
  silt <- v2[["silt"]]
  silt <- silt + siltmin
  silt <- silt - siltmax
  gc()
  sand <- v2[["sand"]]
  sand <- sand + sandmin
  sand <- sand - sandmax
  v3 <- stack(sand, silt, clay)
  names(v3) <- c('sand', 'silt', 'clay')
  writeRaster(v3, filename = paste0(newfldr,"/sum",depthno,"_fix.tif"), overwrite=TRUE,datatype = 'INT1U', options=c("COMPRESS=DEFLATE", "TFW=NO"))
  #return(v3)
}

## Linux parallel list apply
cpus <- min(length(depths),detectCores() - 2)
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
parLapply(cl,depths,try(divide_by_sum))
stopCluster(cl)

## Set up to run function in parallel
# rasterOptions(maxmemory = 6.5e+09,chunksize = 5e+08, memfrac = 0.9)
# beginCluster(30,type='SOCK')

# clusterR(ssc0, overlay,filename= "/mnt/disks/sped/ssc_sum/sum0_fix.tif",
#     overwrite=TRUE, args = list(fun = divide_by_sum),datatype = 'INT1U',progress="text")

# sum5_fix <- clusterR(ssc5, overlay, args = list(fun = divide_by_sum),progress="text")
#filename = "/mnt/disks/sped/ssc_sum/sum5_fix.tif", overwrite=TRUE,datatype = 'INT1U', options=c("COMPRESS=DEFLATE", "TFW=YES")

# app(ssc5, divide_by_sum, cores=120, filename= "/mnt/disks/sped/ssc_sum/sum5_fix.tif", overwrite=TRUE, wopt=list(datatype='INT1U', names = c('sand', 'silt', 'clay')))
# app(ssc15, divide_by_sum, cores=120, filename= "/mnt/disks/sped/ssc_sum/sum15_fix.tif", overwrite=TRUE, wopt=list(datatype='INT1U', names = c('sand', 'silt', 'clay')))
# app(ssc30, divide_by_sum, cores=120, filename= "/mnt/disks/sped/ssc_sum/sum30_fix.tif", overwrite=TRUE, wopt=list(datatype='INT1U', names = c('sand', 'silt', 'clay')))
# app(ssc60, divide_by_sum, cores=120, filename= "/mnt/disks/sped/ssc_sum/sum60_fix.tif", overwrite=TRUE, wopt=list(datatype='INT1U', names = c('sand', 'silt', 'clay')))
# app(ssc100, divide_by_sum, cores=120, filename= "/mnt/disks/sped/ssc_sum/sum100_fix.tif", overwrite=TRUE, wopt=list(datatype='INT1U', names = c('sand', 'silt', 'clay')))
# app(ssc150, divide_by_sum, cores=120, filename= "/mnt/disks/sped/ssc_sum/sum150_fix.tif", overwrite=TRUE, wopt=list(datatype='INT1U', names = c('sand', 'silt', 'clay')))

endCluster()
