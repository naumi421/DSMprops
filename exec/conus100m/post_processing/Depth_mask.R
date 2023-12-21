### Script to implement final masks and rename all SOLUS files to match standards

# Workspace setup
# Install packages if not already installed
required.packages <- c("raster","sp","sf","ggplot2", "doParallel","parallel","mgcv")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
rasterOptions(maxmemory = 2e+08,chunksize = 4e+07, memfrac = 0.1)

## Folders
oldrastfldr <- "/mnt/disks/sped/solus100preds/public_layers"
newfldr <- "/mnt/disks/sped/solus100preds/public_layers_depthmasked"


## Other prediction files
lyrfiles <- list.files(path = oldrastfldr,pattern=".tif",full.names = T, recursive = T)
lyrfiles <- lyrfiles[!grepl("lithic",lyrfiles)]
lyrfiles <- lyrfiles[!grepl("resdept_",lyrfiles)]

## Final mask raster
r_msk <- raster("/mnt/disks/sped/solus100preds/public_layers/anylithicdpt_cm_2D_p.tif")


### Now process the rest of the files
lyrind <- seq(1:length(lyrfiles))

## General file loop
layer_fn <- function(f){
  fname <- lyrfiles[f]
  r <- raster(fname)
  dtype <- dataType(r)
  fbase <- basename(fname)
  depth <-  as.numeric(unlist(strsplit(fbase, "_"))[2])
  r[r_msk < depth] <- NA
  writeRaster(r, filename = paste0(newfldr,"/",fbase),  options=c("COMPRESS=DEFLATE", "TFW=NO"), datatype = dtype)
  gc()
}

## Linux parallel list apply
cpus <- min(120,detectCores() - 2,length(lyrind))
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
parLapply(cl,lyrind,try(layer_fn))
stopCluster(cl)


