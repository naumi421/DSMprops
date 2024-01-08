
# Workspace setup
# Install packages if not already installed
required.packages <- c("raster", "doParallel","gdalUtilities") #"gdalraster"
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
# memory.limit(100000)
rasterOptions(maxmemory = 6.5e+09,chunksize = 5e+08, memfrac = 0.9)

cpus <- min(124,detectCores() - 2)

## Folders
infldr <- "/mnt/disks/sped/solus100preds/notpublic_layers"
outfldr <- "/mnt/disks/sped/solus100preds/notpublic_layers_cog"

## Files and index
infiles <- list.files(path = infldr,pattern=".tif*",full.names = T, recursive = T)
rlist <- 1:length(infiles)

cog_fn <- function(r){
  rast <- raster(infiles[r])
  fullpath <- infiles[r]
  #src_dataset <- system.file(fullpath, package="gdalUtils")
  flnm <- basename(infiles[r])
  dtype <- dataType(rast)
  if(dtype == "INT1U"){dtypegdal <- "Byte"}
  if(dtype == "INT2U"){dtypegdal <- "UInt16"}
  if(dtype == "FLT4S"){dtypegdal <- "Float32"}
  # args <- c("-ot", dtypegdal)
  # args <- c(args, "-of", "COG", "-co", "COMPRESSED=YES")
  # translate(fullpath, file.path(outfldr,flnm), args)
  gdal_translate(fullpath,file.path(outfldr,flnm),of="COG", ot = dtypegdal)
}

## Linux parallel list apply
cpus <- min(length(rlist),detectCores() - 2)
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
parLapply(cl,rlist,try(cog_fn))
stopCluster(cl)
