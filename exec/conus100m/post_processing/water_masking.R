
# Workspace setup
# Install packages if not already installed
required.packages <- c("raster","rgdal", "rasterVis","maptools","RColorBrewer","ggplot2","gridExtra","classInt","RStoolbox","hexbin", "ranger", "parallel", "doParallel" ,"dplyr","Hmisc","viridisLite","DSMprops")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
# memory.limit(100000)
rasterOptions(maxmemory = 6.5e+09,chunksize = 5e+08, memfrac = 0.9)

cpus <- min(124,detectCores() - 2)

## Key Folder Locations
topfolder <- "/mnt/covs/solus_preds/v2tst_gnat_pts"

## Create lists of folders and files to cycle through
# Files
topfiles <- list.files(path = topfolder,pattern=".tif$",full.names = T, recursive = T)
topfiles <- topfiles[grepl("_250k",topfiles)]
topfiles <- topfiles[!grepl("first_test",topfiles)]
topfiles <- topfiles[!grepl("compare",topfiles)]
topfiles <- topfiles[!grepl("analysis",topfiles)]
topfiles <- topfiles[!grepl("models_more",topfiles)]
topfiles <- topfiles[!grepl("results_oldSCD", topfiles)]
topfiles <- topfiles[!grepl("ValPlot", topfiles)]
topfiles <- topfiles[!grepl("ptdensity", topfiles)]
topfiles <- topfiles[!grepl("masked", topfiles)]
# Directories
topdirs <-  list.dirs(topfolder,recursive = F, full.names = T)
topdirs
## Subsetting
topdirs <- topdirs[grepl("_250k",topdirs)]
topdirs <- topdirs[!grepl("first_test",topdirs)]
topdirs <- topdirs[!grepl("analysis",topdirs)]
topdirs <- topdirs[!grepl("models_more",topdirs)]
topdirsnm <- basename(topdirs)

## Now create new subdirectories for all folders for new files to be stored in.
# for(f in topdirs){
#   newdir <- paste0(f,"/masked")
#   dir.create(newdir)
# }


## Index for list apply
rlist <- 1:length(topfiles)

## Raster masking function
mskfn <- function(rast,mask){
    ind <- rast*mask
    ind[ind<0]<-0 # to bring the slighly negative predictions back to zero
    return(ind)
  }


# ## par list apply fn
watermask_fn <- function(g){
  folder <- dirname(topfiles[g])
  filepath <- topfiles[g]
  rast <- raster(filepath)
  filenm <- basename(filepath)
  names(rast) <- "rast"
  newfolder <- paste0(folder,"/masked")
  existingfiles <- list.files(path = newfolder,pattern=".tif$",full.names = T, recursive = F)
  newfile <- paste0(folder,"/masked/",filenm)
  # If statement to only create if file does not already exist.
  if(!newfile %in% existingfiles){
    ## NLCD mask raster
    mask <- raster("/mnt/covs/nlcd19/nlcd19_water_ice_solus100.tif")
    h2ostk <- stack(rast,mask)
    overlay(h2ostk,fun=mskfn,filename=newfile, options=c("COMPRESS=DEFLATE", "TFW=YES"), datatype = dataType(rast),progress = "text")
  }
  gc()
}


## Linux parallel list apply
cpus <- min(length(rlist),detectCores() - 2)
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
parLapply(cl,rlist,try(watermask_fn))
stopCluster(cl)













