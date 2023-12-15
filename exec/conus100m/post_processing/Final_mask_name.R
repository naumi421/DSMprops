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
datafldr <- "/mnt/disks/sped/solus100preds"
claypredfldr <- "/mnt/disks/sped/solus100preds/Predictionsv2/trimmed/Clay_gRPI_250k/normalized"
sandpredfldr <- "/mnt/disks/sped/solus100preds/Predictionsv2/trimmed/Sand_gRPI_250k/normalized"
siltpredfldr <- "/mnt/disks/sped/solus100preds/Predictionsv2/trimmed/Silt_gRPI_250k/normalized"
oldpredfldr <- "/mnt/disks/sped/solus100preds/Predictionsv2/trimmed"
oldrastfldr <- "/mnt/disks/sped/solus100preds/Predictionsv2"
newfldr <- "/mnt/disks/sped/solus100preds/public_layers"

## Normalized sand/silt/clay prediction files
textnfiles <- list.files(path = claypredfldr,pattern=".tif",full.names = T, recursive = F)
textnfiles <- append(textnfiles, list.files(path = sandpredfldr,pattern=".tif",full.names = T, recursive = F))
textnfiles <- append(textnfiles, list.files(path = siltpredfldr,pattern=".tif",full.names = T, recursive = F))

## Other prediction files
predfiles <- list.files(path = oldpredfldr,pattern=".tif",full.names = T, recursive = T)
predfiles <- predfiles[!grepl("normalized",predfiles)]
predfiles <- predfiles[!grepl("ssc_sum",predfiles)]
predfiles <- predfiles[!grepl("claytotal",predfiles)]
predfiles <- predfiles[!grepl("sandtotal",predfiles)]
predfiles <- predfiles[!grepl("silttotal",predfiles)]

## Uncertainty files
uncfiles <- list.files(path = oldrastfldr,pattern=".tif",full.names = T, recursive = T)
uncfiles <- uncfiles[!grepl("trimmed",uncfiles)]
uncfiles <- uncfiles[!grepl("QRFadj",uncfiles)]
uncfiles <- uncfiles[!grepl("QRF.tif",uncfiles)]
uncfiles <- uncfiles[!grepl("QRF_bt.tif",uncfiles)]

## Final mask raster
r_msk <- raster("/mnt/disks/sped/solus100preds/notpub_conus_mask.tif")

## Texture index
textind <- seq(1:length(textnfiles))

## Texture loop
texture_fn <- function(f){
  fname <- textnfiles[f]
  fbase <- basename(fname)
  prop <- unlist(strsplit(fbase, "_"))[1]
  depth <- paste(unlist(strsplit(fbase, "_"))[2],unlist(strsplit(fbase, "_"))[3],sep = "_")
  scalar <- 1 ## For scaled integers
  ftype <- "prediction"
  r <- raster(fname)
  r_out <- r*r_msk
  writeRaster(r_out, filename = paste0(newfldr,"/",fbase),  options=c("COMPRESS=DEFLATE", "TFW=NO"), datatype = dataType(r))
  newrow <- data.frame(filename = fbase, property = prop, depth = depth, filetype = ftype, scalar = scalar)
  return(newrow)
}

## Linux parallel list apply
cpus <- min(124,detectCores() - 2,length(textind))
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
lyr_fn.lst <- parLapply(cl,textind,try(texture_fn))
stopCluster(cl)
lyr_sum_tab <- plyr::rbind.fill(lyr_fn.lst)


### Now process the rest of the files
lyrfiles <- append(predfiles, uncfiles)
lyrind <- seq(1:length(lyrfiles))

## General file loop
layer_fn <- function(f){
  fname <- lyrfiles[f]
  r <- raster(fname)
  fbase <- basename(fname)
  prop <- unlist(strsplit(fbase, "_"))[1]
  depth <- paste(unlist(strsplit(fbase, "_"))[4],unlist(strsplit(fbase, "_"))[5],sep = "_")
  scalar <- gsub("x", "", unlist(strsplit(fbase, "_"))[3])
  if(scalar == 1000 & grepl("relwidth",fbase)){
    r <- r/10
    scalar <- 100
  }
  if(grepl("QRFadj",fbase)){ftype <- "prediction"
  newflnm <- paste0(prop,"_",depth,"_p.tif")}
  if(grepl("95PI_l",fbase)){ftype <- "95% low prediction interval"
  newflnm <- paste0(prop,"_",depth,"_l.tif")}
  if(grepl("95PI_h",fbase)){ftype <- "95% high prediction interval"
  newflnm <- paste0(prop,"_",depth,"_h.tif")}
  if(grepl("relwidth",fbase)){ftype <- "relative prediction interval"
  newflnm <- paste0(prop,"_",depth,"_rpi.tif")}
  r_out <- r*r_msk
  writeRaster(r_out, filename = paste0(newfldr,"/",newflnm),  options=c("COMPRESS=DEFLATE", "TFW=NO"), datatype = dataType(r))
  newrow <- data.frame(filename = newflnm, property = prop, depth = depth, filetype = ftype, scalar = scalar)
  return(newrow)
  gc()
}

## Linux parallel list apply
cpus <- min(120,detectCores() - 2,length(lyrind))
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
lyr_fn.lst2 <- parLapply(cl,lyrind,try(layer_fn))
stopCluster(cl)
lyr_sum_tab2 <- plyr::rbind.fill(lyr_fn.lst2)


## Final layer table
Final_Layer_Tab <- rbind(lyr_sum_tab,lyr_sum_tab2)
## Now Save Table
write.csv(Final_Layer_Tab,paste0(newfldr,"/Final_Layer_Table_",gsub("-","",Sys.Date()),".csv"))


