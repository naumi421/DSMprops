
# Workspace setup
# Install packages if not already installed
required.packages <- c("raster", "doParallel" ,"dplyr","DSMprops")
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
rastfldr <- "/mnt/solus100/Predictionsv2"
newrastfldr <- "/mnt/solus100/Predictionsv2/trimmed"

## Create lists of folders and files to cycle through
# Files
topfiles <- list.files(path = topfolder,pattern="TrainPTS",full.names = T, recursive = T)
topfiles <- topfiles[grepl("_250k",topfiles)]
topfiles <- topfiles[!grepl("first_test",topfiles)]
topfiles <- topfiles[!grepl("compare",topfiles)]
topfiles <- topfiles[!grepl("analysis",topfiles)]
topfiles <- topfiles[!grepl("models_more",topfiles)]
topfiles <- topfiles[!grepl("results_oldSCD", topfiles)]
topfiles <- topfiles[!grepl("ValPlot", topfiles)]
topfiles <- topfiles[!grepl("ptdensity", topfiles)]
topfiles <- topfiles[!grepl("masked", topfiles)]
topfiles <- topfiles[!grepl("AWC", topfiles)]
# Directories
topdirs <-  list.dirs(topfolder,recursive = F, full.names = T)
topdirs
## Subsetting
topdirs <- topdirs[grepl("_250k",topdirs)]
topdirs <- topdirs[!grepl("first_test",topdirs)]
topdirs <- topdirs[!grepl("analysis",topdirs)]
topdirs <- topdirs[!grepl("models_more",topdirs)]
topdirs <- topdirs[!grepl("AWC",topdirs)]
topdirsnm <- basename(topdirs)

## Untrimmed Prediction folders
rastdirs <-  list.dirs(rastfldr,recursive = F, full.names = T)
rastdirs <- rastdirs[!grepl("trimmed",rastdirs)]

#Now create new subdirectories for all folders for new files to be stored in.
# for(f in topdirsnm){
#   newdir <- paste0(newrastfldr,"/",f)
#   dir.create(newdir)
# }

## Deal with restriction and lithic (not depth specific)
restrfiles <- topfiles[grepl("RestrDepth",topfiles)]
topfiles<- topfiles[!grepl("RestrDepth",topfiles)]
mskfn_d <- function(rast){
  ind <- rast
  ind[ind<0]<-0 # to bring the slighly negative predictions back to zero
  ind[ind>201]<-201
  return(ind)
}
restrast <- raster("/mnt/covs/solus_preds/v2tst_gnat_pts/RestrDepth_gRPI_250k/masked/resdept_r_1x_all_cm_2D_QRFadj.tif")
lithrast <- raster("/mnt/covs/solus_preds/v2tst_gnat_pts/RestrDepth_gRPI_250k/masked/anylithicdpt_1x_all_cm_2D_QRFadj.tif")
newrestrast <- calc(restrast,fun=mskfn_d,filename="/mnt/solus100/Predictionsv2/trimmed/RestrDepth_gRPI_250k/resdept_r_1x_all_cm_2D_QRFadj.tif",
                    options=c("COMPRESS=DEFLATE", "TFW=YES"), datatype = dataType(restrast),progress = "text", overwrite = T)
newlithrast <- calc(lithrast,fun=mskfn_d,filename="/mnt/solus100/Predictionsv2/trimmed/RestrDepth_gRPI_250k/anylithicdpt_1x_all_cm_2D_QRFadj.tif",
                    options=c("COMPRESS=DEFLATE", "TFW=YES"), datatype = dataType(lithrast),progress = "text", overwrite = T)

## Check for completed files (for rerunning after errors)
compfiles <- list.files(path = newrastfldr,pattern="QRFadj",full.names = T, recursive = T)
compfiles <- compfiles[!grepl(".tfw",compfiles)]

## Index for list apply
rlist <- 1:length(topfiles)

# ## par list apply fn
trimmask_fn <- function(g){
  folder <- dirname(topfiles[g])
  filepath <- topfiles[g]
  trnpts <- readRDS(filepath)
  filenm <- basename(filepath)
  dirnm <- basename(folder)
  predpath <-  rastdirs[grepl(dirnm, rastdirs)]
  for(pth in predpath){
  if(basename(pth) == dirnm){
    predpath <- pth
  }
  }
  # Bring in raster to trim
  depth <- unlist(strsplit(filenm, "_"))[4]
  depth_str <- paste0(depth, '_cm')
  rastfile <- list.files(path = predpath,pattern=paste0("_",depth_str, "_2D_QRFadj"),full.names = T, recursive = T)
  rastfile <- rastfile[!grepl(".tfw",rastfile)]
  rastfilenm <- basename(rastfile)
  rast <- raster(rastfile)
  names(rast) <- "rast"
  newfolder <- paste0(newrastfldr,"/",dirnm)
  newfile <- paste0(newfolder, "/",rastfilenm)
  if(!newfile %in% compfiles){
    # Now determine raster scalar and min/max values to trim to
    minval <- min(trnpts$prop, na.rm = T)
    maxval <- max(trnpts$prop, na.rm = T)
    if(grepl('1x', rastfilenm)){
      scalar <- 1
    }
    if(grepl('10x', rastfilenm)){
      scalar <- 10
    }
    if(grepl('100x', rastfilenm)){
      scalar <- 100
    }
    if(grepl('1000x', rastfilenm)){
      scalar <- 1000
    }
    # If statement to only create if file does not already exist.
    mskfn <- function(rast){
      ind <- rast/scalar
      ind[ind<minval]<-minval
      ind[ind>maxval]<-maxval
      ind <- ind*scalar
      return(ind)
    }
    calc(rast,fun=mskfn,filename=newfile, options=c("COMPRESS=DEFLATE", "TFW=YES"), datatype = dataType(rast),progress = "text")
  }

  gc()
}


## Linux parallel list apply
cpus <- min(length(rlist),detectCores() - 2)
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
parLapply(cl,rlist,try(trimmask_fn))
stopCluster(cl)


## SOC no data corrections: 10/26/2023
socdir <- topdirs[grepl("SOC_",topdirs)]
socfiles <- list.files(path = socdir,pattern="QRF_bt.tif",full.names = T, recursive = T)
socfiles <- socfiles[grepl("masked",socfiles)]
soctrimfiles <- list.files(path = "/mnt/solus100/Predictionsv2/trimmed/SOC_gRPI_250k",pattern="QRFadj_bt.tif",full.names = T, recursive = T)
soctrainfiles <- list.files(path = socdir,pattern="TrainPTS",full.names = T, recursive = T)

## index for appl
#socind <- 1:length(socfiles)

#soc_fix_fn <- function(si){
for(sdep in c("_0_","_5_","_15_","_30_","_60_","_100_","_150_")){
  soc_orig_file <- socfiles[grepl(sdep,socfiles)]
  soc_orig_file <- raster(soc_orig_file)
  soc_trim_file <- soctrimfiles[grepl(sdep,soctrimfiles)]
  soc_trim_file_new <- paste0(gsub(".tif","",soc_trim_file),"_new.tif")
  soc_trim_file_r <- raster(soc_trim_file)
  trainfile <- soctrainfiles[grepl(sdep,soctrainfiles)]
  trainfile <-  readRDS(trainfile)
  socmax <- max(trainfile$prop)*1000
  newssoctrim_fn <- function(x,y){
    ind <- ifelse(x > -1 & is.na(y), socmax, y)
    return(ind)
  }
  socstk <- stack(soc_orig_file,soc_trim_file_r)
  newsocr <- overlay(socstk,fun=newssoctrim_fn,filename=soc_trim_file_new, options=c("COMPRESS=DEFLATE", "TFW=YES"), datatype = dataType(soc_orig_file),progress = "text", overwrite=TRUE)
}






