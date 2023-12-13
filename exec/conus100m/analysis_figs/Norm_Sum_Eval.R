### Script to evaluate new normalized sand-silt-clay fraction estimates

# Workspace setup
# Install packages if not already installed
required.packages <- c("terra","sp","sf","ggplot2", "doParallel","parallel","mgcv")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)

## Folders
datafldr <- "/mnt/disks/sped/ssc_sum"
clayfldr <- "/mnt/covs/solus_preds/v2tst_gnat_pts/Clay_gRPI_250k"
sandfldr <- "/mnt/covs/solus_preds/v2tst_gnat_pts/Sand_gRPI_250k"
siltfldr <- "/mnt/covs/solus_preds/v2tst_gnat_pts/Silt_gRPI_250k"
newclayfldr <- "/mnt/solus100/Predictionsv2/trimmed/Clay_gRPI_250k/normalized"
newsandfldr <- "/mnt/solus100/Predictionsv2/trimmed/Sand_gRPI_250k/normalized"
newsiltfldr <- "/mnt/solus100/Predictionsv2/trimmed/Silt_gRPI_250k/normalized"

## Sum fix maps
fixfiles <- list.files(path = datafldr,pattern="_fix.tif",full.names = T, recursive = F)
sumfiles <- list.files(path = datafldr,pattern=".tif",full.names = T, recursive = F)
sumfiles <- sumfiles[!grepl("_fix",sumfiles)]
sumfiles <- sumfiles[!grepl("1to1plot",sumfiles)]

## Clay loop
clayfiles <- list.files(path = clayfldr,pattern="AvailPTS",full.names = T, recursive = F)
depths <- c("_0_","_5_","_15_","_30_","_60_","_100_","_150_")
#for(d in depths){
sumeval_fn <- function(d){
  pts <- readRDS(clayfiles[grepl(d,clayfiles)])
  pts <- pts[pts$tid == 'scd',]
  pts_sf <- st_as_sf(pts)
  pts_spatvect <- vect(pts_sf)
  r_fix <- rast(fixfiles[grepl(paste0("sum",gsub("_","",d),"_"),fixfiles)])[[c("clay")]]
  writeRaster(r_fix, filename = paste0(newclayfldr,"/claytotal",d,"cm_p.tif"),  datatype='INT1U')
  r_sum <- rast(sumfiles[grepl(paste0("sum",gsub("_","",d),".tif"),sumfiles)])
  pts_spatvect <- terra::extract(r_sum, pts_spatvect, bind=T)
  pts_spatvect <- pts_spatvect[pts_spatvect$sum < 80,]
  pts_spatvect <- terra::extract(r_fix, pts_spatvect, bind=T)
  model <- paste0("clay ",gsub("_","",d), " cm")
  viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
  scaleFUN <- function(x) round(x,0)
  plot1to1 <- ggplot(data=as.data.frame(pts_spatvect), aes(prop, clay)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1) + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("Measured") + ylab("Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste0("Sum fix ", "clay ", gsub("_","",d), " cm"))
  ggsave(paste(datafldr,'/1to1plot_',"clay",'_',gsub("_","",d),'_cm.tif',sep=""), plot = plot1to1, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  Rsq <- 1-var(pts_spatvect$prop_t - pts_spatvect$clay, na.rm=TRUE)/var(pts_spatvect$prop_t, na.rm=TRUE)
  RMSE <- sqrt(mean((pts_spatvect$prop - pts_spatvect$clay)^2, na.rm=TRUE))
  n <- nrow(as.data.frame(pts_spatvect))
  Bias <- mean(pts_spatvect$prop - pts_spatvect$clay, na.rm=TRUE)/mean(pts_spatvect$prop, na.rm=T)
  newrow <- data.frame(model = model, Rsq = Rsq, RMSE = RMSE, n=n, Bias = Bias)
  return(newrow)
}

## Linux parallel list apply
cpus <- 7
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
cvsum_fn.lst <- parLapply(cl,depths,try(sumeval_fn))
stopCluster(cl)
cv_sum_tab <- plyr::rbind.fill(cvsum_fn.lst)


## Sand loop
sandfiles <- list.files(path = sandfldr,pattern="AvailPTS",full.names = T, recursive = F)
depths <- c("_0_","_5_","_15_","_30_","_60_","_100_","_150_")
#for(d in depths){
sumeval_fn2 <- function(d){
  pts <- readRDS(sandfiles[grepl(d,sandfiles)])
  pts <- pts[pts$tid == 'scd',]
  pts_sf <- st_as_sf(pts)
  pts_spatvect <- vect(pts_sf)
  r_fix <- rast(fixfiles[grepl(paste0("sum",gsub("_","",d),"_"),fixfiles)])[[c("sand")]]
  writeRaster(r_fix, filename = paste0(newsandfldr,"/sandtotal",d,"cm_p.tif"),  datatype='INT1U')
  r_sum <- rast(sumfiles[grepl(paste0("sum",gsub("_","",d),".tif"),sumfiles)])
  pts_spatvect <- terra::extract(r_sum, pts_spatvect, bind=T)
  pts_spatvect <- pts_spatvect[pts_spatvect$sum < 80,]
  pts_spatvect <- terra::extract(r_fix, pts_spatvect, bind=T)
  model <- paste0("sand ",gsub("_","",d), " cm")
  viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
  scaleFUN <- function(x) round(x,0)
  plot1to1 <- ggplot(data=as.data.frame(pts_spatvect[pts_spatvect$sum < 80,]), aes(prop, sand)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1) + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("Measured") + ylab("Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste0("Sum fix ", "sand ", gsub("_","",d), " cm"))
  plot1to1
  ggsave(paste(datafldr,'/1to1plot_',"sand",'_',gsub("_","",d),'_cm.tif',sep=""), plot = plot1to1, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  Rsq <- 1-var(pts_spatvect$prop_t - pts_spatvect$sand, na.rm=TRUE)/var(pts_spatvect$prop_t, na.rm=TRUE)
  RMSE <- sqrt(mean((pts_spatvect$prop - pts_spatvect$sand)^2, na.rm=TRUE))
  n <- nrow(as.data.frame(pts_spatvect))
  Bias <- mean(pts_spatvect$prop - pts_spatvect$sand, na.rm=TRUE)/mean(pts_spatvect$prop, na.rm=T)
  newrow <- data.frame(model = model, Rsq = Rsq, RMSE = RMSE, n=n, Bias = Bias)
  return(newrow)
}

## Linux parallel list apply
cpus <- 7
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
cvsum_fn.lst2 <- parLapply(cl,depths,try(sumeval_fn2))
stopCluster(cl)
cv_sum_tab2 <- plyr::rbind.fill(cvsum_fn.lst2)


## Silt loop
siltfiles <- list.files(path = siltfldr,pattern="AvailPTS",full.names = T, recursive = F)
depths <- c("_0_","_5_","_15_","_30_","_60_","_100_","_150_")
#for(d in depths){
sumeval_fn3 <- function(d){
  pts <- readRDS(siltfiles[grepl(d,siltfiles)])
  pts <- pts[pts$tid == 'scd',]
  pts_sf <- st_as_sf(pts)
  pts_spatvect <- vect(pts_sf)
  r_fix <- rast(fixfiles[grepl(paste0("sum",gsub("_","",d),"_"),fixfiles)])[[c("silt")]]
  writeRaster(r_fix, filename = paste0(newsiltfldr,"/silttotal",d,"cm_p.tif"),  datatype='INT1U')
  r_sum <- rast(sumfiles[grepl(paste0("sum",gsub("_","",d),".tif"),sumfiles)])
  pts_spatvect <- terra::extract(r_sum, pts_spatvect, bind=T)
  pts_spatvect <- pts_spatvect[pts_spatvect$sum < 80,]
  pts_spatvect <- terra::extract(r_fix, pts_spatvect, bind=T)
  model <- paste0("silt ",gsub("_","",d), " cm")
  viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
  scaleFUN <- function(x) round(x,0)
  plot1to1 <- ggplot(data=as.data.frame(pts_spatvect), aes(prop, silt)) +
    stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1) + #xlim(-5,105) + ylim(-5,105) +
    theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
    xlab("Measured") + ylab("Prediction") + scale_fill_gradientn(name = "Count", trans = "log", colours = rev(viri),labels=scaleFUN) +
    ggtitle(paste0("Sum fix ", "silt ", gsub("_","",d), " cm"))
  ggsave(paste(datafldr,'/1to1plot_',"silt",'_',gsub("_","",d),'_cm.tif',sep=""), plot = plot1to1, device = "tiff", dpi = 600, limitsize = TRUE, width = 6, height = 5, units = 'in',compression = c("lzw"))
  Rsq <- 1-var(pts_spatvect$prop_t - pts_spatvect$silt, na.rm=TRUE)/var(pts_spatvect$prop_t, na.rm=TRUE)
  RMSE <- sqrt(mean((pts_spatvect$prop - pts_spatvect$silt)^2, na.rm=TRUE))
  n <- nrow(as.data.frame(pts_spatvect))
  Bias <- mean(pts_spatvect$prop - pts_spatvect$silt, na.rm=TRUE)/mean(pts_spatvect$prop, na.rm=T)
  newrow <- data.frame(model = model, Rsq = Rsq, RMSE = RMSE, n=n, Bias = Bias)
  return(newrow)
}

## Linux parallel list apply
cpus <- 7
cl <- makeCluster(cpus, type="FORK")
registerDoParallel(cl)
cvsum_fn.lst3 <- parLapply(cl,depths,try(sumeval_fn3))
stopCluster(cl)
cv_sum_tab3 <- plyr::rbind.fill(cvsum_fn.lst3)

## Final Sum Evaluation table
SumFix_Eval_Tab <- rbind(cv_sum_tab,cv_sum_tab2,cv_sum_tab3)
## Now Save Table
write.csv(SumFix_Eval_Tab,paste0(datafldr,"/texture_sum_fix_eval_",gsub("-","",Sys.Date()),".csv"))

## Summarize Sum maps: 1 million random points
randpts <- readRDS(paste("/mnt/solus100/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final","/CONUS_random_gRPIsamp.rds",sep=""))
randpts_sf <- st_as_sf(randpts)
randpts_spatvect <- vect(randpts_sf)
for(f in sumfiles){
  r_ext <- rast(f)
  randpts_spatvect <- terra::extract(r_ext, randpts_spatvect, bind=T)
}

randpts_df <- as.data.frame(randpts_spatvect)
summary(randpts_df)
