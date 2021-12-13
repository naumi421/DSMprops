### Script to check RFE selection

## Folder with rfe files
rfefoldr <- "/nvme1/HYBconus100m/Gypsum_gRPI"

## RFE files
rfefiles <- list.files(path = rfefoldr,pattern="rf.RFE",full.names = T, recursive = F)

rfe0 <-readRDS(rfefiles[grepl("_0_cm",rfefiles)])
idx_rfe0 <- which(rfe0$rsq > 0.98*max(rfe0$rsq))
idx_rfe0 <- rfe0[min(idx_rfe0),]$var_idx

rfe5 <- readRDS(rfefiles[grepl("_5_cm",rfefiles)])
idx_rfe5 <- which(rfe5$rsq > 0.98*max(rfe5$rsq))
idx_rfe5 <- rfe5[min(idx_rfe5),]$var_idx

rfe15 <- readRDS(rfefiles[grepl("_15_cm",rfefiles)])
idx_rfe15 <- which(rfe15$rsq > 0.98*max(rfe15$rsq))
idx_rfe15 <- rfe15[min(idx_rfe15),]$var_idx

rfe30 <- readRDS(rfefiles[grepl("_30_cm",rfefiles)])
idx_rfe30 <- which(rfe30$rsq > 0.98*max(rfe30$rsq))
idx_rfe30 <- rfe30[min(idx_rfe30),]$var_idx

rfe60 <- readRDS(rfefiles[grepl("_60_cm",rfefiles)])
idx_rfe60 <- which(rfe60$rsq > 0.98*max(rfe60$rsq))
idx_rfe60 <- rfe60[min(idx_rfe60),]$var_idx

rfe100 <- readRDS(rfefiles[grepl("_100_cm",rfefiles)])
idx_rfe100 <- which(rfe100$rsq > 0.98*max(rfe100$rsq))
idx_rfe100 <- rfe100[min(idx_rfe100),]$var_idx

rfe150 <- readRDS(rfefiles[grepl("_150_cm",rfefiles)])
idx_rfe150 <- which(rfe150$rsq > 0.98*max(rfe150$rsq))
idx_rfe150 <- rfe150[min(idx_rfe150),]$var_idx

