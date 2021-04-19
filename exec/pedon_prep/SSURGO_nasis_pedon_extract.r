###################################
### Script to extract property dat from gSSURGO
### to NASIS series transect observations

## Packages
required.packages <- c("raster", "sp", "rgdal","doParallel", "plyr", "dplyr", "ncdf4","maptools", "rgeos","stats","spdep","sf","aqp")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
##Raster settings
rasterOptions(maxmemory = 2e+09, chunksize = 2e+08)

## Key folders
resultfolder <- "/push/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final"
ptsfolder <- "/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map/2020_Pedons"
ssurgo_fgdb <- "/push/gSSURGO20/gSSURGO_CONUS.gdb"
ssurgofolder <- "/push/gSSURGO20"
ssurgohuc_fldr <- "/push/gSSURGO20/ssurgo_huc"
loopresults <- "/push/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final/loopresults"
hucfolder <- "/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map"

### Explore gSSURGO GDB and open necessary files
# Vignettes to read all GDB files
fc_list <- ogrListLayers(ssurgo_fgdb) # list layers
polys <- readOGR(ssurgo_fgdb, "MUPOLYGON") # Takes 24 hrs to load with SSURGO for all USA
saveRDS(polys,paste(ssurgofolder,"/gSSURGO20_conus_mupolys.rds",sep=""))
polys <- readRDS(paste(ssurgofolder,"/gSSURGO20_conus_mupolys.rds",sep=""))
polys@data[] <- lapply(polys@data, function(x) if (is.factor(x)) as.character(x) else {x}) # factors to char
polys.proj <- projection(polys)

### Break SSURGO up into chunks to process ###########
## Load map clip boundary (if needed) 
polybounds <- readOGR(hucfolder, "Huc4_conus")
polybounds@data[] <- lapply(polybounds@data, function(x) if (is.factor(x)) as.character(x) else {x}) # factors to char
polybounds <- spTransform(polybounds, polys.proj)
for(h in 1:nrow(polybounds@data)){
  polybound <- polybounds[h,]
  ssurgochnk <- polys[polybound,] # Clipping down to bite size chunks (huc4s)
  saveRDS(ssurgochnk,paste(ssurgohuc_fldr,"/hucnum_",h,"_gSSURGO20_mupolys.rds",sep=""))
  print(paste("huc",h, "was done at", Sys.time(), sep=" "))
}

## List of unique mukeys in area
mukeys.df <- sf::st_read(dsn = ssurgo_fgdb, layer = "mapunit")
mukeys.df[] <- lapply(mukeys.df, function(x) if (is.factor(x)) as.character(x) else {x})
saveRDS(mukeys.df, paste(ssurgofolder,"/mukeysdf_conus.rds", sep=""))
mukeys.df <- readRDS(paste(ssurgofolder,"/mukeysdf_conus.rds", sep=""))
mukeys <- unique(mukeys.df$mukey)

## Also load necessary data tables and subset to study area
## Component table
comp.df <- sf::st_read(dsn = ssurgo_fgdb, layer = "component")
comp.df[] <- lapply(comp.df, function(x) if (is.factor(x)) as.character(x) else {x})
comps <- comp.df[comp.df$mukey %in% mukeys,]
comps$compname <- toupper(comps$compname)
# Clean up compnames for matching
comps$compname = gsub(" VARIANT","", comps$compname)
comps$compname = gsub(" TAXADJUNCT","", comps$compname)
comps$compname = gsub(" TAXAJUNCT","", comps$compname)
comps$compname = gsub(" FAMILY","", comps$compname)
comps$compname = gsub(" LIKE","", comps$compname)
comps$compname = gsub("-LIKE","", comps$compname)
comps$compname = gsub("-SIMILAR","", comps$compname)
comps$compname = gsub(" SIMILAR","", comps$compname)
comps$compname = gsub(" TAX.","", comps$compname)
comps$compname = gsub(" TAX","", comps$compname)
comps$compname = gsub(" ERODED","", comps$compname)
comps$compname = gsub(", ERODED","", comps$compname)
saveRDS(comps, paste(ssurgofolder,"/comps_conus.rds", sep=""))
comps <- readRDS(paste(ssurgofolder,"/comps_conus.rds", sep=""))
cokeys <- unique(comps$cokey)


#### Load in NASIS points
## Soil Profile Collection Import: just site table in this case
load(paste(ptsfolder,"/nasis_sites_20210325.RData",sep="")) # all nasis obs incl just site taxa: s
## Load offcial component-linked pedon that does not need to be linked
linked <- read.csv("/home/tnaum/OneDrive/USGS/NCSS/DSM_Focus_team/Natl_map/2020_Pedons/copedon-link-cokey.csv", stringsAsFactors = F)
excl_lst <- which(!as.character(s$peiid) %in% as.character(unique(linked$peiid))) # index of sites not included in linked pedons
s.ids <- s$peiid[excl_lst] # IDs to select just pedons not included in the linked pedon table
# ## Organize data for extraction
s.sites <- s[,c("peiid","taxonname_recent","x_std","y_std")]
shp.pts.all <- na.omit(s.sites)
coordinates(shp.pts.all) <- ~ x_std + y_std
# TODO Need to update to proj6
proj4string(shp.pts.all) <- CRS("+init=epsg:4326") 
## Ssave as gkg
outname <- paste(ptsfolder,"/nasis_pedons_20210325_spatial.gpkg",sep="")
outname.shp <- paste(ptsfolder,"/",sep="")
writeOGR(shp.pts.all,dsn=outname,layer="nasis_pedons_20210325_spatial.gpkg",driver="GPKG",overwrite_layer = T)
writeOGR(shp.pts.all,dsn=outname.shp,layer="nasis_pedons_20210325_shp",driver="ESRI Shapefile",overwrite_layer = T)
## Select pts to use in spatial extrations
shp.pts <- shp.pts.all[shp.pts.all$peiid %in% s.ids,] # subset to non matched pedons
## Clear out some objects to keep memory clear
rm(s,s.sites)
gc()

## Clean up pedon fields and compnames for matching
shp.pts@data[] <- lapply(shp.pts@data, function(x) if (is.factor(x)) as.character(x) else {x}) # factors to char
shp.pts@data$taxonname <- toupper(shp.pts@data$taxonname_recent)
## Clean up pt compnames for matching
shp.pts@data$taxonname = gsub(" VARIANT","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(" TAXADJUNCT","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(" TAXAJUNCT","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(" FAMILY","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(" LIKE","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub("-LIKE","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub("-SIMILAR","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(" SIMILAR","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(" TAX.","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(" TAX","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(" ERODED","", shp.pts@data$taxonname)
shp.pts@data$taxonname = gsub(", ERODED","", shp.pts@data$taxonname)
## Project
polyex <- readRDS(paste(ssurgohuc_fldr,"/",list.files(ssurgohuc_fldr,pattern="rds$")[1],sep=""))
polys.proj <- projection(polyex)
rm(polyex)
gc()
shp.pts <- spTransform(shp.pts, polys.proj)

## Create list of observation point ids (pid) for all observations
shp.pts@data$pid <- seq_len(nrow(shp.pts@data))
pids <- shp.pts@data$pid

## Function for spatially matching pedons to SSURGO components by 1) selecting home and adjacent map units
## and 2) matching to a cleaned component name first checking home mu then adjacents.
spatial.extfn <- function(pid){ ### list apply function using pedon pid to parallelize
  obs <- hucshp.pts[hucshp.pts$pid %in% pid,]
  obs_buff <- gBuffer(obs,width=2000)
  obs_df <- as.data.frame(obs)
  obs_df$mtchtype <- NA
  polysbuf <- try(raster::intersect(obs_buff, polys))
  polyint <- try(sp::over(obs,polysbuf))
  comps1 <- try(comps[comps$mukey %in% polyint$MUKEY,])
  comp_mtch1 <- try(comps1[comps1$compname %in% obs_df$taxonname,][1,])
  # Check to see if point lined up with any polys
  if(is.data.frame(comp_mtch1)){
    # Now check for a match in home map unit using the comppct_r (always in SSURGO for populated components)
    if(!is.na(comp_mtch1$comppct_r)){
      obs_df$mtchtype <- "home"
      obs_df <- try(cbind(obs_df,comp_mtch1))
      try(return(obs_df))
      try(rm(obs_df,polyint,polysbuf,obs,obs_buff,comp_mtch1))
      #gc()
      # If no match in home map unit, check adjacent polygons within buffer
    } else {
      polyid <- try(polyint$polyid)
      poly1 <- try(polys[polys$polyid %in% polyid,])
      polyneigh_mtx <- try(gTouches(polysbuf, poly1, byid = TRUE, returnDense = FALSE))
      neigh.list <- c()
      for(i in seq(1:length(polyneigh_mtx))){
        neigh <- polyneigh_mtx[i]
        if(length(neigh[[1]])>0){
          neigh.list <- append(neigh.list,names(neigh), after = length(neigh.list))
          #print(i)
        }
      }
      polys2 <- try(polysbuf[rownames(polysbuf@data) %in% neigh.list,])
      mukeysadj <- try(c(polys2$MUKEY))
      comps_chk <- try(comps[comps$mukey %in% mukeysadj,]) 
      comp_mtch2 <- try(comps_chk[comps_chk$compname %in% obs_df$taxonname,][1,])
      # Now force a cbind and test to see if a match was made
      if(!is.na(comp_mtch2$comppct_r)){
      #if(ncol(obs_df)>5){
        obs_df$mtchtype <- "adjacent"
        obs_df <- try(cbind(obs_df,comp_mtch2))
        try(return(obs_df))
        try(rm(obs_df,polys2,comps_chk, polyneigh_mtx, poly1,mukeysadj,polyint,polyid,polysbuf,obs,obs_buff,comp_mtch2,comp_mtch1,neigh.list))
        gc()
      } else {
        try(rm(obs_df,polys2,comps_chk, polyneigh_mtx, poly1,mukeysadj,polyint,polyid,polysbuf,obs,obs_buff,comp_mtch1,neigh.list))
        #gc()
      }
    }
  } else {
    try(rm(obs_df,polyint,polysbuf,obs,obs_buff,comp_mtch1,comps1))
    #gc()
  }
}


## Parallel list apply
shucs <- list.files(path = ssurgohuc_fldr,pattern=".rds$", full.names = T, recursive = F)
Sys.time()
for(h in 1:length(shucs)){
  pretime <- Sys.time()
  polys <- readRDS(shucs[h])
  #hucshp.pts <- shp.pts[polys,] ## too slow
  hucshp.pts <- crop(shp.pts, extent(polys))
  hucshp.pts@data <- hucshp.pts@data[,c('pid','taxonname')]
  brks_list_for <- hucshp.pts@data$pid
  ## Poly IDs
  polys@data$polyid <- seq_len(nrow(polys@data))
  cpus <- detectCores() - 2
  cl <- makeCluster(cpus, type="FORK")
  registerDoParallel(cl)
  obs_df_list <- parLapply(cl,brks_list_for,try(spatial.extfn))
  stopCluster(cl)
  obs_all_df <- do.call("rbind", obs_df_list)
  #obs_all_dfc <- na.omit(obs_all_df)
  #setwd(loopresults)
  saveRDS(obs_all_df, paste(loopresults,"/obs_all_df_",h,".rds",sep=""))
  rm(obs_df_list,obs_all_df)
  gc()
  posttime <- Sys.time()
  runtime <- posttime - pretime
  print(paste(h, "was done at", posttime,"in",runtime, sep=" "))
}

## Loop through all result files and bind back together into one dataset
optfiles <- list.files(path = loopresults, pattern=".rds$", full.names = T, recursive = F)
obs_df_all <- readRDS(optfiles[1])
obs_df_all <- obs_df_all[!is.na(obs_df_all$comppct_r),]
for(t in 2:length(optfiles)){
  newtrial <- readRDS(optfiles[t])
  newtrial <- newtrial[!is.na(newtrial$comppct_r),]
  obs_df_all <- rbind(obs_df_all,newtrial)
  rm(newtrial)
  print(paste("Done with ",t))
}


##### Now join back to original pedon ids and put into soil profile collection ###
## Clean up and join spatially matched pedons
obs_df_links <- data.frame(pid=obs_df_all[!duplicated(obs_df_all$pid),c("pid")])
obs_df_links <- left_join(obs_df_links,shp.pts@data[,c("peiid","taxonname_recent","pid")],by="pid")
obs_df_all$compname_clean <- obs_df_all$compname
obs_df_all$compname <- NULL # Clean up compname
obs_df_all$taxonname <- NULL # Clean up compname
obs_df_links <- left_join(obs_df_links,obs_df_all[!duplicated(obs_df_all$pid),],by="pid") 
## Reattach original (uncleaned) component name
obs_df_links <- left_join(obs_df_links,comp.df[,c("cokey","compname")],by="cokey")

## Trim down to just desired columns
spatialpedoncomps <- obs_df_links[,2:116]

## Save matched pedon-component data for just the spatial matching
saveRDS(spatialpedoncomps, paste(resultfolder,"/NASIS_spatial_component_match_ssurgo20.rds",sep=""))

## Clean up explicitly linked pedons
linkcomps <- linked[!duplicated(linked$peiid),] # Grab one of the comp matches for each pedon - maybe a better way to do??
linkcomps$mtchtype <- "direct" # Noting that these were explicity component matches
linkcomps$cokey <- as.character(linkcomps$cokey)
linkcompsjn <- left_join(linkcomps[,c("cokey","peiid","mtchtype")],comp.df,by="cokey")
linkcompsjn$peiid <- as.character(linkcompsjn$peiid)
pts.df <- as.data.frame(shp.pts.all)
pts.df$peiid <- as.character(pts.df$peiid)
linkcompsjn <- left_join(linkcompsjn,pts.df[,c("peiid","taxonname_recent","x_std","y_std")],by="peiid")
## Create a clean compname for comparability with the spatially matched pedons
linkcompsjn$compname_clean <- toupper(linkcompsjn$compname)
## Clean up pt cleaned compnames for back compatibility
linkcompsjn$compname_clean = gsub(" VARIANT","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(" TAXADJUNCT","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(" TAXAJUNCT","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(" FAMILY","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(" LIKE","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub("-LIKE","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub("-SIMILAR","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(" SIMILAR","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(" TAX.","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(" TAX","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(" ERODED","", linkcompsjn$compname_clean)
linkcompsjn$compname_clean = gsub(", ERODED","", linkcompsjn$compname_clean)
# Eliminate link pedons without spatial coordinates
linkcompsjn <- linkcompsjn[!is.na(linkcompsjn$x_std),]

## Merge explicitly linked comps with spatially linked comps
allpedoncomps <- rbind(spatialpedoncomps,linkcompsjn)

## Save all pedon component matches (including explicit matches)
saveRDS(allpedoncomps, paste(resultfolder,"/NASIS_all_component_match_ssurgo20.rds",sep=""))

## Open SSURGO horizon table to join with components
horiz.df <- sf::st_read(dsn = ssurgo_fgdb, layer = "chorizon")
horiz.df[] <- lapply(horiz.df, function(x) if (is.factor(x)) as.character(x) else {x})

## Now select horizons and merge to create soil profile collection object
pedonhoriz <- inner_join(horiz.df,allpedoncomps[,c("cokey","peiid")],by="cokey")
pedoncomp.spc <- pedonhoriz
depths(pedoncomp.spc) <- peiid ~ hzdept_r + hzdepb_r # Create SPC
allpedoncomps$peiid <- as.character(allpedoncomps$peiid)
pedoncomp.spc@site <- left_join(pedoncomp.spc@site, allpedoncomps, by ="peiid") # Add component data
proj4string(pedoncomp.spc) <- polys.proj # Set spatial projection

## Save final soil profile collection
saveRDS(pedoncomp.spc, paste(resultfolder,"/NASIS_all_component_horizon_match_SPC_ssurgo20.rds",sep=""))
## Also Save as gpkg and shp
pts.proj <- proj4string(pedoncomp.spc)
pts <- pedoncomp.spc@site
coordinates(pts) <- ~ x_std + y_std
sp.pts <- as(pts,'SpatialPointsDataFrame')
projection(sp.pts) <- pts.proj
outname <- paste(resultfolder,"/NASIS_all_component_match_ssurgo20.gpkg",sep="")
outname.shp <- paste(resultfolder,"/",sep="")
writeOGR(sp.pts,dsn=outname,layer="NASIS_all_component_match_ssurgo20.gpkg",driver="GPKG")
writeOGR(sp.pts,dsn=outname.shp,layer="NASIS_all_component_match_ssurgo20_shp",driver="ESRI Shapefile",overwrite_layer = T)




