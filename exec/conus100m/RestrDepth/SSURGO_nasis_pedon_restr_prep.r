###################################
### Script to extract property data from gSSURGO, nasis pedons and SCD pedons
### for soil restriction depth and depth to lithic contact

## Packages
required.packages <- c("raster", "sp", "rgdal","doParallel", "plyr", "dplyr", "ncdf4","maptools", "rgeos","stats","spdep","sf","aqp","RSQLite")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
##Raster settings
rasterOptions(maxmemory = 2e+09, chunksize = 2e+08)

## Key folders
resultfolder <- "/mnt/covs/solus_preds/v2tst_gnat_pts/RestrDepth_gRPI_250k"
ptsfolder <- "/mnt/solus100/2020_Pedons"
ssurgo_fgdb <- "/media/nped/gSSURGO20/gSSURGO_CONUS.gdb"
ssurgofolder <- "/media/nped/gSSURGO20"
nasisextfldr <- "/mnt/solus100/NASIS_SSURGO_Extracts/NASIS20_SSURGO20_ext_final"


### Explore gSSURGO GDB and open necessary files
# Vignettes to read all GDB files
fc_list <- ogrListLayers(ssurgo_fgdb) # list layers
# polys <- readOGR(ssurgo_fgdb, "MUPOLYGON") # Takes 24 hrs to load with SSURGO for all USA
# saveRDS(polys,paste(ssurgofolder,"/gSSURGO20_conus_mupolys.rds",sep=""))
polys <- readRDS(paste(ssurgofolder,"/gSSURGO20_conus_mupolys.rds",sep=""))
polys@data[] <- lapply(polys@data, function(x) if (is.factor(x)) as.character(x) else {x}) # factors to char
polys.proj <- projection(polys)

## CONUS boundary to clip psuedopoints
polybound <- readOGR("/media/sped/GIS_data/US_Census_500k/cb_2020_us_state_500k", "conus_bound")
polybound <- spTransform(polybound, polys.proj)

## List of unique mukeys in area
# mukeys.df <- sf::st_read(dsn = ssurgo_fgdb, layer = "mapunit")
# mukeys.df[] <- lapply(mukeys.df, function(x) if (is.factor(x)) as.character(x) else {x})
# saveRDS(mukeys.df, paste(ssurgofolder,"/mukeysdf_conus.rds", sep=""))
mukeys.df <- readRDS(paste(ssurgofolder,"/mukeysdf_conus.rds", sep=""))
mukeys <- unique(mukeys.df$mukey)

## Also load necessary data tables and subset to study area
## Component table
# comp.df <- sf::st_read(dsn = ssurgo_fgdb, layer = "component")
# comp.df[] <- lapply(comp.df, function(x) if (is.factor(x)) as.character(x) else {x})
# comps <- comp.df[comp.df$mukey %in% mukeys,]
# comps$compname <- toupper(comps$compname)
# # Clean up compnames for matching
# comps$compname = gsub(" VARIANT","", comps$compname)
# comps$compname = gsub(" TAXADJUNCT","", comps$compname)
# comps$compname = gsub(" TAXAJUNCT","", comps$compname)
# comps$compname = gsub(" FAMILY","", comps$compname)
# comps$compname = gsub(" LIKE","", comps$compname)
# comps$compname = gsub("-LIKE","", comps$compname)
# comps$compname = gsub("-SIMILAR","", comps$compname)
# comps$compname = gsub(" SIMILAR","", comps$compname)
# comps$compname = gsub(" TAX.","", comps$compname)
# comps$compname = gsub(" TAX","", comps$compname)
# comps$compname = gsub(" ERODED","", comps$compname)
# comps$compname = gsub(", ERODED","", comps$compname)
# saveRDS(comps, paste(ssurgofolder,"/comps_conus.rds", sep=""))
comps <- readRDS(paste(ssurgofolder,"/comps_conus.rds", sep=""))
cokeys <- unique(comps$cokey)

## Pick out rock outcrop components (from Colby Brungard and Shawn Salley)
# subset out bedrock type components
CMU.5 <- comps[grepl(comps$compname, pattern = "rock outcrop", ignore.case = TRUE),]
CMU.7 <- comps[grepl(comps$compname, pattern = "outcrop", ignore.case = TRUE),]
CMU.9 <- comps[grepl(comps$compname, pattern = "bedrock", ignore.case = TRUE),]
CMU.10 <- do.call("rbind", list(CMU.5, CMU.7, CMU.9))
CMU.11 <- unique(CMU.10)
CMU.12 <- CMU.11[!(CMU.11$compkind %in% c("Series", "Taxon above family", "Family", "Taxadjunct", "Variant")), ]
CMU.13 <- CMU.12[CMU.12$comppct_r > 49,]
CMU.14 <- CMU.13[!is.na(CMU.13$compname),]
outcrop_mukeys <- unique(CMU.14$mukey)

## Outcrop polys
outcr_polys <- polys[polys$MUKEY %in% outcrop_mukeys,]

####### Now create a random sample of points in outcrop polygons
outcr_polypts <- spsample(outcr_polys, 10000, type = 'random')
outcr_polypts <- SpatialPointsDataFrame(outcr_polypts, data.frame(row.names=row.names(outcr_polypts), ID=1:length(outcr_polypts)))

## Check
plot(polybound)
plot(outcr_polypts, add =T)
outcr_polypts <- outcr_polypts[polybound,] # trim down
saveRDS(outcr_polypts, paste0(resultfolder,"/outcr_polypts.rds"))

############################### Now label existing pedons with restriction depths

####### NASIS
## Load nasis pedons pre-matched to components (including explicit matches)
nasispts <- readRDS(paste0(nasisextfldr,"/NASIS_all_component_match_ssurgo20.rds"))
nas_cokeys <- unique(nasispts$cokey)
# nas_comps <- comps[comps$cokey %in% mukeys,] #?

## SSURGO restriction table to join with nasis pts by comp key
restr.df <- sf::st_read(dsn = ssurgo_fgdb, layer = "corestrictions")
restr.df[] <- lapply(restr.df, function(x) if (is.factor(x)) as.character(x) else {x})
restrs <- restr.df[restr.df$cokey %in% nas_cokeys,]
## Order restrs to get first restriction for each category
restrs.all <- restrs[with(restrs,order(cokey,resdept_r)),]
# Updated method to join that leaves in comps with no restr layer
nasispts_rest <- cbind(restrs.all[match(nasispts$cokey, restrs.all$cokey),], nasispts)
## Adds max obs depth values of 201 where no restriction exists (as per sdvattribute nullratingreplacementvalue)
nasispts_rest$resdept_r <- ifelse(is.na(nasispts_rest$resdept_r), 201, nasispts_rest$resdept_r)
nasispts_rest$resdept_r <- ifelse(nasispts_rest$resdept_r > 201, 201, nasispts_rest$resdept_r) # REclass a couple >201 to 201 for consistent right censorship
saveRDS(nasispts_rest,paste0(resultfolder,"/nasis_restrs.rds"))

## Depth paralithic or bedrock lithics contact
restrs_anylithic <- subset(restrs.all, restrs.all$reskind == "Paralithic bedrock"|restrs.all$reskind=="Lithic bedrock")
restrs_anylithic <- restrs_anylithic[with(restrs_anylithic,order(cokey,resdept_r)),]
nasispts_anylithic <- cbind(restrs_anylithic[match(nasispts$cokey, restrs_anylithic$cokey),], nasispts)
nasispts_anylithic$anylithicdpt <- ifelse(is.na(nasispts_anylithic$resdept_r), 201, nasispts_anylithic$resdept_r)
nasispts_anylithic$anylithicdpt <- ifelse(nasispts_anylithic$anylithicdpt > 201, 201, nasispts_anylithic$anylithicdpt)
saveRDS(nasispts_anylithic,paste0(resultfolder,"/nasis_anylithic.rds"))



######## SCD
## load scd pts and horizon data to determin restriction depth
## RSQlite workflow form https://github.com/ncss-tech/gsp-sas/blob/master/lab_data.Rmd
con <- dbConnect(RSQLite::SQLite(), paste0(ptsfolder,"/KSSL-data.sqlite"))
(ldm_names <- dbListTables(con))
ldm <- lapply(c("NCSS_Layer","NCSS_Site_Location","NCSS_Pedon_Taxonomy"), function(x) dbReadTable(con , x))
names(ldm) <- c("NCSS_Layer","NCSS_Site_Location","NCSS_Pedon_Taxonomy")
dbDisconnect(con)

## Now wrangle tables: From code provided by Colby Brungard
h <- ldm$NCSS_Layer

### Depth to restricting layer
# All those with some form of root restricting layer. Once does have to be a bit careful with this because this does catch old horizon designations of ir which is 's' in todays taxonomy, and it catches horizons labeled with 'and' in the word.
r <- h[grep('[mrdx]', h$hzn_desgn, ignore.case=TRUE),]
# Remove a few odd horizon designations that dont make sense or aren't root restricting
table(r$hzn_desgn)
r1 <- r[!(r$hzn_desgn == 'Bt/E(Btpart' | r$hzn_desgn== 'Bt/E(Epart' | r$hzn_desgn== 'Bt/Cr1' | r$hzn_desgn== '2Cdkng' | r$hzn_desgn== 'Cr/Bt') | r$hzn_desgn == 'Bt/Cr2',]
# Only the first row of each pedon (i.e., the upper boundary of the first horizon indicating a root restricting layer)
r2 <- r1[!duplicated(r1$site_key),]
summary(r2$Depth)
hist(r2$Depth,breaks = c(0,50,100,150,200,300,400,500,1768),freq=T)
# Assign depth
r2$Depth = r2$hzn_top
# Now remove all pedons with root restricting layers from the whole dataset
h2 <- subset(h, !(h$site_key %in% r2$site_key))
# Get the max recorded depth of pedons without a recorded root restricting layer
h3 <- ddply(h2, "site_key", function(x) max(x$hzn_bot))
names(h3)[2] <- 'Depth'
# check that there is only one depth value for each pedon h3[duplicated(h3$pedon_key),] # no dups 2/28/22
# Remove the pedons for which I couldn't determine a depth class, then set the remaining depth values to 201 (which is essentially a flag indicating 201+ cm and should be useful for survival analysis) because these depths (even if > 201 cm) only record the lower depth to which the pedon was excavated and are not physically meaningful.
h4 <- h3[complete.cases(h3),]
## These are pedons which had no restriction flag
summary(h4$Depth)
hist(h4$Depth,breaks = c(0,50,100,150,200,300,400,500,2590),freq=T)
# However, many have max depths much less than 201, so not appropriate to call all very deep
h4 <-  h4[h4$Depth > 100,] # Include all that are at least deep due to concerns that some are partially sampled pedons
# Join both datasets
rrl <- rbind(r2[,c("site_key","Depth")], h4)
summary(rrl$Depth)
rrl$resdept_r <- ifelse(rrl$Depth > 201, 201, rrl$Depth) ##
summary(rrl$resdept_r)
saveRDS(rrl[,c("site_key","resdept_r")], paste0(resultfolder,"/scd_restrs_v2.rds"))

## Depth to lithic contact
### Depth to restricting layer
# All those with some form of lithic layer. Once does have to be a bit careful with this because this does catch old horizon designations of ir which is 's' in todays taxonomy, and it catches horizons labeled with 'and' in the word.
r <- h[grep('[r]', h$hzn_desgn, ignore.case=TRUE),]
# Remove a few odd horizon designations that dont make sense or aren't root restricting
table(r$hzn_desgn)
r1 <- r[!(r$hzn_desgn == 'Bt/E(Btpart' | r$hzn_desgn== 'Bt/E(Epart' | r$hzn_desgn== 'Bt/Cr1' | r$hzn_desgn== '2Cdkng' | r$hzn_desgn== 'Cr/Bt') | r$hzn_desgn == 'Bt/Cr2',]
# Only the first row of each pedon (i.e., the upper boundary of the first horizon indicating a root restricting layer)
r2 <- r1[!duplicated(r1$site_key),]
# Assign depth
r2$Depth = r2$hzn_top
# Now remove all pedons with lithic layers from the whole dataset
h2 <- subset(h, !(h$site_key %in% r2$site_key))
# Get the max recorded depth of pedons without a recorded lithic layer
h3 <- ddply(h2, "site_key", function(x) max(x$hzn_bot))
names(h3)[2] <- 'Depth'
# check that there is only one depth value for each pedon h3[duplicated(h3$pedon_key),] # no dups 2/28/22
# Remove the pedons for which I couldn't determine a depth class, then set the remaining depth values to 201
# (which is essentially a flag indicating 201+ cm and should be useful for survival analysis) because these depths (even if > 201 cm) only record the lower depth to which the pedon was excavated and are not physically meaningful.
h4 <- h3[complete.cases(h3),]
summary(h4$Depth)
hist(h4$Depth,breaks = c(0,50,100,150,200,300,400,500,2590),freq=T)
h4 <-  h4[h4$Depth > 100,]
# Join both datasets
rrl <- rbind(r2[,c("site_key","Depth")], h4)
summary(rrl$Depth)
rrl$Depth <- ifelse(rrl$Depth > 201, 201, rrl$Depth) ##
summary(rrl$Depth)
saveRDS(rrl[,c("site_key","Depth")], paste0(resultfolder,"/scd_anylithic_v2.rds"))
