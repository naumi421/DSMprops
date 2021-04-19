#DSM Soil Properties Team
# weighting levels by geocodesource and date range work flow
# Dave White 03/31/21


setwd("D:/DSM_focus_team/PropertiesTeam/conusPedonData")

#Bring in site date
load("nasis_sites_20210325.RData")
names(s)

# remove any observations where utmeasting, utmnorthing, x, y is all null thus removing any sites without location values, if all 4 of these are null there is no valid location data
s.sub <- s[!with(s, is.na(utmeasting) & is.na(utmnorthing) & is.na(x) & is.na(y)),]


# investigate geocodesource
# checking levels of geocodesource
levels(as.factor(s.sub$geocoordsource))

# replace the NA values in geocodesource with not populated
s.sub$geocoordsource[is.na(s.sub$geocoordsource)] <- "not populated"

levels(as.factor(s.sub$geocoordsource))


# investigate obsdate and obsdatekind
levels(as.factor(s.sub$obsdate))
levels(as.factor(s$obsdatekind))


# subset based on geocodesource and obsdate
# first remove cols not needed
names(s.sub)
s.sub <- s.sub[,c(4,11,15:23)]
names(s.sub)

# check if na is false
nrow(s.sub) # number of total observations
length(subset(is.na(s.sub$geocoordsource), is.na(s.sub$geocoordsource)=="FALSE")) # number of observations that do not have NA for geocoordsource
length(subset(is.na(s.sub$obsdate), is.na(s.sub$obsdate)=="FALSE")) # number of observations that do not have NA for obsdate

# subset
s1 <- subset(s.sub, s.sub$geocoordsource=="auto-populated from survey grade gps")
s2 <- subset(s.sub, s.sub$geocoordsource=="auto-populated from gps")
s3 <- subset(s.sub, s.sub$geocoordsource=="imported from gps")
s4 <- subset(s.sub, s.sub$geocoordsource=="manually entered from post-validation")
s5 <- subset(s.sub, s.sub$geocoordsource=="manually entered from gps")
levels(as.factor(s.sub$geocoordsource))
# make sure obsdate is in date format for subsetting
s.sub$obsdate <- as.Date(s.sub$obsdate, format="%Y-%m-%d")
  
s6 <- subset(s.sub, geocoordsource %in% c('estimated from other source', 'not populated', 'unknown') & s.sub$obsdate >= "2010-01-01")
s7 <- subset(s.sub, geocoordsource %in% c('estimated from other source', 'not populated', 'unknown') & s.sub$obsdate < "2010-01-01" & s.sub$obsdate >= "2005-01-01")
s8 <- subset(s.sub, s.sub$geocoordsource %in% c('estimated from other source', 'not populated', 'unknown') & s.sub$obsdate < "2005-01-01")


# combine into one data.frame and assign wt_level
s.wt <- rbind(data.frame(s1[1], wt = 1),
              data.frame(s2[1], wt = 2),
              data.frame(s3[1], wt = 3),
              data.frame(s4[1], wt = 4),
              data.frame(s5[1], wt = 5),
              data.frame(s6[1], wt = 6),
              data.frame(s7[1], wt = 7),
              data.frame(s8[1], wt = 8))



saveRDS(s.wt, "geocode_weighting.RDS")

#####################

# condensed version if pre-processing steps are done, removing duplicates, removing null location data etc

#s.sub is the data.frame of site observations containing the peiid, geocodesource, and obsdate cols

# make sure obsdate is in date format for subsetting
s.sub$obsdate <- as.Date(s.sub$obsdate, format="%Y-%m-%d")

# subset
s1 <- subset(s.sub, s.sub$geocoordsource=="auto-populated from survey grade gps")
s2 <- subset(s.sub, s.sub$geocoordsource=="auto-populated from gps")
s3 <- subset(s.sub, s.sub$geocoordsource=="imported from gps")
s4 <- subset(s.sub, s.sub$geocoordsource=="manually entered from post-validation")
s5 <- subset(s.sub, s.sub$geocoordsource=="manually entered from gps")
s6 <- subset(s.sub, geocoordsource %in% c('estimated from other source', 'not populated', 'unknown') & s.sub$obsdate >= "2010-01-01")
s7 <- subset(s.sub, geocoordsource %in% c('estimated from other source', 'not populated', 'unknown') & s.sub$obsdate < "2010-01-01" & s.sub$obsdate >= "2005-01-01")
s8 <- subset(s.sub, s.sub$geocoordsource %in% c('estimated from other source', 'not populated', 'unknown') & s.sub$obsdate < "2005-01-01")


# combine into one data.frame and assign wt_level
s.wt <- rbind(data.frame(s1[1], wt = 1),
              data.frame(s2[1], wt = 2),
              data.frame(s3[1], wt = 3),
              data.frame(s4[1], wt = 4),
              data.frame(s5[1], wt = 5),
              data.frame(s6[1], wt = 6),
              data.frame(s7[1], wt = 7),
              data.frame(s8[1], wt = 8))

