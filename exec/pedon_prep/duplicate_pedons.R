library(aqp)
library(sf)


# count lab pedons ----
load(file = "G:/Box/Box Sync/data/LDM-compact_20200709.RData")
lp <- ldm_bs; rm(ldm_bs)

lp_s <- site(lp)
lp_s <- subset(lp_s, !duplicated(pedlabsampnum))
nrow(lp_s)



# count field pedons ----
load(file = "C:/workspace2/nasis_pedons_20201112.RData")
fp <- spc


# check for duplicate horizons
fp_h <- horizons(fp)
fp_h_s <- aggregate(cbind(hzname, hzdept, hzdepb, texture, claytotest, fragvoltot, d_hue, m_hue) ~ peiid, data = fp_h, paste, na.action = na.pass)
fp_h_s$Z <- with(fp_h_s, paste(hzdept, hzdepb))
# fp_h_s$Z <- with(fp_h_s, paste(hzname, hzdept, hzdepb, texture, claytotest, fragvoltot, d_hue, m_hue))
sum(!duplicated(fp_h_s$idx))


# check of pedons within 1-meter of each other
fp_s <- {
  site(fp) ->.;
  transform(., 
            x_std = ifelse(is.na(x_std), x, x_std),
            y_std = ifelse(is.na(y_std), y, y_std)
  ) ->.;
}
fp_s_sf <- subset(fp_s, complete.cases(x_std, y_std))
fp_s_sf <- st_as_sf(
  fp_s_sf,
  coords = c("x_std", "y_std"),
  crs = 4326
)
fp_s_sf <- st_transform(fp_s_sf, crs = 5070)
fp_s_sf <- as.data.frame(cbind(fp_s_sf["peiid"], round(st_coordinates(fp_s_sf))))
fp_s_sf$geometry <- NULL

fp_s <- merge(fp_s, fp_s_sf, by = "peiid", all.x = TRUE)
fp_s <- merge(fp_s, fp_h_s,  by = "peiid", all.x = TRUE)

fp_s <- with(fp_s, fp_s[order(siteiid, pedlabsampnum, rcasiteid), ])

# remove duplicates
fp_s_sub <- subset(
  fp_s, 
  !duplicated(pedlabsampnum, incomparables = NA)
  & !duplicated(paste(pedlabsampnum, X, Y, Z))
  )

idx <- !is.na(fp_s_sub$pedlabsampnum) &
  nchar(fp_s_sub$pedlabsampnum) >= 6 & 
  fp_s_sub$pedlabsampnum != "Unknown"
nrow(fp_s_sub[idx,  ])
nrow(fp_s_sub[!idx, ])



# soil series ----
ss <- soilDB::get_soilseries_from_NASIS()
ss_n <- with (ss, data.frame(
  n_established = length(soilseriesstatus[soilseriesstatus == "established"]),
  n_inactive    = length(soilseriesstatus[soilseriesstatus == "inactive"]),
  n_tentative   = length(soilseriesstatus[soilseriesstatus == "tentative"]),
  n_benchmarks  = sum(benchmarksoilflag),
  n_statsgo     = sum(statsgoflag)
  ))
ss_n <- stack(ss_n)[c(2, 1)]
names(ss_n) <- c("soil series", "count")
knitr::kable(ss_n)
