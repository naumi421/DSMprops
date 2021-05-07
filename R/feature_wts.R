### Function to create case weights to balance sample bias in feature space.
# Assumes data preparation of rasters based on
# https://github.com/ncss-tech/dsm-properties/blob/main/balance_sample_weights.R
# lines 19-83, 5/3/2021
# TODO Make this function more generalizable to new raster inputs.
# This is very specific to Ramcharan et alk 2018 data https://doi.org/10.2136/sssaj2017.04.0122


#' Title Feature space weights generator
#' @description Function to weight a sample set to better match a reference set of raster environmental variables
#' for use in model case weights of in weighted sub-sampling schemes. The function is currently set up to run with specific
#' variables from Ramcharan et al., 2018 SSSAJ. https://doi.org/10.2136/sssaj2017.04.0122.
#' Rasters at https://nrcs.box.com/s/p3dcg9lyw6ocmmlhtobm04ticax4thj5 as of 5/3/2021.
#' @param ref_df data.frame of pre-extracted raster values sampled from the broader data sets to serve as a reference for
#' weighting the observed samples.
#' @param obs_df data.frame of observed field samples with pre-extracted values of the rasters included in the reference dataset
#' @param vars character vector of the variable names being used for the weights
#' @param probs list of probability breaks used for each variable in the weighting
#'
#' @return
#' @export
feat_wts <- function(ref_df, obs_df, vars, probs){

  vars_l <- list(
    # DEM = "DEMNED6",
    SLOPE = "SLPNED6",
    POS   = "POSNED6",
    # NLCD  = "NLCD116",
    # PM    = "PMTGSS7",
    EVI   = paste0("EX", 1:6, "MOD5"),
    PPT   = paste0("P",  formatC(1:12, width = 2, flag = "0"), "PRI5"),
    TEMP  = paste0("T",  formatC(1:12, width = 2, flag = "0"), "PRI5")
  )

  # transform pts and rasters ----
  #obs_df <- as.data.frame(pts)
  ppt_idx  <- names(obs_df) %in% vars_l$PPT[-1]
  temp_idx <- names(obs_df) %in% vars_l$TEMP
  evi_idx  <- names(obs_df) %in% vars_l$EVI
  obs_df$PPT  <- apply(obs_df[ppt_idx],  1, sum,  na.rm = TRUE)
  obs_df$TEMP <- apply(obs_df[temp_idx], 1, mean, na.rm = TRUE)
  obs_df$EVI  <- apply(obs_df[evi_idx],  1, mean, na.rm = TRUE)
  # obs_df$PM   <- as.factor(obs_df$PMTGSS7)

  #ref_df    <- as.data.frame(ref_df)
  ppt_idx   <- names(ref_df) %in% vars_l$PPT[-1]
  temp_idx  <- names(ref_df) %in% vars_l$TEMP
  evi_idx   <- names(ref_df) %in% vars_l$EVI
  ref_df$PPT  <- apply(ref_df[ppt_idx],  1, sum, na.rm = TRUE)
  ref_df$TEMP <- apply(ref_df[temp_idx], 1, mean, na.rm = TRUE)
  ref_df$EVI  <- apply(ref_df[evi_idx],  1, mean, na.rm = TRUE)
  # ref_df$PM <- as.factor(ref_df$PMTGSS7)
  # ref_df$NLCD <- as.factor(ref_df$NLCD116)

  # tidy variables
  ref_df <- ref_df[vars]
  obs_df <- obs_df[vars]

  # count inputs
  n_ref   <- ncol(ref_df)
  n_obs   <- ncol(obs_df)
  n_vars  <- length(vars)

  if (!is.list(probs)) {
    probs <- list(probs)[rep(1, n_ref)]
    n_probs <- n_ref
  } else n_probs <- length(probs)

  # add id
  ref_df$id <- 1:nrow(ref_df)
  obs_df$id <- 1:nrow(obs_df)

  # check
  test <- all.equal(n_ref, n_obs, n_vars, n_probs)
  if (!test) stop("ref_df and obs_df must contain all vars and be the same length")

  # tabulate ref
  brks <- lapply(vars, function(x) {
    if (is.numeric(ref_df[, x])) {
      quantile(ref_df[, x], probs = probs[[x]], na.rm = TRUE)
    } else levels(as.factor(ref_df[, x]))
  })
  names(brks) <- vars

  # Split up ref into quantile classes
  ref_brks <- ref_df
  ref_brks[1:n_ref] <- lapply(vars, function(x) {
    if (is.numeric(ref_df[, x])) {
      # brks <- quantile(ref_df[x], probs = probs[[x]], na.rm = TRUE)
      as.integer(cut(ref_df[, x], breaks = unique(brks[[x]]), include.lowest = TRUE))
    } else as.integer(ref_df[, x])
  })
  ref_brks$source <- "ref"
  ref_brks[1:n_ref] <- lapply(vars, function(x) {
    formatC(ref_brks[, x], width = 2, flag = "0")
  })
  ref_brks$interval <- apply(ref_brks[vars], 1, paste0, collapse = "-")

  # Compute reference probabilities for each interval
  ref_pct <- lapply("interval", function(x) {
    temp           <- prop.table(table(ref_brks[x])) * 100
    temp           <- as.data.frame.table(temp)
    temp$Var1      <- as.character(temp$Var1)
    names(temp)[1] <- "interval"
    temp$var <- x
    temp$source <- "ref"
    return(temp)
  })
  ref_pct <- do.call("rbind", ref_pct)

  # tabulate obs classes from brks
  obs_brks <- obs_df
  obs_brks[1:n_obs] <- lapply(vars, function(x) {
    if (is.numeric(obs_df[, x])) {
      as.integer(cut(obs_df[, x], breaks = unique(brks[[x]]), include.lowest = TRUE))
    } else as.integer(obs_df[, x])
  })
  obs_brks$source <- "obs"
  obs_brks[1:n_obs] <- lapply(vars, function(x) {
    formatC(obs_brks[, x], width = 2, flag = "0")
  })
  obs_brks$interval <- apply(obs_brks[vars], 1, paste0, collapse = "-")

  # Calculate observation probabilities for breaks
  obs_pct <- lapply("interval", function(x) {
    temp           <- prop.table(table(obs_brks[x])) * 100
    temp           <- as.data.frame.table(temp)
    temp$Var1      <- as.character(temp$Var1)
    names(temp)[1] <- "interval"
    temp$var <- x
    temp$source <- "obs"
    return(temp)
  })
  obs_pct <- do.call("rbind", obs_pct)

  # compute weights
  df <- merge(ref_pct, obs_pct, by = c("interval"), all.x = TRUE)

  cw <- function(ref, obs) {
    # pct <- 1 - round(obs / ref, 4) + 1
    pct <- ref / obs # or 1 / (obs / ref)
    return(pct)
  }

  df$wts <- cw(df$Freq.x, df$Freq.y)
  # df$wts2 <- df$wts + min(df$wts, na.rm = TRUE) * -1
  # df$wts2 <- df$wts2 / max(df$wts2, na.rm = TRUE)
  # df$wts <- ifelse(df$wts < 0.01, 0.01, df$wts)
  # df$wts2 <- df$wts / max(df$wts, na.rm = TRUE)

  # merge results and tidy
  obs_brks <- merge(obs_brks, df[c("interval", "wts")], by = "interval")
  obs_df   <- merge(obs_df, obs_brks[c("id", "wts")], by  = "id", all.x = TRUE)
  obs_df   <- obs_df[order(obs_df$id), ]
  ref_df$id <- NULL

  return(list(wts = obs_df$wts, brks = brks))

}
