

#' @title Extraction of covariate values to pedon locations to generate a regression matrix.
#' @description Query values of covariate rasters for the location of each pedon being used to train predictive models. It is implemented as a parallel list
#' apply to speed up implementation. The parallel implementation varies by operating system.
#'
#' @param sp SpatialPointsDataFrame object of all pedons
#' @param gridlist Vector list of full path names of all covariate grids. Can be created using list.files using pattern= and full.names=T.
#' @param os Operating system. Supports "windows" or "linux"
#' @param nthreads Number of logical processors to use in parallel functions
#'
#' @return a data.frame with same original fields as sp plus the extracted values of rasters as new columns
#' @export
#'
parPTextr <- function(sp, gridlist, os = "windows",nthreads = detectCores() - 1){

  if(nthreads == 0){nthreads <- 1}

  if(projection(sp)!=projection(raster(gridlist[1]))){
    message('reprojecting points to match grid projections')
    sp <- spTransform(sp, CRS(projection(raster(gridlist[1]))))
  }

  if(os == "linux"){
    cl <- makeCluster(nthreads, type="FORK")
    registerDoParallel(cl)
    ov.lst <- parLapply(cl,gridlist,function(i){try( raster::extract(raster(i), sp) )})
    stopCluster(cl)
    ov.lst <- as.data.frame(ov.lst)
    names(ov.lst) = tools::file_path_sans_ext(basename(gridlist))
    ov.lst$DID <- seq.int(nrow(ov.lst))
    sp$DID <- seq.int(nrow(sp))
    pts.ext <- merge(as.data.frame(sp),ov.lst, by="DID")
  }

  if(os == "windows"){
    sfInit(parallel=TRUE, cpus=nthreads)
    sfExport("sp", "gridlist")
    sfLibrary(raster)
    sfLibrary(rgdal)
    ov.lst <- sfLapply(gridlist, function(i){try( raster::extract(raster(i), sp) )})
    snowfall::sfStop()
    ov.lst <- as.data.frame(ov.lst)
    names(ov.lst) = tools::file_path_sans_ext(basename(gridlist))
    ov.lst$DID <- seq.int(nrow(ov.lst))
    sp$DID <- seq.int(nrow(sp))
    pts.ext <- merge(as.data.frame(sp),ov.lst, by="DID")
  }

  return(pts.ext)
}
