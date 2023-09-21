#' Create random points for VI modeling
#'
#' This function generates a user-defined number of randomly distributed points within a study area at an appropriate distance from the study area boundary. The study area extent is defined by the input digital terrain model (DTM) and digital surface model (DSM). These points are intended for use as training and validation points in the VisiMod visibility modeling workflow.
#' 
#' @details
#' * This is the second suggested function in the VisiMod workflow, following [VisiMod::prep_dems()]. 
#' * `dtm` and `dsm` SpatRasters can be defined using the terra library. They should have the same coordinate system, resolution, extent, and origin.
#' * The more points generated, the more robust the modeling procedure will be; however, more points will also increase processing time for subsequent functions in the workflow. The default is set to 200, which should balance model performance and processing time. We do not recommend generating more than 1000 points.
#' * As maximum VI radius (`max_vi_rad`) increases, processing time will increase exponentially. However, this also depends on the spatial resolution of the input `dtm`/`dsm`. For example, a 200m radius with a 1m resolution is functionally the same as a 400m distance with a 2m resolution, in terms of processing time. We do not recommend attempting this workflow at distances beyond 2000x the input resolution.. 
#' @param dtm SpatRaster. Digital terrain model at finest resolution available
#' @param dsm SpatRaster. Digital surface model at finest resolution available
#' @param num_pts Numeric. The number of points generated
#' @param max_vi_rad Numeric. The maximum radius (in meters) to which VI will be calculated. Defaults to 500x the resolution of the input `dtm`.
#' @return data.frame with columns 'x' and 'y'
#' @export
#' @examples
#' # read in dtm and dsm
#' dtm <- rast("dtm.tif")
#' dsm <- rast("dsm.tif")
#'
#' # generate 1,000 points with a buffer of 1,000 meters
#' gen_pts(dtm, dsm, 1000, 1000)

gen_pts <- function(dtm, dsm, num_pts = 200, max_vi_rad = terra::res(dtm)[1] * 500){
  
  # print message
  message(paste0(Sys.time(), " gen_pts() has begun"))
  
  # first compareGeom, it will return an error if the rasters do not match
  cg <- terra::compareGeom(dtm, dsm, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)
  if(cg==FALSE){
    stop("\nDTM and DSM geometry not comparable, check crs, extent, and resolution.\n")
  }
  
  # second make sure that max_vi_rad is small enough to leave sample area
  e <- terra::ext(dtm)
  wid <- abs(e[2] - e[1])
  hgt <- abs(e[4] - e[3])
  min.dim <- min(c(wid, hgt))
  if (min.dim <= 2 * max_vi_rad){
    stop("\nYour study area is too small for that max_vi_rad. Choose larger study area or smaller max_vi_rad before proceeding.\n")
  }
  
  # get a single value version of the raster
  message(paste0(Sys.time(), "   Creating mask of full data extent..."))
  r1 <- dtm
  r1 <- terra::ifel(!is.na(r1), 1, NA)
  
  # turn raster to polygon
  p1 <- terra::as.polygons(r1)
  
  # buffer the max distance in case of wedge
  max_vi_rad <- max_vi_rad + 10
  
  # buffer polygon by distance
  message(paste0(Sys.time(), "   Creating sample area via negative buffer of mask..."))
  b <- terra:::buffer(p1, -max_vi_rad)
  
  # sample points within the buffered polygon
  message(paste0(Sys.time(), "   Generating points within sample area..."))
  pts <- terra::spatSample(x = b, size = num_pts)
  
  # get data.frame of coordinates
  pts_df <- terra::geom(pts, df=TRUE)
  df <- pts_df |>
    dplyr::select(x, y)
  
  # return data.frame
  message(paste0(Sys.time(), " gen_pts() is complete"))
  return(df)
  
}
