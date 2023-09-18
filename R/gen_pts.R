#' Create random points for VI modeling
#'
#' Function that generates a user-defined number of randomly distributed points within a study area at an appropriate distance from the study area boundary. Study area extent defined by input digital terrain model (DTM) and digital surface model (DSM). These points are intended for use as training and validation points in the VisiMod visibility modeling workflow.
#' @details
#' * `dtm` and `dsm` SpatRasters can be defined using the terra library. They should have the same coordinate system, resolution, extent, and origin.
#' * The more points generated, the more robust the modeling procedure will be; however, more points will also increase total processing time for subsequent functions in the workflow. The default is set to 200, which should achieve a nice balance between model performance and processing time. We do not recommend generating more than 1000 points.
#' * The maximum VI distance will also affect processing time quite dramatically. As distance increases, processing time will increase exponentially. We do not recommend attempting this workflow at distances beyond 3000 m. 
#' @param dtm SpatRaster. Digital terrain model at finest resolution available
#' @param dsm SpatRaster. Digital surface model at finest resolution available
#' @param num_pts Numeric. The number of points generated
#' @param max_vis_dist Numeric. The maximum distance (in meters) to which VI will be calculated, to account for study area boundary
#' @return data.frame with columns 'x' and 'y'
#' @export
#' @examples
#' # read in dtm and dsm
#' dtm <- rast("dtm.tif")
#' dsm <- rast("dsm.tif")
#'
#' # generate 1,000 points with a buffer of 1,000 meters
#' generate_pts(dtm, dsm, 1000, 1000)

gen_pts <- function(dtm, dsm, num_pts = 200, max_vis_dist){
  
  # first make sure dtm and dsm are equivalent (alwayyys)
  cg <- terra::compareGeom(dtm, dsm, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)
  if (cg==TRUE){
    
    # get a single value version of the raster
    r1 <- dtm
    r1 <- terra::ifel(!is.na(r1), 1, NA)

    # turn raster to polygon
    p1 <- terra::as.polygons(r1)
    
    # buffer the max distance in case of wedge
    max_vis_dist <- max_vis_dist + 10
    
    # buffer polygon by distance
    b <- terra:::buffer(p1, -max_vis_dist)

    # sample points within the buffered polygon
    pts <- terra::spatSample(x = b, size = num_pts)

    # get data.frame of coordinates
    pts_df <- terra::geom(pts, df=TRUE)
    df <- pts_df |>
      dplyr::select(x, y)
    
    # return data.frame
    return(df)
    
  } else {
    
    stop("Raster geometry not comparable, check crs, extent, and resolution.")
    
  }
  
}
