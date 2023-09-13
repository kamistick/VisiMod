#' Create random points for VI modeling
#'
#' Function that generates a user-defined number of points at an appropriate distance from study area boundary
#' @param dtm SpatRaster. Digital terrain model at finest resolution available
#' @param dsm SpatRaster. Digital surface model at finest resolution available
#' @param num_pts Numeric. The number of points generated
#' @param max_vis_dist Numeric. The maximum distance (in meters) to which VI will be calculated, to account for study area boundary
#' @return dataframe with columns 'x' and 'y'
#' @export
#' @examples
#' # read in dtm and dsm
#' dtm <- rast("dtm.tif")
#' dsm <- rast("dsm.tif")
#'
#' # generate 1,000 points with a buffer of 1,000 meters
#' generate_pts(dtm, dsm, 1000, 1000)

generate_pts <- function(dtm, dsm, num_pts, max_vis_dist){
  # first make sure dtm and dsm are equivalent (alwayyys)
  cg <- terra::compareGeom(dtm, dsm, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)
  if (cg==TRUE){
    # get a single value version of the raster
    r1 <- dtm
    # r1[!is.na(r1)] <- 1
    r1 <- ifel(!is.na(r1), 1, NA)

    # turn raster to polygon
    p1 <- terra::as.polygons(r1)
    
    # buffer the max distance in case of wedge
    max_vis_dist <- max_vis_dist + 10
    
    # buffer polygon by distance
    b <- terra:::buffer(p1, -max_vis_dist)

    # sample points within the buffered polygon
    pts <- terra::spatSample(x = b, size = num_pts)

    pts_df <- terra::geom(pts, df=TRUE)

    df <- pts_df %>%
      dplyr::select(x, y)
    return(df)
  } else {(stop("Raster geometry not comparable, check crs, extent, and resolution."))
  }
}
