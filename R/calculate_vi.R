#' Calculate visibility index at points
#'
#' Function that uses both a digital terrain model (DTM) and digital surface model (DSM) to calculate the visibility index (VI) from a series of point locations. VI represents the proportion of area visible from a given point relative to the total area within a viewing distance of interest, ranging from 0 (no visibility) to 1 (complete visibility). This function enables the calculation of omnidirectional VI (in 360 degrees surrounding each point) or directional VI (within a viewing "wedge" defined by an azimuth and angular field of view). Directional visibility can be calculated in two ways: (1) in a singular, specific viewing direction/azimuth; or (2) with each point having a randomly assigned viewing direction/azimuth. The former is more useful for building VI predictive models in a singular direction; the latter can be used to build directionally independent models that can be applied for the prediction of visibility in any direction. Irrespective of the VI type, users define one or more distances or viewing radii within which VI is calculated. VI values are appended to input points are intended for use as training and validation in the VisiMod visibility modeling workflow.
#' @details
#' * `dtm` and `dsm` SpatRasters can be defined using the terra library. They should have the same coordinate system, resolution, extent, and origin.
#' * `pts` can be defined using the `gen_pts()` function within the VisiMod library. But, they can also be created through many other means. For example, one could derive a data.frame of x-y coordinate pairs from a SpatVector using `terra::crds()` or from an sf object using `sf::st_coordinates()`. However, using `gen_pts()` is advantageous, as it ensures that points are (1) in the same coordinate system as `dtm` and `dsm`; and (2) are at least `vi_rad` from the edge of the study area.
#' * `vi_rad` should not exceed the `max_vis_dist` defined in the previous `gen_pts()` step in the VisiMod workflow. `vi_rad` will affect processing time quite dramatically. As distance increases, processing time will increase exponentially.
#' @param dtm SpatRaster. Digital terrain model at finest resolution available.
#' @param dsm SpatRaster. Digital surface model at finest resolution available.
#' @param pts dataframe.  The point locations for which vi will be calculated. At a minimum, 'x' and 'y' columns are required. An 'azimuth' column is required if vi_type == 'directional_random'.
#' @param vi_type Character. Defines which type of VI calculation you would like to conduct. One of: "omnidir" (omnidirectional, 360 degree), "directional_single" (a single specified view direction for each point), "directional_random" (viewing directions will be randomly assigned to each point).
#' @param vi_rad Numeric. Viewing radius (in meters) from each point to which VI will be calculated, in meters. Can be defined as a single numeric value or a vector of numeric radii if you want to perform a multiscale analysis.
#' @param vi_fov Numeric. Defines the angular field of view, in degrees, of the directional wedge used for the VI calculation. Only used if vi_type == "directional_single" | vi_type == "directional_random". Values must be > 0 and < 360.
#' @param vi_azi Numeric. Defines the azimuth, or central viewing direction, in degrees, of the directional wedge used for VI calculation.  Only used if vi_type == "directional_single" | vi_type == "directional_random". Values must be 0-360.
#' @return data.frame with columns 'x', 'y', 'vi_x' for x in vi_rad, and 'azimuth' if vi_type = 'directional_random'
#'
#' @export
#' @examples
#' # get your dtm and dsm
#' dsm <- rast("dsm.tif")
#' dtm <- rast("dtm.tif")
#'
#' # get your points
#' my_points <- generate_pts(dtm, dsm, 100, 1000)
#'
#' # calculate vi
#' df <- calculate_vi(dtm, dsm, my_points, "omnidir", c(500, 1000))

calculate_vi <- function(dtm, dsm, pts, vi_type, vi_rad, vi_fov=180, vi_azi = 0){
  # first compareGeom, it will return an error if the rasters do not match
  cg <- terra::compareGeom(dtm, dsm, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)
  if(cg==FALSE){
    stop("DTM and DSM geometry not comparable, check crs, extent, and resolution.")
  }
  
  # check that an x and y column exist
  if ("x" %in% colnames(pts) == FALSE | "y" %in% colnames(pts) == FALSE){
    stop("Input dataframe must have an 'x' and 'y' column.")
  }
  
  # then have dif options for vis_types
  if (vi_type == "omnidir"){
    # go through points, calculate vi
    for (i in 1:nrow(pts)){
      if(i %% 25 == 0 | i == 1){
        print(paste0("Starting point ", as.character(i), " at ", Sys.time(), " ..."))
      }
      # isolate the point we are working with
      x <- pts[i,"x"]
      y <- pts[i,"y"]

      # make it a point
      pt <- sf::st_point(c(x,y))
      pt_sv <- terra::vect(pt)
      terra::crs(pt_sv) <- terra::crs(dtm)

      # crop your dsm AND DTM to a BUFFER with radius of interest
      buf <- terra::buffer(pt_sv, max(vi_rad))
      dsm_buf_crop <- terra::crop(dsm, buf, "near", mask=TRUE)
      dtm_buf_crop <- terra::crop(dtm, buf, "near", mask=TRUE)

      # SET DSM PIXEL WHERE OBSERVER IS TO THE DTM VALUE AT THAT POINT
      cell <- terra::cells(dtm_buf_crop, pt_sv)
      rowcol <- terra::rowColFromCell(dtm_buf_crop, cell[1,2])

      # get DTM value at that point/pixel
      dtm_at <- dtm_buf_crop[rowcol[1,1], rowcol[1,2]]

      # SET DSM value at that spot
      dsm_buf_crop[rowcol[1,1], rowcol[1,2]] <- dtm_at[1,1]

      # run viewshed
      vs <- terra::viewshed(x = dsm_buf_crop, loc = c(x,y), observer = 1.7, target = 0, curvcoef = 0.85714, output = "yes/no")

      for (dist in vi_rad){
        # new buffer of appropriate distance
        buf_dist <- terra::buffer(pt_sv, dist)

        # get number of visible pixels within WEDGE
        zonal_sum <- terra::zonal(vs, buf_dist, "sum")

        # calculate vi
        vi <- zonal_sum[1,1] / ((pi * (dist^2)))
        col_name <- paste0("vi_", as.character(dist))
        pts[i, col_name] <- vi
      }
    }
  } else if (vi_type == "directional_single"){
    # calculate vi per point
    for (i in 1:nrow(pts)){
      if(i %% 25 == 0 | i == 1){
        print(paste0("Starting point ", as.character(i), " at ", Sys.time(), " ..."))
      }

      # set your single azimuth
      azim <- vi_azi

      # isolate the point we are working with
      x <- pts[i,"x"]
      y <- pts[i,"y"]

      # make it a point
      pt <- sf::st_point(c(x,y))
      pt_sv <- terra::vect(pt)
      terra::crs(pt_sv) <- terra::crs(dtm)

      # draw your wedge (at your maximum view distance)
      poly <- wedge(x,y, max(vi_rad), azim, vi_fov)

      # crop your dsm AND DTM to a BUFFER with radius of interest
      buf <- terra::buffer(pt_sv, max(vi_rad+10))
      dsm_buf_crop <- terra::crop(dsm, buf, "near", mask=TRUE)
      dtm_buf_crop <- terra::crop(dtm, buf, "near", mask=TRUE)

      # SET DSM PIXEL WHERE OBSERVER IS TO THE DTM VALUE AT THAT POINT
      cell <- terra::cells(dtm_buf_crop, pt_sv)
      rowcol <- terra::rowColFromCell(dtm_buf_crop, cell[1,2])

      # get DTM value at that point/pixel
      dtm_at <- dtm_buf_crop[rowcol[1,1], rowcol[1,2]]

      # SET DSM value at that spot
      dsm_buf_crop[rowcol[1,1], rowcol[1,2]] <- dtm_at[1,1]

      # run viewshed
      vs <- terra::viewshed(x = dsm_buf_crop, loc = c(x,y), observer = 1.7, target = 0, curvcoef = 0.85714, output = "yes/no")

      for (dist in vi_rad){
        # new wedge of appropriate distance
        poly_dist <- wedge(x,y, dist, azim, vi_fov)

        # get number of visible pixels within WEDGE
        zonal_sum <- terra::zonal(vs, poly_dist, "sum")

        # calculate vi
        vi <- zonal_sum[1,1] / ((pi * (dist^2))*(vi_fov/360))
        col_name <- paste0("vi_", as.character(dist))
        pts[i, col_name] <- vi
      }
    }

  } else if (vi_type == "directional_random") {
    # calculate vi per point
    pts_og <- pts
    for (i in 1:nrow(pts)){
      if(i %% 25 == 0 | i == 1){
        print(paste0("Starting point ", as.character(i), " at ", Sys.time(), " ..."))
      }

      # get your random azimuth:
      # first check if it exists, if it does grab from col if not make it and add to col
      if ("azimuth" %in% colnames(pts_og)){
        azim <- pts[i,"azimuth"]
      } else {
        azim <- sample(359,1)
        pts[i,"azimuth"] <- azim
      }

      # isolate the point we are working with
      x <- pts[i,"x"]
      y <- pts[i,"y"]

      # make it a point
      pt <- sf::st_point(c(x,y))
      pt_sv <- terra::vect(pt)
      terra::crs(pt_sv) <- terra::crs(dtm)

      # draw your wedge (at your maximum view distance)
      poly <- wedge(x,y, max(vi_rad), azim, vi_fov)

      # crop your dsm AND DTM to a BUFFER with radius of interest
      buf <- terra::buffer(pt_sv, max(vi_rad)+10)
      dsm_buf_crop <- terra::crop(dsm, buf, "near", mask=TRUE)
      dtm_buf_crop <- terra::crop(dtm, buf, "near", mask=TRUE)

      # SET DSM PIXEL WHERE OBSERVER IS TO THE DTM VALUE AT THAT POINT
      cell <- terra::cells(dtm_buf_crop, pt_sv)
      rowcol <- terra::rowColFromCell(dtm_buf_crop, cell[1,2])

      # get DTM value at that point/pixel
      dtm_at <- dtm_buf_crop[rowcol[1,1], rowcol[1,2]]

      # SET DSM value at that spot
      dsm_buf_crop[rowcol[1,1], rowcol[1,2]] <- dtm_at[1,1]

      # run viewshed
      vs <- terra::viewshed(x = dsm_buf_crop, loc = c(x,y), observer = 1.7, target = 0, curvcoef = 0.85714, output = "yes/no")

      for (dist in vi_rad){
        poly_dist <- wedge(x,y, dist, azim, vi_fov)
        # get number of visible pixels within WEDGE
        zonal_sum <- terra::zonal(vs, poly_dist, "sum")

        # calculate vi
        vi <- zonal_sum[1,1] / ((pi * (dist^2))*(vi_fov/360))
        col_name <- paste0("vi_", as.character(dist))
        pts[i, col_name] <- vi
      }
    }
  } else {
    stop("vi_type must be one of: 'omnidir', 'directional_single', or 'directional_random'.")
  }
  return(pts)
}
