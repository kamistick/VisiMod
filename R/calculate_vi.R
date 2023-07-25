#' Calculate visibility index at points
#'
#' Uses both a dtm and dsm to calculate visibility index (visible area/ total area) at points.
#' @param dtm SpatRaster digital terrain model at finest resolution available
#' @param dsm SpatRaster digital surface model at finest resolution available
#' @param pts A dataframe of the points (with an 'x' and 'y' column) for which vi will be calculated
#' @param vi_type which type of VI calculation would you like to conduct? One of: omnidir (omnidirectional, 360 degree), directional_single (a single specified view direction), directional_random (multiple random view directions will be sampled)
#' @param vi_dist distance to which vi will be calculated (can be a list of distances if you want to look at multiple (e.g. c(100,200,300)))
#' @param vi_fov optional argument, with a default value of 180, if vi_type "directional_single" or "directional_random" are selected. Specifies the field of view for the vi calculation
#' @param vi_azi optional argument, with a default value of 0, if vi_type "directional_single" specify the direction of vi in degrees values between 0 and 359 (e.g. North = 0, South = 180, East = 90, West = 270)
#' @return dataframe with columns 'x', 'y', 'vi_x' for x in vi_dist, and 'azimuth' if vi_type = 'directional_random'
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

calculate_vi <- function(dtm, dsm, pts, vi_type, vi_dist, vi_fov=180, vi_azi = 0){
  # first compareGeom, it will return an error if the rasters do not match
  cg <- terra::compareGeom(dtm, dsm, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)
  if(cg==FALSE){
    stop("DTM and DSM geometry not comparable, check crs, extent, and resolution.")
  }
  # then have dif options for vis_types
  if (vi_type == "omnidir"){
    # go through points, calculate vi
    for (i in 1:nrow(pts)){
      if(i %% 25 == 0 | i == 1){
        print(paste0("Starting point ", as.character(i), " at ", Sys.time(), " ..."))
      }
      # isolate the point we are working with
      x <- pts[i,1]
      y <- pts[i,2]

      # make it a point
      pt <- sf::st_point(c(x,y))
      pt_sv <- terra::vect(pt)
      terra::crs(pt_sv) <- terra::crs(dtm)

      # crop your dsm AND DTM to a BUFFER with radius of interest
      buf <- terra::buffer(pt_sv, max(vi_dist))
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

      for (dist in vi_dist){
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
      x <- pts[i,1]
      y <- pts[i,2]

      # make it a point
      pt <- sf::st_point(c(x,y))
      pt_sv <- terra::vect(pt)
      terra::crs(pt_sv) <- terra::crs(dtm)

      # draw your wedge (at your maximum view distance)
      poly <- wedge(x,y, max(vi_dist), azim, vi_fov)

      # crop your dsm AND DTM to a BUFFER with radius of interest
      buf <- terra::buffer(pt_sv, max(vi_dist+10))
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

      for (dist in vi_dist){
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
      poly <- wedge(x,y, max(vi_dist), azim, vi_fov)

      # crop your dsm AND DTM to a BUFFER with radius of interest
      buf <- terra::buffer(pt_sv, max(vi_dist)+10)
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

      for (dist in vi_dist){
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
