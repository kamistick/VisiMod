#' Calculate visibility index at points
#'
#' This Function uses both a digital terrain model (DTM) and digital surface model (DSM) to calculate the visibility index (VI) from a series of point locations. VI represents the proportion of area visible from a given point relative to the total area within a viewing distance of interest, ranging from 0 (no visibility) to 1 (complete visibility). 
#' 
#' This function enables the calculation of omnidirectional VI (in 360 degrees surrounding each point) or directional VI (within a viewing "wedge" defined by an azimuth and angular field of view). 
#' 
#' Directional visibility can be calculated in two ways: (1) in a singular, specific viewing direction/azimuth; or (2) with each point having a randomly assigned viewing direction/azimuth. The former is more useful for building VI predictive models in a singular direction; the latter can be used to build directionally independent models that can be applied for the prediction of visibility in any direction. Irrespective of the VI type, users define one or more distances or viewing radii within which VI is calculated. VI values are appended to input points are intended for use as training and validation in the VisiMod visibility modeling workflow.
#' 
#' @details
#' * This is the third suggested function in the VisiMod workflow, following (1) [VisiMod::prep_dems()] and (2) [VisiMod::gen_pts()]. 
#' * `dtm` and `dsm` SpatRasters can be defined using the terra library. They should have the same coordinate system, resolution, extent, and origin. They must be written on disk, they cannot only be held in memory. 
#' * `pts` can be defined using the [VisiMod::gen_pts()] function within the VisiMod library. But, they can also be created through many other means. For example, one could derive a data.frame of x-y coordinate pairs from a SpatVector using `terra::crds()` or from an sf object using `sf::st_coordinates()`. However, using [VisiMod::gen_pts()] is advantageous, as it ensures that points are (1) in the same coordinate system as `dtm` and `dsm`; and (2) are at least `vi_rad` from the edge of the study area.
#' * `vi_rad` should not exceed the `max_vi_rad` defined in the previous [VisiMod::gen_pts()] step in the VisiMod workflow. `vi_rad` will affect processing time quite dramatically. As distance increases, processing time will increase exponentially.
#' * Note that this function is parallelized and can leverage as many cores as your computer has available to speed up processing. As with all parallel processing in R, however, there is an overhead cost associated with setting up parallel operations. So, for small numbers of input points (`pts`) and/or short viewing radii (`view_rad`), using many cores may not speed up your processing significantly.
#' * Due to the parallel nature of this algorithm, there are rare situations in which errors occur in the calculation of VI from one (or more) of the input `pts`. We hypothesize that this may be the result of the same `dtm` and `dsm` being read into several cores at once, causing a conflict. In testing, we have found that by simply retrying a calculation one or more times, the errors can be resolved. The function argument `n_retry` enables the user to define how many times the calculation should be attempted on a point that returns an error before simply omitting it from further consideration. If the error persists, note that the number of inputs `pts` may not match the number of rows in the returned data.frame.
#' @param dtm SpatRaster. Digital terrain model at finest resolution available.
#' @param dsm SpatRaster. Digital surface model at finest resolution available.
#' @param pts data.frame. The point locations for which vi will be calculated. At a minimum, 'x' and 'y' columns are required. An 'azimuth' column is required if vi_type == 'directional_random'.
#' @param vi_type Character. Defines which type of VI calculation you would like to conduct. One of: "omnidir" (omnidirectional, 360 degree), "directional_single" (a single specified view direction for each point), "directional_random" (viewing directions will be randomly assigned to each point).
#' @param vi_rad Numeric. Viewing radius (in meters) from each point to which VI will be calculated, in meters. Can be defined as a single numeric value or a vector of numeric radii if you want to perform a multiscale analysis.
#' @param vi_fov Numeric. Defines the angular field of view, in degrees, of the directional wedge used for the VI calculation. Only used if vi_type == "directional_single" | vi_type == "directional_random". Values must be > 0 and < 360.
#' @param vi_azi Numeric. Defines the azimuth, or central viewing direction, in degrees, of the directional wedge used for VI calculation.  Only used if vi_type == "directional_single" | vi_type == "directional_random". Values must be 0-360.
#' @param cores Numeric. Defines the number of cores you would like to use for parallel processing. Defaults to half of the cores on your machine.
#' @param n_retry Numeric. Defines the number of times you would like to attempt to retry a VI calculation before omitting the point from further analysis. See Details.
#' @return data.frame with columns 'x', 'y', 'vi_x' for x in vi_rad, and 'azimuth' if `vi_type == 'directional_random'`
#'
#' @export
#' @examples
#' # get your dtm and dsm
#' dsm <- rast("dsm.tif")
#' dtm <- rast("dtm.tif")
#'
#' # check dtm and dsm
#' pd <- prep_dems(dtm, dsm, "C:/temp/dtm_filled.tif", "C:/temp/dsm_filled.tif")
#'
#' # get your points
#' my_points <- gen_pts(dtm, dsm, 100, 1000)
#' 
#' # calculate vi
#' my_points <- calc_vi(dtm, dsm, my_points, "directional_single", c(500, 1000), 90, 90, 4L, 5L)


calc_vi <- function(dtm, dsm, pts, vi_type, vi_rad, vi_fov=180, vi_azi = 0, 
                    cores = floor(parallel::detectCores()/2), n_retry = 5L){
  
  # create clean time printing function
  prttm <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # print message
  message(paste0(prttm(), " calc_vi() has begun"))
  
  # first compareGeom, it will return an error if the rasters do not match
  cg <- terra::compareGeom(dtm, dsm, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)
  if(cg==FALSE){
    stop("\nDTM and DSM geometry not comparable, check crs, extent, and resolution.\n")
  }
  
  # check that an x and y column exist
  if ("x" %in% colnames(pts) == FALSE | "y" %in% colnames(pts) == FALSE){
    stop("\nInput data.frame must have an 'x' and 'y' column.\n")
  }
  
  # get source locations of dsm and dtm files (needed for reading into each core)
  dsm.file <- terra::sources(dsm)
  dtm.file <- terra::sources(dtm)
  
  #---------------------------------omnidir------------------------------------#
  
  if (vi_type == "omnidir"){
    
    # begin parallelization
    clust <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(clust)
    `%dopar%` <- foreach::`%dopar%`
    
    # begin foreach loop
    result <- foreach::foreach(i = 1:nrow(pts),
                               .packages = c("terra", "sf"),
                               .combine = "rbind",
                               .inorder = T) %dopar% {
                                 
                                 # begin retry while loop
                                 success <- 0
                                 attempts <- 0
                                 while ((success == 0) && (attempts < n_retry)){
                                   try({
                                     attempts <- attempts + 1
                                     
                                     # read in dtm and dsm
                                     dtm <- terra::rast(dtm.file)
                                     dsm <- terra::rast(dsm.file)
                                     
                                     # isolate the point we are working with
                                     x <- pts[i,"x"]
                                     y <- pts[i,"y"]
                                     
                                     # make it a point
                                     pt <- sf::st_point(c(x,y))
                                     pt_sv <- terra::vect(pt)
                                     terra::crs(pt_sv) <- terra::crs(dtm)
                                     
                                     # crop your dsm and dtm to a buffer with radius of interest
                                     buf <- terra::buffer(pt_sv, max(vi_rad))
                                     dsm_buf_crop <- terra::crop(dsm, buf, "near", mask=TRUE)
                                     dtm_buf_crop <- terra::crop(dtm, buf, "near", mask=TRUE)
                                     
                                     # get the row/col for the the pixel in which the observer falls
                                     cell <- terra::cells(dtm_buf_crop, pt_sv)
                                     rowcol <- terra::rowColFromCell(dtm_buf_crop, cell[1,2])
                                     
                                     # get dtm value at that point/pixel
                                     dtm_at <- dtm_buf_crop[rowcol[1,1], rowcol[1,2]]
                                     
                                     # set dsm value at that spot
                                     dsm_buf_crop[rowcol[1,1], rowcol[1,2]] <- dtm_at[1,1]
                                     
                                     # run viewshed
                                     vs <- terra::viewshed(x = dsm_buf_crop, loc = c(x,y), 
                                                           observer = 1.7, target = 0, 
                                                           curvcoef = 0.85714, output = "yes/no")
                                     
                                     # create data.frame of results
                                     new_pts <- data.frame(x = x, y = y)
                                     
                                     # loop through viewing radii
                                     for (rad in vi_rad){
                                       
                                       # new buffer of appropriate distance
                                       buf_rad <- terra::buffer(pt_sv, rad)
                                       
                                       # get number of visible pixels within wedge
                                       zonal_sum <- terra::zonal(vs, buf_rad, "sum")
                                       
                                       # calculate vi
                                       vi <- zonal_sum[1,1] / ((pi * (rad^2)))
                                       
                                       # append vi column to results data.frame
                                       col_name <- paste0("vi_", rad)
                                       new_pts[,col_name] <- vi
                                       
                                     }
                                     
                                     # return data frame
                                     return(new_pts)
                                     success <- 1
                                     
                                   })
                                   
                                 }
                                 
                               }
    
    # end parallelization
    parallel::stopCluster(clust)
    
    #---------------------------directional_single-------------------------------#
    
  } else if (vi_type == "directional_single"){
    
    # begin parallelization
    clust <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(clust)
    `%dopar%` <- foreach::`%dopar%`
    
    # begin foreach loop
    result <- foreach::foreach(i = 1:nrow(pts),
                               .packages = c("terra", "sf", "VisiMod"),
                               .combine = "rbind",
                               .inorder = T) %dopar% {
                                 
                                 # begin retry while loop
                                 success <- 0
                                 attempts <- 0
                                 while ((success == 0) && (attempts < n_retry)){
                                   try({
                                     attempts <- attempts + 1
                                     
                                     # read in dtm and dsm
                                     dtm <- terra::rast(dtm.file)
                                     dsm <- terra::rast(dsm.file)
                                     
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
                                     poly <- VisiMod::wedge(x,y, max(vi_rad), azim, vi_fov)
                                     
                                     # crop your dsm and dtm to a buffer with radius of interest
                                     buf <- terra::buffer(pt_sv, max(vi_rad + 10))
                                     dsm_buf_crop <- terra::crop(dsm, buf, "near", mask=TRUE)
                                     dtm_buf_crop <- terra::crop(dtm, buf, "near", mask=TRUE)
                                     
                                     # get the row/col for the the pixel in which the observer falls
                                     cell <- terra::cells(dtm_buf_crop, pt_sv)
                                     rowcol <- terra::rowColFromCell(dtm_buf_crop, cell[1,2])
                                     
                                     # get dtm value at that point/pixel
                                     dtm_at <- dtm_buf_crop[rowcol[1,1], rowcol[1,2]]
                                     
                                     # set dsm value at that spot
                                     dsm_buf_crop[rowcol[1,1], rowcol[1,2]] <- dtm_at[1,1]
                                     
                                     # run viewshed
                                     vs <- terra::viewshed(x = dsm_buf_crop, loc = c(x,y), 
                                                           observer = 1.7, target = 0, 
                                                           curvcoef = 0.85714, output = "yes/no")
                                     
                                     # create data.frame of results
                                     new_pts <- data.frame(x = x, y = y)
                                     
                                     # loop through viewing radii
                                     for (rad in vi_rad){
                                       
                                       # new wedge of appropriate distance
                                       poly_rad <- VisiMod::wedge(x,y, rad, azim, vi_fov)
                                       
                                       # get number of visible pixels within wedge
                                       zonal_sum <- terra::zonal(vs, poly_rad, "sum")
                                       
                                       # calculate vi
                                       vi <- zonal_sum[1,1] / ((pi * (rad^2))*(vi_fov/360))
                                       
                                       # append vi column to results data.frame
                                       col_name <- paste0("vi_", rad)
                                       new_pts[,col_name] <- vi
                                       
                                     }
                                     
                                     # return data frame
                                     return(new_pts)
                                     success <- 1
                                     
                                   })
                                   
                                 }
                                 
                               }
    
    # end parallelization
    parallel::stopCluster(clust)
    
    #---------------------------directional_random-------------------------------#
    
  } else if (vi_type == "directional_random"){
    
    # get original pts
    pts_og <- pts
    
    # begin parallelization
    clust <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(clust)
    `%dopar%` <- foreach::`%dopar%`
    
    # begin foreach loop
    result <- foreach::foreach(i = 1:nrow(pts),
                               .packages = c("terra", "sf", "VisiMod"),
                               .combine = "rbind",
                               .inorder = T) %dopar% {
                                 
                                 # begin retry while loop
                                 success <- 0
                                 attempts <- 0
                                 while ((success == 0) && (attempts < n_retry)){
                                   try({
                                     attempts <- attempts + 1
                                     
                                     # get your random azimuth:
                                     # first check if it exists, if it does grab from col 
                                     # if not make it and add to col
                                     if ("azimuth" %in% colnames(pts_og)){
                                       azim <- pts[i,"azimuth"]
                                     } else {
                                       azim <- sample(359,1)
                                       pts[i,"azimuth"] <- azim
                                     }
                                     
                                     # read in dtm and dsm
                                     dtm <- terra::rast(dtm.file)
                                     dsm <- terra::rast(dsm.file)
                                     
                                     # isolate the point we are working with
                                     x <- pts[i,"x"]
                                     y <- pts[i,"y"]
                                     
                                     # make it a point
                                     pt <- sf::st_point(c(x,y))
                                     pt_sv <- terra::vect(pt)
                                     terra::crs(pt_sv) <- terra::crs(dtm)
                                     
                                     # draw your wedge (at your maximum view distance)
                                     poly <- VisiMod::wedge(x,y, max(vi_rad), azim, vi_fov)
                                     
                                     # crop your dsm and to a buffer with radius of interest
                                     buf <- terra::buffer(pt_sv, max(vi_rad) + 10)
                                     dsm_buf_crop <- terra::crop(dsm, buf, "near", mask=TRUE)
                                     dtm_buf_crop <- terra::crop(dtm, buf, "near", mask=TRUE)
                                     
                                     # get the row/col for the the pixel in which the observer falls
                                     cell <- terra::cells(dtm_buf_crop, pt_sv)
                                     rowcol <- terra::rowColFromCell(dtm_buf_crop, cell[1,2])
                                     
                                     # get dtm value at that point/pixel
                                     dtm_at <- dtm_buf_crop[rowcol[1,1], rowcol[1,2]]
                                     
                                     # set dsm value at that spot
                                     dsm_buf_crop[rowcol[1,1], rowcol[1,2]] <- dtm_at[1,1]
                                     
                                     # run viewshed
                                     vs <- terra::viewshed(x = dsm_buf_crop, loc = c(x,y), 
                                                           observer = 1.7, target = 0, 
                                                           curvcoef = 0.85714, output = "yes/no")
                                     
                                     # create data.frame of results
                                     new_pts <- data.frame(x = x, y = y, azimuth = azim)
                                     
                                     # loop through viewing radii
                                     for (rad in vi_rad){
                                       
                                       # new wedge of appropriate distance
                                       poly_rad <- VisiMod::wedge(x,y, rad, azim, vi_fov)
                                       
                                       # get number of visible pixels within wedge
                                       zonal_sum <- terra::zonal(vs, poly_rad, "sum")
                                       
                                       # calculate vi
                                       vi <- zonal_sum[1,1] / ((pi * (rad^2))*(vi_fov/360))
                                       
                                       # append vi column to results data.frame
                                       col_name <- paste0("vi_", rad)
                                       new_pts[,col_name] <- vi
                                       
                                     }
                                     
                                     # return data frame
                                     return(new_pts)
                                     success <- 1
                                     
                                   })
                                   
                                 }
                                 
                               }
    
    # end parallelization
    parallel::stopCluster(clust)
    
  } else {
    
    stop("\nvi_type must be one of: 'omnidir', 'directional_single', or 'directional_random'.\n")
    
  }
  
  # return result
  message(paste0(prttm(), " calc_vi() is complete"))
  return(result)
  
}