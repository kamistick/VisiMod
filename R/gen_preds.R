#' Generate predictors for modeling visibility index
#'
#' The fourth suggested function in the VisiMod workflow, following (1) [VisiMod::prep_dems()], (2) [VisiMod::gen_pts()], and (3) [VisiMod::calc_vi()]. This function generates and extracts a suite of visibility predictor pixel values at each input point location. The predictor variables were designed to capture the dominant landscape controls of visibility in a wildland (undeveloped) environment that can be readily quantified using a combination of a DTM and DSM. They include canopy cover, canopy height, elevation, slope, derivative of slope, terrain curvature (plan and profile), cosine and sine of aspect, and slope x cosine and sine of aspect. Predictor values are extracted both locally (i.e., representing the terrain/vegetation characteristic of the pixel within which each point lies) and focally (i.e., within a focal neighborhood of pixels surrounding the point). In the case of focally derived predictors, the focal mean and standard deviation are calculated for each landscape variable within a neighborhood with radius `x` pixels for `x %in% c(2, 4, 6, 8, 16, 32)`. If the `vi_type == 'omnidir'`, the focal neighborhood will be a circle surrounding each point. If the `vi_type == 'directional_single' | vi_type == 'directional_random'`, the focal neighborhood will be a 'wedge' defined by the central viewing azimuth (`vi_azi`) and the angular field of view (`vi_fov`).
#' @details
#' * `dtm` and `dsm` SpatRasters can be defined using the terra library. They should have the same coordinate system, resolution, extent, and origin.
#' * `pts` should come from the successful execution of the previous function in the VisiMod workflow: [VisiMod::calc_vi()].
#' * `vi_type`, `vi_fov`, and `vi_azi` should match those used in your previous call to the [VisiMod::calc_vi()] function. In the VisiMod workflow, this will ensure that your visibility modeling response variable(s) (VI within a spatial scope of interest) will correspond with the spatial scope of your predictor variables.
#' * Note that if `vi_type == 'directional_random'`, this function will not be able to save a multiband raster with spatially-exhaustive ('wall-to-wall') representations of the predictor variables, as the directional predictor values are calculated according to the randomly-assigned azimuths of each input point. This is useful for building a robust, directionally-independent model that can generate directional VI predictions, but to apply that model to generate a VI map in a specific direction will require rerunning this function with `vi_type == 'directional_single'`. For most applications, it will be more useful and efficient to set `vi_type` to `'directional_single'`.
#' * Predictors are generated at a user-defined coarser resolution than the input DTM/DSM resolution to facilitate broad scale, efficient VI modeling and mapping. Defining the aggregation factor (`agg_fact`) to be too low will result in very large predictor files, consume larger proportions of RAM, and slow down the VisiMod workflow. Conversely, setting `agg_fact` too high, although computationally more efficient, may result in predictors that are too spatially coarse to capture sufficient landscape structural detail capable of predicting visibility. In testing, we have found predictors with spatial resolutions of approximately 10m to produce a desirable balance between these extremes.
#' * For the focal predictors, by default, the following pixel radii are used to define your focal anlaysis neighborhood: 2, 4, 6, 8, 16, and 32 pixels. It is possible that a circular or wedge neighborhood around a given point may extend beyond the bounds of your input DTM/DSM, particularly at the larger radii (i.e., 16- or 32-pixels). Whether this happens or not will be dependent upon two factors: (1) the `max_vi_rad` supplied to the function [VisiMod::gen_pts()]; and (2) the `agg_fact` specified in this function. For example, if you defined your `max_vi_rad` to be 200m, that means some of your sample points may be 200m from the edge of your DTM/DSM. If your DTM/DSM were 1m resolution, and you specified an `agg_fact` of 10, your predictor layers would have a 10m resolution. A 32-pixel radius circle or wedge would have a radius of 320m, meaning the neighborhood would reach beyond the extent of your study area. In situations like this, the focal predictors would all return NAs to ensure that incomplete neighborhoods aren't affecting the quality of the predictor data extracted at each point. So, you must think ahead about the spatial relationships between study area size, `max_vi_rad`, input and output resolutions to ensure that your data are complete.
#' 
#' @param dtm SpatRaster. Digital terrain model at finest resolution available.
#' @param dsm SpatRaster. Digital surface model at finest resolution available.
#' @param pts data.frame. The point locations for which vi will be calculated. At a minimum, 'x' and 'y' columns are required. An 'azimuth' column is required if `vi_type == 'directional_random'`. One or more 'vi' columns should also be present if this is run after [VisiMod::calc_vi()].
#' @param vi_type Character. Defines which type of VI calculation you would like to conduct. One of: 'omnidir' (omnidirectional, 360 degree), 'directional_single' (a single specified view direction for each point), 'directional_random' (viewing directions will be randomly assigned to each point).
#' @param vi_fov Numeric. Defines the angular field of view, in degrees, of the directional wedge used for the VI calculation. Only used if `vi_type == 'directional_single' | vi_type == 'directional_random'`. Values must be > 0 and < 360.
#' @param vi_azi Numeric. Defines the azimuth, or central viewing direction, in degrees, of the directional wedge used for VI calculation. Only used if `vi_type == 'directional_single'`. Values must be 0-360.
#' @param agg_fact Numeric. Defines the aggregation factor that determines the output resolution of your predictors (and your eventual VI map if you complete the VisiMod workflow). A multiplier applied to the input resolution of your DTM/DSM. Must be an integer.
#' @param save  Boolean. Defines whether or not you would like to save a multiband raster of predictors. Default is TRUE. Only used if `vi_type == 'omnidir' | vi_type == 'single_directional'`.
#' @param save_dir Character. If `save == TRUE`, this specifies the directory where the multiband raster of predictors should be saved. Only used if `vi_type == 'omnidir' | vi_type == 'directional_single'`. Default is your working directory.
#' @return If `vi_type == 'omnidir' | vi_type == 'directional_single'`, a list with two items: (1) `pred_pts` is a data.frame with predictor data associated with each input point; and (2) `pred_rast` is a SpatRaster multiband raster with all of the predictor data. If `save == TRUE` a TIF will also be written to file. If `vi_type == 'directional_random'`, then only a data.frame with predictor point values will be returned.
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
#'
#' # generate predictors
#' preds <- gen_preds(dtm, dsm, my_points, "directional_single", 90, 90, T, "C:/temp")

gen_preds <- function(dtm, dsm, pts, vi_type, vi_fov=180, vi_azi=0, agg_fact = 10L,
                      save = TRUE, save_dir = getwd()){
  
  # print message
  message(paste0(Sys.time(), " gen_preds() has begun"))
  
  # suppress progress bars
  terra::terraOptions(progress=0)
  
  # first compareGeom, it will return an error if the rasters do not match
  cg <- terra::compareGeom(dtm, dsm, crs=TRUE, ext=TRUE, rowcol=TRUE, res=TRUE)
  if(cg==FALSE){
    stop("\nDTM and DSM geometry not comparable, check crs, extent, and resolution.\n")
  }
  
  # make sure vi_type is one of the three options
  if (!vi_type %in% c("omnidir", "directional_single", "directional_random")){
    stop("\nvi_type must be one of: 'omnidir', 'directional_single', or 'directional_random'.\n")
  }
  
  # print message
  message(paste0(Sys.time(), "   Starting predictor generation..."))
  message(paste0(Sys.time(), "   Running local predictors..."))
  
  # generate canopy height model, fill negatives with zeros
  chm <- dsm-dtm
  chm <- terra::ifel(chm < 0, 0, chm)
  
  # aggregate dtm and chm by agg_fact
  dtm_agg <- terra::aggregate(dtm, agg_fact, mean)
  chm_agg <- terra::aggregate(chm, agg_fact, mean)
  
  ##############################################################################
  # omnidir and directional_single
  ##############################################################################
  
  #---------------local terrain predictors
  
  # elevation
  elevation <- dtm_agg
  names(elevation) <- "elevation"
  
  # slope
  slope <- terra::terrain(dtm_agg, "slope")
  names(slope) <- "slope"
  
  # slope derivative
  slope_derivative <- terra::terrain(slope, "slope")
  names(slope_derivative) <- "slope_derivative"
  
  # curvature
  curvature <- spatialEco::curvature(dtm_agg, "total")
  names(curvature) <- "curvature"
  
  # curvature_plan
  curvature_plan <- spatialEco::curvature(dtm_agg, "planform")
  names(curvature_plan) <- "curvature_plan"
  
  # curvature_prof
  curvature_prof <- spatialEco::curvature(dtm_agg, "profile")
  names(curvature_prof) <- "curvature_prof"
  
  # aspect
  aspect <- terra::terrain(dtm_agg, "aspect")
  names(aspect) <- "aspect"
  
  # aspect transformations -- depend upon vi_type
  if (vi_type == "omnidir") {
    
    # aspect_sin
    aspect_sin <- sin(aspect * pi/180)
    names(aspect_sin) <- "aspect_sin"
    
    # aspect_cos
    aspect_cos <- cos(aspect * pi/180)
    names(aspect_cos) <- "aspect_cos"
    
    # slope_aspect_sin
    slope_aspect_sin <- slope * aspect_sin
    names(slope_aspect_sin) <- "slope_aspect_sin"
    
    # slope_aspect_cos
    slope_aspect_cos <- slope * aspect_cos
    names(slope_aspect_cos) <- "slope_aspect_cos"
    
  } else if (vi_type == "directional_single") {
    
    # normalize aspect relative to viewing direction
    aspect <- aspect - vi_azi
    
    # aspect_sin
    aspect_sin <- sin(aspect * pi/180)
    names(aspect_sin) <- "aspect_sin"
    
    # aspect_cos
    aspect_cos <- cos(aspect * pi/180)
    names(aspect_cos) <- "aspect_cos"
    
    # slope_aspect_sin
    slope_aspect_sin <- slope * aspect_sin
    names(slope_aspect_sin) <- "slope_aspect_sin"
    
    # slope_aspect_cos
    slope_aspect_cos <- slope * aspect_cos
    names(slope_aspect_cos) <- "slope_aspect_cos"
    
  }
  
  #---------------local vegetation predictors
  
  # canopy height (ch)
  ch <- chm_agg
  names(ch) <- "ch"
  
  # canopy cover (cc)
  tnt <- chm # tnt = tree/non-tree, where tree = chm > 2m
  tnt <- terra::ifel(tnt > 2, 1, 0)
  cc <- terra::aggregate(tnt, agg_fact, mean)
  names(cc) <- "cc"
  
  # local stack
  stack.local <- c(elevation, slope, slope_derivative, 
                   curvature, curvature_plan, curvature_prof,
                   aspect_sin, aspect_cos, slope_aspect_sin, slope_aspect_cos,
                   ch, cc)
  
  #---------------focal predictors, omnidir and directional_single
  
  if (vi_type == "omnidir" | vi_type == "directional_single") {
    
    # loop through focal radii
    for (pred_rad in c(2,4,6,8,16,32)) {
      
      # print message
      message(paste0(Sys.time(), "   Running focal predictors at ",  pred_rad,  "px radius..."))
      
      # generate focal neighborhoods (matrices) for omnidir
      if (vi_type == "omnidir") {
        
        # add small amount to radius
        buff_dist <- pred_rad + 0.5
        
        # create two buffers, outer and inner
        outer_buf <- terra::buffer(terra::vect(cbind(0, 0),crs=terra::crs(dtm)),buff_dist)
        inner_buf <- terra::buffer(terra::vect(cbind(0,0),crs=terra::crs(dtm)),buff_dist/2)
        
        # create annulus (donut difference between outer and inner buffer)
        annulus <- terra::erase(outer_buf, inner_buf)
        
        # create template raster for rasterization
        r <- terra::rast(xmin = -buff_dist, xmax = buff_dist, 
                         ymax = buff_dist, ymin = -buff_dist, 
                         ncols = buff_dist * 2, nrows = buff_dist * 2, vals = 1)
        
        # create raster and matrix of outer buffer (for full circle focal area)
        full_ras <- terra::rasterize(outer_buf, r)
        full_mat <- as.matrix(full_ras, wide = TRUE)
        
        # create raster and matrix of annulus (for tpi calcs)
        annul_ras <- terra::rasterize(annulus, r)
        annul_mat <- as.matrix(annul_ras, wide = TRUE)
        
        # generate focal neighborhoods (matrices) for directional_single
      } else {
        
        # add small amount to radius
        buff_dist <- pred_rad + 0.5
        
        # create two wedges, outer and inner
        outer_wdg <- VisiMod::wedge(0, 0, buff_dist, vi_azi, vi_fov)
        inner_buf <- terra::buffer(terra::vect(cbind(0,0),crs=terra::crs(dtm)),buff_dist/2)
        
        # create annulus wedge (donut difference between outer and inner wedges)
        annulus <- terra::erase(outer_wdg, inner_buf)
        
        # create template raster for rasterization
        r <- terra::rast(xmin = -buff_dist, xmax = buff_dist, 
                         ymax = buff_dist, ymin = -buff_dist, 
                         ncols = buff_dist * 2, nrows = buff_dist * 2, vals = 1)
        
        # create raster and matrix of outer wedge (for full wedge focal area)
        full_ras <- terra::rasterize(outer_wdg, r, touches = T)
        full_mat <- as.matrix(full_ras, wide = TRUE)
        
        # create raster and matrix of annulus wedge (for tpi calcs)
        annul_ras <- terra::rasterize(annulus, r, touches = T)
        annul_mat <- as.matrix(annul_ras, wide = TRUE)
        
      }
      
      #---------------focal terrain predictors, mean
      
      # elevation_mean_x
      elevation_mean <- terra::focal(elevation, full_mat, mean)
      names(elevation_mean) <- paste0("elevation_mean_", pred_rad)
      
      # slope_mean_x
      slope_mean <- terra::focal(slope, full_mat, mean)
      names(slope_mean) <- paste0("slope_mean_", pred_rad)
      
      # slope_derivative_mean_x
      slope_derivative_mean <- terra::focal(slope_derivative, full_mat, mean)
      names(slope_derivative_mean) <- paste0("slope_derivative_mean_", pred_rad)
      
      # curvature_mean_x
      curvature_mean <- terra::focal(curvature, full_mat, mean)
      names(curvature_mean) <- paste0("curvature_mean_", pred_rad)
      
      # curvature_plan_mean_x
      curvature_plan_mean <- terra::focal(curvature_plan, full_mat, mean)
      names(curvature_plan_mean) <- paste0("curvature_plan_mean_", pred_rad)
      
      # curvature_prof_mean_x
      curvature_prof_mean <- terra::focal(curvature_prof, full_mat, mean)
      names(curvature_prof_mean) <- paste0("curvature_prof_mean_", pred_rad)
      
      # aspect_sin_mean_x
      aspect_sin_mean <- terra::focal(aspect_sin, full_mat, mean)
      names(aspect_sin_mean) <- paste0("aspect_sin_mean_", pred_rad)
      
      # aspect_cos_mean_x
      aspect_cos_mean <- terra::focal(aspect_cos, full_mat, mean)
      names(aspect_cos_mean) <- paste0("aspect_cos_mean_", pred_rad)
      
      # slope_aspect_sin_mean_x
      slope_aspect_sin_mean <- terra::focal(slope_aspect_sin, full_mat, mean)
      names(slope_aspect_sin_mean) <- paste0("slope_aspect_sin_mean_", pred_rad)
      
      # slope_aspect_cos_mean_x
      slope_aspect_cos_mean <- terra::focal(slope_aspect_cos, full_mat, mean)
      names(slope_aspect_cos_mean) <- paste0("slope_aspect_cos_mean_", pred_rad)
      
      #---------------focal vegetation predictors, mean
      
      # ch_mean_x
      ch_mean <- terra::focal(ch, full_mat, mean)
      names(ch_mean) <- paste0("ch_mean_", pred_rad)
      
      # cc_mean_x
      cc_mean <- terra::focal(cc, full_mat, mean)
      names(cc_mean) <- paste0("cc_mean_", pred_rad)
      
      #---------------focal terrain predictors, sd
      
      # elevation_sd_x
      elevation_sd <- terra::focal(elevation, full_mat, sd)
      names(elevation_sd) <- paste0("elevation_sd_", pred_rad)
      
      # slope_sd_x
      slope_sd <- terra::focal(slope, full_mat, sd)
      names(slope_sd) <- paste0("slope_sd_", pred_rad)
      
      # slope_derivative_sd_x
      slope_derivative_sd <- terra::focal(slope_derivative, full_mat, sd)
      names(slope_derivative_sd) <- paste0("slope_derivative_sd_", pred_rad)
      
      # curvature_sd_x
      curvature_sd <- terra::focal(curvature, full_mat, sd)
      names(curvature_sd) <- paste0("curvature_sd_", pred_rad)
      
      # curvature_plan_sd_x
      curvature_plan_sd <- terra::focal(curvature_plan, full_mat, sd)
      names(curvature_plan_sd) <- paste0("curvature_plan_sd_", pred_rad)
      
      # curvature_prof_sd_x
      curvature_prof_sd <- terra::focal(curvature_prof, full_mat, sd)
      names(curvature_prof_sd) <- paste0("curvature_prof_sd_", pred_rad)
      
      # aspect_sin_sd_x
      aspect_sin_sd <- terra::focal(aspect_sin, full_mat, sd)
      names(aspect_sin_sd) <- paste0("aspect_sin_sd_", pred_rad)
      
      # aspect_cos_sd_x
      aspect_cos_sd <- terra::focal(aspect_cos, full_mat, sd)
      names(aspect_cos_sd) <- paste0("aspect_cos_sd_", pred_rad)
      
      # slope_aspect_sin_sd_x
      slope_aspect_sin_sd <- terra::focal(slope_aspect_sin, full_mat, sd)
      names(slope_aspect_sin_sd) <- paste0("slope_aspect_sin_sd_", pred_rad)
      
      # slope_aspect_cos_sd_x
      slope_aspect_cos_sd <- terra::focal(slope_aspect_cos, full_mat, sd)
      names(slope_aspect_cos_sd) <- paste0("slope_aspect_cos_sd_", pred_rad)
      
      #---------------focal vegetation predictors, sd
      
      # ch_sd_x
      ch_sd <- terra::focal(ch, full_mat, sd)
      names(ch_sd) <- paste0("ch_sd_", pred_rad)
      
      # cc_sd_x
      cc_sd <- terra::focal(cc, full_mat, sd)
      names(cc_sd) <- paste0("cc_sd_", pred_rad)
      
      #---------------topographic position index (tpi)
      
      # tpi_x
      tpi <- elevation - terra::focal(x = elevation, w = annul_mat, fun = mean)
      names(tpi) <- paste0("tpi_", pred_rad)
      
      # focal radius stack
      stack.focal.radius <- c(elevation_mean, slope_mean, slope_derivative_mean, 
                              curvature_mean, curvature_plan_mean, curvature_prof_mean,
                              aspect_sin_mean, aspect_cos_mean, slope_aspect_sin_mean, slope_aspect_cos_mean,
                              ch_mean, cc_mean,
                              elevation_sd, slope_sd, slope_derivative_sd, 
                              curvature_sd, curvature_plan_sd, curvature_prof_sd,
                              aspect_sin_sd, aspect_cos_sd, slope_aspect_sin_sd, slope_aspect_cos_sd,
                              ch_sd, cc_sd,
                              tpi)
      
      # add to combined focal stack
      if (pred_rad == 2){
        stack.focal <- stack.focal.radius
      } else {
        stack.focal <- c(stack.focal, stack.focal.radius)
      }
      
    }
    
    # combine local and focal stacks
    stack <- c(stack.local, stack.focal)
    
    # save it to file, if selected    
    if (save == TRUE) {
      if (vi_type == "omnidir") {
        out_name <- "predictor_raster_stack_f360.tif"
        terra::writeRaster(stack, file.path(save_dir, out_name), overwrite = TRUE)
      } else if (vi_type == "directional_single") {
        out_name <- paste0("predictor_raster_stack_f", vi_fov, "_a", vi_azi, ".tif")
        terra::writeRaster(stack, file.path(save_dir, out_name), overwrite = TRUE)
      } 
    }
    
    # extract pixel values at points
    xy_only <- pts |> dplyr::select(c(x, y))
    vals <- terra::extract(stack, xy_only, method = "bilinear", bind = TRUE)
    df <- as.data.frame(vals, geom = "XY")
    final_df <- merge(df, pts)
    
    # return stack and df
    message(paste0(Sys.time(), " gen_preds() is complete"))
    return(list(pred_pts = final_df, pred_rast = stack))
    
    ##############################################################################
    # directional_random
    ##############################################################################
    
  } else if (vi_type == "directional_random"){
    
    # start by extracting all local variables the easy way
    xy_only <- pts |> dplyr::select(c(x,y))
    vals <- terra::extract(stack.local, xy_only, method = "bilinear", bind=TRUE)
    df <- as.data.frame(vals, geom="XY")
    final_df <- merge(df, pts)
    pts <- final_df
    
    # next make sure you have an azimuth column before jumping into the loop
    if ("azimuth" %in% colnames(pts) == FALSE){
      stop("\nInput data.frame must have an 'azimuth' column when using vi_type 'directional_random'.\n")
    }
    
    # for this one you have to do pt by pt for focal
    for (i in 1:nrow(pts)){
      
      # print an update every 25 points
      if(i %% 25 == 0 | i == 1){
        message(paste0(Sys.time(), "   Starting point ", i, "..."))
      }
      
      # isolate point and azimuth
      x <- pts[i, "x"]
      y <- pts[i, "y"]
      az <- pts[i, "azimuth"]
      
      # make it a point
      pt <- sf::st_point(c(x, y))
      pt_sv <- terra::vect(pt)
      terra::crs(pt_sv) <- terra::crs(dtm)
      
      # normalize aspect relative to viewing direction
      aspect_norm <- aspect - az
      
      # aspect_sin
      aspect_sin <- sin(aspect_norm * pi/180)
      names(aspect_sin) <- "aspect_sin"
      
      # aspect_cos
      aspect_cos <- cos(aspect_norm * pi/180)
      names(aspect_cos) <- "aspect_cos"
      
      # slope_aspect_sin
      slope_aspect_sin <- slope * aspect_sin
      names(slope_aspect_sin) <- "slope_aspect_sin"
      
      # slope_aspect_cos
      slope_aspect_cos <- slope * aspect_cos
      names(slope_aspect_cos) <- "slope_aspect_cos"
      
      # extract those values
      aspect_sin_val <- terra::extract(aspect_sin, pt_sv, method = "bilinear")
      aspect_cos_val <- terra::extract(aspect_cos, pt_sv, method = "bilinear")
      slope_aspect_sin_val <- terra::extract(slope_aspect_sin, pt_sv, method = "bilinear")
      slope_aspect_cos_val <- terra::extract(slope_aspect_cos, pt_sv, method = "bilinear")
      
      # add to df
      pts[i, "aspect_sin"] <- aspect_sin_val[1, 2]
      pts[i, "aspect_cos"] <- aspect_cos_val[1, 2]
      pts[i, "slope_aspect_sin"] <- slope_aspect_sin_val[1, 2]
      pts[i, "slope_aspect_cos"] <- slope_aspect_cos_val[1, 2]
      
      # create list of predictor rasts and names to loop through for focal extraction
      predlist <- list(names = c("elevation", "slope", "slope_derivative", 
                                 "curvature", "curvature_plan", "curvature_prof",
                                 "aspect_sin", "aspect_cos", "slope_aspect_sin", "slope_aspect_cos", 
                                 "ch", "cc"), 
                       rasts = list(elevation, slope, slope_derivative,
                                    curvature, curvature_plan, curvature_prof,
                                    aspect_sin, aspect_cos, slope_aspect_sin, slope_aspect_cos, 
                                    ch, cc))
      
      # loop through focal radii
      for (pred_rad in c(2,4,6,8,16,32)) {
        
        # convert pred_rad to meters using resolution and agg_fact
        conv <- terra::res(dtm)[1] * agg_fact
        
        # add small amount to radius
        buff_dist <- pred_rad * conv + 0.5 * conv
        
        # create two wedges, outer and inner
        outer_wdg <- VisiMod::wedge(x, y, buff_dist, az, vi_fov)
        terra::crs(outer_wdg) <- terra::crs(dtm)
        inner_buf <- terra::buffer(pt_sv, buff_dist/2)
        
        # create annulus wedge (donut difference between outer and inner wedges)
        annulus <- terra::erase(outer_wdg, inner_buf)
        
        # loop through predictors
        for (k in 1:length(predlist$names)) {
          
          # get SpatRaster and name
          pred_ras <- predlist$rasts[[k]]
          name <- predlist$names[k]
          
          # get focal mean and sd
          ext_mean <- terra::extract(pred_ras, outer_wdg, 
                                     mean, touches = T)
          ext_sd <- terra::extract(pred_ras, outer_wdg, 
                                   sd, touches = T)
          
          # get column names for mean and sd
          colname_mean <- paste0(name, "_mean_", pred_rad)
          colname_sd <- paste0(name, "_sd_", pred_rad)
          
          # add values to data.frame
          pts[i, colname_mean] <- ext_mean[1, 2]
          pts[i, colname_sd] <- ext_sd[1, 2]
          
        }
        
        # handle tpi separately
        wa_elev <- terra::extract(elevation, annulus, 
                                  mean, exact = TRUE)
        elev_val <- pts[i, "elevation"]
        tpi <- elev_val - wa_elev[1, 2]
        pred_col <- paste0("tpi_", pred_rad)
        pts[i, pred_col] <- tpi
        
      }
      
    }
    
    # rename your df and return it
    final_df <- pts
    message(paste0(Sys.time(), " gen_preds() is complete"))
    return(final_df)
    
  } 
  
}