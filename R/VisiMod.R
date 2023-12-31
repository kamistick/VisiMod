#' Model and Map VI 
#'
#' This function takes a digital terrain model (DTM) and a digital surface model (DSM) and builds a predictive model to map modeled visibility index (VI) across the study area. This function sequentially combines the full suite of VisiMod functions ((1) [VisiMod::prep_dems()] (2) [VisiMod::gen_pts()], (3) [VisiMod::calc_vi()], (4) [VisiMod::gen_preds()], (5) [VisiMod::mod_vi()], and [VisiMod::map_vi()]) to return a wall-to-wall SpatRaster of VI values.       
#' 
#' @details
#' * This function assumes default parameters for some of the internal functions ([VisiMod::gen_preds()] and [VisiMod::mod_vi()]). For [VisiMod::gen_preds()] a default aggregation factor of 10 is assumed. For [VisiMod::mod_vi()] `cross_validate` and `tune` are assumed to be `FALSE`. See respective function details for more. If users desire more control over these parameters they should execute the workflow in the order laid out in this function's description. 
#' * `dtm` and `dsm` SpatRasters can be defined using the terra library. They should have the same coordinate system, resolution, extent, and origin.
#' * The more points (`num_pts`) generated, the more robust the modeling procedure will be; however, more points will also increase processing time for subsequent functions in the workflow. The default is set to 200, which should balance model performance and processing time. We do not recommend generating more than 1000 points.
#' * Processing time will increase exponentially with increasing distances (`dist`). However, this also depends on the spatial resolution of the input `dtm`/`dsm`. For example, looking at 200m with a 1m resolution is functionally the same as a 400m distance with a 2m resolution, in terms of processing time. We do not recommend attempting this workflow at distances beyond 2000x the input resolution.
#' * Note that this function is parallelized and can leverage as many cores as your computer has available to speed up processing. As with all parallel processing in R, however, there is an overhead cost associated with setting up parallel operations. So, for small numbers of input points (`num_pts`) and/or short viewing distance (`dist`), using many cores may not speed up your processing significantly.
#' 
#' 
#' @param dtm SpatRaster. A digital terrain model at finest resolution available
#' @param dsm SpatRaster. A digital surface model at finest resolution available
#' @param num_pts Numeric. Number of points used to train the model 
#' @param dist Numeric. distance to which VI is calculated. If multiple distances should be considered use c().
#' @param vi_type Character. Defines which type of VI calculation you would like to conduct. One of: "omnidir" (omnidirectional, 360 degree), "directional_single" (a single specified view direction for each point), "directional_random" (viewing directions will be randomly assigned to each point).
#' @param vi_fov Numeric. Defines the angular field of view, in degrees, of the directional wedge used for the VI calculation. Only used if vi_type == "directional_single". Values must be > 0 and < 360
#' @param vi_azi Numeric. Defines the azimuth, or central viewing direction, in degrees, of the directional wedge used for VI calculation.  Only used if vi_type == "directional_single". Values must be >= 0 and < 360.
#' @param save_dir Character. The directory where intermediary and output files will be saved. Default is your working directory.
#' @param cores Numeric. The number of cores used for parallel processing of VI calculation, modeling, and mapping. The default number of cores is half of the cores on your machine.
#' @return Returns a list of SpatRasters of mapped VI with a potential range of 0-1. If only one combination of distance, fov, and azimuth are considered, the list will have a length of 1 and the map of VI can be accessed using VisiMod_fxn_output[[1]] where VisiMod_fxn_output <- VisiMod(...).
#' 
#' @export
#' @examples
#' # get your dtm and dsm
#' dsm <- rast("dsm.tif")
#' dtm <- rast("dtm.tif")
#' 
#' # model and map VI
#' outras <- VisiMod(dtm, dsm, 500, 500, "directional_single", 120, c(0,120,240),  "C:\\proj1\\", 4)
#' 
#' # create a 3 band raster to display in RGB
#' stack <- c(outras[[1]], outras[[2]], outras[[3]])
#' writeRaster(stack, "C:\\proj1\\vi_map.tif")


VisiMod <- function(dtm, dsm, num_pts, dist, vi_type, vi_fov=180, vi_azi=0, save_dir = getwd(), cores = floor(parallel::detectCores()/2)){
  
  # create clean time printing function
  prttm <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  message(paste0(prttm(), " Preparing rasters..."))
  pdems <- prep_dems(dtm, dsm, paste0(save_dir, "dtm_filled.tif"), paste0(save_dir, "dsm_filled.tif"))
  dtm <- pdems$dtm
  dsm <- pdems$dsm 
  
  message(paste0(prttm(), " Generating points..."))
  gpt <- gen_pts(dtm, dsm, num_pts, dist)
  
  # now we split off b/c omni will be handled differently 
  if(vi_type=="omnidir"){
    # calculate_vi 
    message(paste0(prttm(), " Calculating VI..."))
    cv_df <- calc_vi(dtm, dsm, gpt, vi_type, dist, vi_fov, vi_azi, cores)
    
    message(paste0(prttm(), " Generating predictors..."))
    gpd <- gen_preds(dtm, dsm, gpt, vi_type, vi_fov, vi_azi, 10L, save=TRUE, save_dir)
    #preds <- rast(paste0(save_dir, "\\predictor_raster_stack_f360.tif"))
    
    rass <- c()
    for (d in dist){
      message(paste0(prttm(), " Modeling..."))
      dist_col <- paste0("vi_", d)
      cv_df_dist <- cv_df %>%
        dplyr::select(c("x", "y", dist_col))
      df <- merge(gpd$pred_pts, cv_df_dist, by=c("x", "y"))
      mod <- mod_vi(df, d, cross_validate = FALSE, tune = F, cores)
      message(paste0(prttm(), " Mapping..."))
      vimap <- map_vi(mod$ranger_mod, gpd$pred_rast, cores, T, 
                      file.path(save_dir, paste0("vi_d", d, "_f360.tif")))
      names(vimap) <- paste0("vi_d", d, "_f360")
      rass <- c(rass, vimap)
    }
    
    return(rass)
    
  } else if (vi_type =="directional_single"){
    rass <- c()
    for (fov in vi_fov){
      for(az in vi_azi){
        # calculate_vi
        message(paste0(prttm(), " Calculating VI..."))
        cv_df <- calc_vi(dtm, dsm, gpt, vi_type, dist, fov, az, cores)
        
        message(paste0(prttm(), " Generating predictors..."))
        gpd <- gen_preds(dtm, dsm, gpt, vi_type, fov, az, 10L, save=TRUE, save_dir)
        
        for (d in dist){
          message(paste0(prttm(), " Modeling..."))
          dist_col <- paste0("vi_", d)
          cv_df_dist <- cv_df %>%
            dplyr::select(c("x", "y", dist_col))
          df <- merge(gpd$pred_pts, cv_df_dist, by=c("x", "y"))
          mod <- mod_vi(df, d, cross_validate = FALSE, tune = F, cores)
          message(paste0(prttm(), " Mapping..."))
          vimap <- map_vi(mod$ranger_mod, gpd$pred_rast, cores, T, 
                          file.path(save_dir, paste0("vi_d", d, "_f", fov, "_a", az, ".tif")))
          names(vimap) <- paste0("vi_d", d, "_f", fov, "_a", az)
          rass <- c(rass, vimap)
        }
      }
    }
    
    return(rass)
  } else {
    stop("vi_type must be 'omnidir' or 'single_directional.")
  }
  
}