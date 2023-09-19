#' Model and Map VI
#'
#' Full workflow for mapping VI using random forest modeling.  
#' @param dtm SpatRaster. A digital terrain model at finest resolution available
#' @param dsm SpatRaster. A digital surface model at finest resolution available
#' @param num_pts Numeric. Number of points used to train the model 
#' @param dist Numeric. distance to which VI is calculated. If multiple distances should be considered use c().
#' @param vi_type Character. Defines which type of VI calculation you would like to conduct. One of: "omnidir" (omnidirectional, 360 degree), "directional_single" (a single specified view direction for each point), "directional_random" (viewing directions will be randomly assigned to each point).
#' @param vi_fov Numeric. Defines the angular field of view, in degrees, of the directional wedge used for the VI calculation. Only used if vi_type == "directional_single". Values must be > 0 and < 360
#' @param vi_azi Numeric. Defines the azimuth, or central viewing direction, in degrees, of the directional wedge used for VI calculation.  Only used if vi_type == "directional_single". Values must be >= 0 and < 360.
#' @param save_dir Character. The directory where intermediary files will be saved.
#' @return A SpatRaster, or list of SpatRasters, of mapped VI with a potential range of 0-1.
#' 
#' @export
#' @examples
#' # get your dtm and dsm
#' dsm <- rast("dsm.tif")
#' dtm <- rast("dtm.tif")
#' 
#' # model and map VI
#' outras <- modMapVI(dtm, dsm, 1000, 500, "directional_single", 120, c(0,120,240),  "C:\\proj1\\")
#' 
#' # create a 3 band raster to display in RGB
#' stack1 <- terra::`add<-`(outras[[1]], outras[[2]])
#' multiband_ras <- terra::`add<-`(stack1, outras[[3]])
#' writeRaster(multiband_ras, "C:\\proj1\\vi_map.tif")


VisiMod <- function(dtm, dsm, num_pts, dist, vi_type, vi_fov, vi_azi, save_dir){
  
  message(paste0(Sys.time(), ": Generating points..."))
  gpt <- gen_pts(dtm, dsm, num_pts, dist)
  
  # now we split off b/c omni will be handled differently 
  if(vi_type=="omnidir"){
    # calculate_vi 
    message(paste0(Sys.time(), ": Calculating VI..."))
    cv_df <- calc_vi(dtm, dsm, gpt, vi_type, dist, vi_fov, vi_azi)
    
    message(paste0(Sys.time(), ": Generating predictors..."))
    gpd <- gen_preds(dtm, dsm, gpt, vi_type, vi_fov, vi_azi, save=TRUE, save_dir)
    
    message(paste0(Sys.time(), ": Modeling..."))
    df <- merge(gpd, cv_df, by=c("x", "y"))
    mod <- mod_vi(df, cross_validate = FALSE, tune = tune)
    
    message(paste0(Sys.time(), ": Mapping..."))
    preds <- rast(paste0(save_dir, "\\predictor_raster_stack_f360.tif"))
    vimap <- map_vi(mod, preds, 50)
    return(vimap)
    
  } else if (vi_type =="directional_single"){
    # i think i want to "initialize" my raster "stack"
    rass <- c()
    for (fov in vi_fov){
      for(az in vi_azi){
        # calculate_vi
        message(paste0(Sys.time(), ": Calculating VI..."))
        cv_df <- calc_vi(dtm, dsm, gpt, vi_type, dist, fov, az)
        
        message(paste0(Sys.time(), ": Generating predictors..."))
        gpd <- gen_preds(dtm, dsm, gpt, vi_type, fov, az, save=TRUE, save_dir)
        preds <- rast(paste0(save_dir, "\\predictor_raster_stack_f", as.character(fov), "_a", as.character(az), ".tif"))
        
        for (d in dist){
          message(paste0(Sys.time(), ": Modeling..."))
          dist_col <- paste0("vi_", as.character(d))
          cv_df_dist <- cv_df %>%
            dplyr::select(c("x", "y", dist_col))
          df <- merge(gpd, cv_df_dist, by=c("x", "y"))
          mod <- mod_vi(df, cross_validate = FALSE, tune = tune)
          message(paste0(Sys.time(), ": Mapping..."))
          vimap <- map_vi(mod, preds, 50)
          names(vimap) <- paste0("vi_d", as.character(d), "_f", as.character(fov), "_a", as.character(az))
          rass <- c(rass, vimap)
        }
      }
    }
    
    return(rass)
  } else {
    stop("vi_type must be 'omnidir' or 'single_directional.")
  }
  
}