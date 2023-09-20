#' Map VI across a landscape
#' 
#' The sixth and final suggested function in the VisiMod workflow, following (1) [VisiMod::prep_dems()] (2) [VisiMod::gen_pts()], (3) [VisiMod::calc_vi()], (4) [VisiMod::gen_preds()], and (5) [VisiMod::mod_vi()]. This function takes a random forest model (generated using [VisiMod::mod_vi()]), representing the predictive relationship between visibility index (VI) (calculated at a series of point locations using [VisiMod::calc_vi()]), and a suite of visibility predictor layers representing terrain and vegetation structure (derived using [VisiMod::gen_preds()]), and applies this model to the prediction of VI across an entire study area of interest. Successful execution of this function will result in a spatially exhaustive ("wall-to-wall") map, where each pixel represents an estimate of VI based on a spatial scope of interest (viewing distance, and direction and field of view, if mapping directional VI).
#' @details
#' * The spatial resolution, coordinate system, origin, and extent of the resulting map will be based entirely on those of your multiband predictors raster.
#' * It is recommended to be somewhat conservative with your definition of `cores`. Using a large proportion of your computer's cores can cause a memory overload, as a large amount of information has to be stored in memory for each core.
#' * Note that, by default, if `save == TRUE` and you supply an existing filename to `out_file`, it will be overwritten.
#' 
#' @param mod ranger. A random forest model generated using [ranger::ranger()] aimed at predicting VI.
#' @param predictors SpatRaster. A multiband SpatRaster of all predictors used for VI modeling.
#' @param cores Numeric. Defines the number of cores you would like to use for parallel processing. Defaults to half of the cores on your machine.
#' @param save Boolean. Defines whether or not you would like to save the resulting VI prediction map to file. Default is TRUE.
#' @param out_file Character. The full file path (including the ".tif" extension) of your output VI map.
#' @return A SpatRaster containing pixel-level VI predictions. 
#' @export
#' @examples
#' 
#' # get your dtm and dsm
#' dsm <- rast("dsm.tif")
#' dtm <- rast("dtm.tif")
#' 
#' # check dtm and dsm
#' pd <- prep_dems(dtm, dsm, "C:/temp/dtm_filled.tif", "C:/temp/dsm_filled.tif")
#' 
#' # get your points
#' my_points <- generate_pts(dtm, dsm, 100, 1000)
#' 
#' # calculate vi
#' my_points <- calculate_vi(dtm, dsm, my_points, "directional_single", c(500, 1000), 90, 90, 4L, 5L)
#' 
#' # generate predictors
#' preds <- gen_preds(dtm, dsm, my_points, "directional_single", 90, 90, T, "C:/temp")
#' 
#' # model vi
#' mod <- model_vi(preds, 500, T, T, 5L)
#' 
#' # map vi
#' rf <- mod$ranger_mod
#' predictors <- preds$pred_rast
#' mymap <- map_vi(rf, predictors, 5L, T, "C:/temp/vi_map.tif")
#' plot(mymap)

map_vi <- function(mod, predictors, cores = floor(parallel::detectCores()/2), 
                   save = T, out_file = file.path(getwd(), "vi.tif")){
  
  # print message
  message(paste0(Sys.time(), " map_vi() has begun"))
  
  # check to make sure all predictor variables from the model are present within the predictor raster
  mod_names <- mod$forest$independent.variable.names
  pred_names <- names(predictors)
  if(all(pred_names %in% mod_names)==FALSE){
    stop("Independent variable names from model do not match predictor raster names. The names in the SpatRaster object should exactly match those expected by the model.")
  }
  
  # define prediction function
  pred.fun <- function(model, ...){
    predict(model, ..., num.threads = cores)$predictions
  }
  
  # generate prediction map
  message(paste0(Sys.time(), "   Generating prediction map..."))
  pred_ras <- terra::predict(predictors, mod, fun = pred.fun, na.rm=TRUE)
  
  # save, if desired
  if (save == T){
    message(paste0(Sys.time(), "   Writing to file..."))
    terra::writeRaster(pred_ras, out_file, overwrite = T)
  }
  
  # print message
  message(paste0(Sys.time(), " map_vi() is complete"))
  
  # return the prediction map
  return(pred_ras)
  
}