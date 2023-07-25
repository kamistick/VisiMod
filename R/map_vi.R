#' Map VI across a landscape
#' 
#' This function takes in a model and a stack of predictor rasters and models VI across each pixel.
#' 
#' @param mod random forest model 
#' @param predictors a multiband raster of all predictors
#' @param percent_threads percent of threads that should be used in prediction, default is 50. 
#' 
#' @export
#' @examples
#' 
#' # get your dtm and dsm
#' dsm <- rast("dsm.tif")
#' dtm <- rast("dtm.tif")
#'
#' # get your points
#' my_points <- generate_pts(dtm, dsm, 100, 1000)
#'
#' # calculate vi
#' df <- calculate_vi(dtm, dsm, my_points, "omnidir", 500)
#'
#' # generate predictors
#' df_preds <- gen_preds(dtm, dsm, df, "omnidir", save = TRUE)
#'
#' # model
#' mod <- model_vi(df_preds, TRUE, FALSE, 80)
#'
#' rf <- mod[[2]]
#' predictors <- rast("predictor_raster_stack.tif"))
#' 
#' mymap <- map_vi(rf, predictors, 50)
#' plot(mymap)

map_vi <- function(mod, predictors, percent_threads){
  useThreads <- floor(detectCores() * (percent_cores/100))
  mod_names <- mod$forest$independent.variable.names
  pred_names <- names(predictors)
  if(all(pred_names %in% mod_names)==FALSE){
    stop("Independent variable names from model do not match predictor raster names. The names in the Raster object should exactly match those expected by the model.")
  }
  pred.fun <- function(model, ...){
    
    predict(model, ..., num.threads = useThreads)$predictions
  }
  
  pred_ras <- terra::predict(predictors, mod, fun = pred.fun, na.rm=TRUE)
  
  return(pred_ras)
}