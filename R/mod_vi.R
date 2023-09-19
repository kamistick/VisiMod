#' Model VI using a random forest
#'
#' The fifth suggested function in the VisiMod workflow, folowing (1) `prep_dems()`, (2) `gen_pts()`, (3) `calculate_vi()`, and (4) `gen_preds()`. This function models the visibility index (VI) using random forests, driven by a data.frame containing sample points with known VI values (response variable), derived from viewshed analysis, and a suite of point-level predictor variables aimed at capturing the dominant landscape controls of visibility in a wildland (undeveloped) environment. It relies on the efficient, parallelized implementation of random forests from the ranger library, and provides the capacity to undergo an automatic tuning process provided by the tuneRanger library. The user can specify whether to perform a 4-fold spatial cross-validation procedure, or to simply use all of the input data to train the model.
#' @details
#' * The model will use all columns in `df`, except 'x' and 'y', to build a predictive model. So, if your data.frame contains additional columns (i.e., if you did not strictly follow the recommended VisiMod function sequence to arrive at a set of sample points), you should either (a) know that these columns will be included as potential predictors of VI; or (b) remove them prior to running the function
#' * If you opted to generate VI estimates at multiple viewing distances in the `calculate_vi()` step, your data.frame will likely have multiple 'vi_x' columns, where x is equal to one or more 'vi_rad' values. However, this function can only be used to model one viewing radius at a time, that you specify with `vi_rad`. If other VI columns are present, they will be automatically removed prior to modeling.
#' * Note that the cross-validation procedure defined by `cross_validate` uses a simple extent-based approach to split the study area up into four quadrants. It takes the study area's extent, finds the vertical and horizontal midpoints of that extent, creates four equally sized rectangles, and determines which points fall in each rectangle. If your study area is rectangular, you can expect that a similar number of points will fall in each of the four spatial folds, with minor variation occurring due to the assumed randomness of the initial point placement driven by `gen_pts()`. However, if your study area is not rectangular, it is possible that the points may be very unevenly distributed throughout the folds. It is even possible that one or more folds could contain no points, in extreme cases. In situations like this, if cross-validation is desired, it is recommended to carry out your own cross-validation procedure.
#' * The `ranger()` function upon which this modeling procedure is based will not run successfully with the presence of NA values in either the response variable or any of the predictor variables. By default, `model_vi()` removes rows with NAs prior to modeling. Users should either (a) be sure that no NAs are present within `df` prior to modeling; or (b) understand that the number of sample points used in the modeling procedure may be less than the total number of input points if NAs are present.
#' 
#' @param df data.frame. Should contain 'x' and 'y' columns (resulting from `gen_pts()`), a 'vi' column (resulting from `calculate_vi()`), and several predictor columns (resulting from `gen_preds()`). See Details for information on the presence of additional columns in your data.frame.
#' @param vi_rad Numeric. Defines the radius at which you intend to model VI. This value should match at least one of those used in the execution of `calculate_vi()`, and as a result, there should be a column in `df` named `vi_x` where `x == vi_rad`.
#' @param cross_validate Boolean. Defines whether or not you would like to perform a 4-fold spatial cross-validation procedure to assess model performance. If TRUE, the sample points are split into four equally-sized, extent-based spatial quadrants, and models are iteratively trained using points from 3/4 of the quadrants and applied to the prediction of the fourth. If FALSE, all sample points are used to train the model.
#' @param tune Boolean. Defines whether or not you would like to tune the model's hyperparameters (mtry, min.node.size, and sample.fraction) using the `tuneRanger()` function.
#' @param num_cores Numeric. Defines the number of cores you would like to use for parallel processing. Defaults to half of the cores on your machine.
#' @return A list with three items: (1) `df_pred_obs` is a data.frame containing model predictions vs true VI observations; (2) `ranger_mod` is a ranger model generated from the full dataset; and (3) `perf_mets` is a list containing basic model performance metrics comparing the predicted vs. observed VI values.
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
#' my_points <- generate_pts(dtm, dsm, 100, 1000)
#' 
#' # calculate vi
#' my_points <- calculate_vi(dtm, dsm, my_points, "directional_single", c(500, 1000), 90, 90, 4L, 5L)
#' 
#' # generate predictors
#' preds <- gen_preds(dtm, dsm, my_points, "directional_single", 90, 90, T, "C:/temp")
#' 
#' # model
#' mod <- model_vi(preds, 500, T, T, 5L)

mod_vi <- function(df, vi_rad, cross_validate = FALSE, tune = FALSE, num_cores = floor(parallel::detectCores()/2)){
  
  # print message
  message(paste0(Sys.time(), " model_vi() has begun"))
  
  # check for presence of x and y columns -- needed if cross_validate == T
  contains_xy <- FALSE
  if ("x" %in% colnames(df) & "y" %in% colnames(df)) {
    contains_xy <- TRUE
  }
  if (contains_xy == FALSE & cross_validate == TRUE) {
    stop("Input data.frame must contain separate 'x' and 'y' columns with coordinates for spatial cross-validation.")
  }
  
  # check for presence of distance-specific vi column defined by vi_rad
  vi_col <- paste0("vi_", vi_rad)
  if (!vi_col %in% colnames(df)){
    stop("A vi column at the specified vi_rad ('vi_", vi_rad, "') does not exist within df")
  }
  
  # remove other vi columns, if present, to avoid including them in the model
  oth_vi_cols <- colnames(df)[grep("vi_", colnames(df))]
  oth_vi_cols <- oth_vi_cols[oth_vi_cols != vi_col]
  if (length(oth_vi_cols) > 0){
    keep_cols <- colnames(df)[!colnames(df) %in% oth_vi_cols]
    df <- df[,keep_cols]
  }
  
  # if azimuth is a column remove that too
  if ("azimuth" %in% colnames(df)){
    df <- df |>
      dplyr::select(-c(azimuth))
  }
  
  # get your vi column name & rename for use in ranger formula
  colnames(df)[colnames(df) == vi_col] <- "vi"

  # remove any nas
  df <- df |> na.omit()

  #----------------------modeling without cross-validation
  
  if (cross_validate == FALSE){
    
    # print message
    message(paste0(Sys.time(), " modeling without cross-validation"))
    
    # remove x and y columns, since they're not needed without cv
    if(contains_xy == TRUE) df_noxy <- df |> dplyr::select(-c(x,y))
    
    #-------------------with tuning
    
    if(tune == TRUE){
      
      # tune by making task, tuning, and getting tune vars
      message(paste0(Sys.time(), " tuning the model"))
      rf_task <- mlr::makeRegrTask(data = df_noxy, target = "vi")
      tuned <- tuneRanger::tuneRanger(rf_task, num.threads = num_cores,
                                      show.info = getOption("mlrMBO.show.info", F))
      mtry_val <- tuned$recommended.pars$mtry
      min_node_val <- tuned$recommended.pars$min.node.size
      sample_fraction_val <- tuned$recommended.pars$sample.fraction
      
      # build full model
      message(paste0(Sys.time(), " building the full model"))
      rf_all <- ranger::ranger(formula = vi ~ ., data = df_noxy, 
                               mtry = mtry_val, importance = "permutation", 
                               min.node.size = min_node_val, sample.fraction = sample_fraction_val, 
                               num.threads = num_cores)
    
    #-------------------without tuning

    } else {
      
      # build full model
      message(paste0(Sys.time(), " building the full model"))
      rf_all <- ranger::ranger(formula = vi ~ ., data = df_noxy, 
                               importance = "permutation", num.threads = num_cores)
      
    }
    
    # compare predictions vs. observations
    message(paste0(Sys.time(), " assessing model performance"))
    predictions <- rf_all$predictions
    df_pred_obs <- data.frame(observed = df_noxy$vi,
                              predicted = predictions)
    
    # get performance metrics and put them in a list
    lm <- lm(predicted ~ observed, df_pred_obs)
    rmse <- Metrics::rmse(df_pred_obs[, "observed"], df_pred_obs[, "predicted"])
    nrmse <- rmse/(max(df_pred_obs[, "observed"]) - min(df_pred_obs[, "predicted"]))
    r2 <- summary(lm)$adj.r.squared
    perf_mets <- list(r2 = r2, rmse = rmse, nrmse = nrmse)
    
    # return pred vs. obs, the full model, and the performance metrics
    message(paste0(Sys.time(), " model_vi() is complete"))
    return(list(df_pred_obs = df_pred_obs, ranger_mod = rf_all, perf_mets = perf_mets))
    
  }
  
  #----------------------modeling with cross-validation
  
  if (cross_validate == TRUE){
    
    # print message
    message(paste0(Sys.time(), " modeling with cross-validation"))
    
    # create copy of df (df_cv will retain fold info, df will not)
    df_cv <- df 
    
    # get midpoints for spatial cv
    mid_x <- ((max(df$x) + min(df$x))/2)
    mid_y <- ((max(df$y) + min(df$y))/2)
    
    # assign folds
    df_cv$fold <- NA
    df_cv$fold[df_cv$x < mid_x & df_cv$y >= mid_y] <- 1
    df_cv$fold[df_cv$x >= mid_x & df_cv$y >= mid_y] <- 2
    df_cv$fold[df_cv$x < mid_x & df_cv$y < mid_y] <- 3
    df_cv$fold[df_cv$x >= mid_x & df_cv$y < mid_y] <- 4

    # plot the points to show spatial folds
    par(mar = c(5,5,1,1), las = 1)
    plot(df_cv$x / 1000, df_cv$y / 1000, pch = 16, 
         col = df_cv$fold, xlab = "x (km)", ylab = "y (km)")
    grid()
    abline(h = mid_y / 1000)
    abline(v = mid_x / 1000)
    legend("topleft", "Fold 1", bty = "n", cex = 1.5, x.intersp = 0)
    legend("topright", "Fold 2", bty = "n", cex = 1.5, x.intersp = 0)
    legend("bottomleft", "Fold 3", bty = "n", cex = 1.5, x.intersp = 0)
    legend("bottomright", "Fold 4", bty = "n", cex = 1.5, x.intersp = 0)
  
    # remove x and y cols and rows with nas
    df_noxy <- df |> dplyr::select(-c(x, y))
    df_cv <- df_cv |> dplyr::select(-c(x, y))
    
    #-------------------full model with tuning
    
    if(tune==TRUE){
      
      # tune by making task, tuning, and getting tune vars
      message(paste0(Sys.time(), " tuning the model"))
      rf_task <- mlr::makeRegrTask(data = df_cv, target = "vi")
      tuned <- tuneRanger::tuneRanger(rf_task, num.threads = num_cores,
                                      show.info = getOption("mlrMBO.show.info", F))
      mtry_val <- tuned$recommended.pars$mtry
      min_node_val <- tuned$recommended.pars$min.node.size
      sample_fraction_val <- tuned$recommended.pars$sample.fraction
      
      # build full model
      message(paste0(Sys.time(), " building the full model"))
      rf_all <- ranger::ranger(formula = vi ~ ., data = df_noxy, 
                               mtry = mtry_val, importance = "permutation", 
                               min.node.size = min_node_val, sample.fraction = sample_fraction_val, 
                               num.threads = num_cores)
      
    #-------------------full model without tuning

    } else {
      
      # build full model
      message(paste0(Sys.time(), " building the full model"))
      rf_all <- ranger::ranger(formula = vi ~ ., data = df_noxy, 
                               importance = "permutation", num.threads = num_cores)
      
    }
    
    # loop through folds for cv
    for (fold in seq(1,4)){
      
      # print message
      message(paste0(Sys.time(), " modeling fold #", fold))
      
      # define training and validation data
      train_df <- df_cv[df_cv$fold != fold,]
      valid_df <- df_cv[df_cv$fold == fold,]
      
      #-----------------cv model with tuning
      
      if (tune == TRUE){
        
        # build model
        rf <- ranger::ranger(formula = vi ~ ., data = train_df, 
                             mtry = mtry_val, importance = "permutation", 
                             min.node.size = min_node_val, sample.fraction = sample_fraction_val, 
                             num.threads = num_cores)
        
        # compare predictions vs. observations
        predictions <- predict(rf, valid_df, num.threads = num_cores)$predictions
        df_pred_obs_fold <- data.frame(observed = valid_df$vi,
                                       predicted = predictions)
        
      #-----------------cv model without tuning
        
      } else {
        
        # build model
        rf <- ranger::ranger(formula = vi ~ ., data = train_df, 
                             importance = "permutation", num.threads = num_cores)
        predictions <- predict(rf, valid_df, num.threads = num_cores)$predictions

        # compare predictions vs. observations
        df_pred_obs_fold <- data.frame(observed = valid_df$vi,
                                       predicted = predictions)
        
      }
      
      # compile predictions vs. observations among folds
      if (fold == 1){
        
        df_pred_obs <- df_pred_obs_fold
        
      } else {
        
        df_pred_obs <- rbind(df_pred_obs, df_pred_obs_fold)
        
      }
      
    }
    
    # get performance metrics and put them in a list
    message(paste0(Sys.time(), " assessing model performance"))
    lm <- lm(predicted ~ observed, df_pred_obs)
    rmse <- Metrics::rmse(df_pred_obs[, "observed"], df_pred_obs[, "predicted"])
    nrmse <- rmse/(max(df_pred_obs[, "observed"]) - min(df_pred_obs[, "predicted"]))
    r2 <- summary(lm)$adj.r.squared
    perf_mets <- list(r2 = r2, rmse = rmse, nrmse = nrmse)
    
    # return pred vs. obs, the full model, and the performance metrics
    message(paste0(Sys.time(), " model_vi() is complete"))
    return(list(df_pred_obs = df_pred_obs, ranger_mod = rf_all, perf_mets = perf_mets))
  
  }

}