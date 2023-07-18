#' Model VI using a random forest
#'
#' Models VI froma dataframe that includes predictor values.
#' @param df dataframe with predictor variables, target vi values, and an 'x' and 'y' column with coordinates
#' @param cross_validate default is FALSE. True will execute a 4 fold spatial cross validation. False will train on all available points. If true, the predictions that are returned are predictions for folds made by models build excluding the respective fold. For example, predictions made on fold one are generated from models built using folds 2, 3, and 4.
#' @param tune default is FALSE. TRUE will tune a full model and use tuned hyperparameters (mtry, minimum node size, sample fraction) for modeling.
#' @param percent_cores The percent of cores you'd like to use for modeling and tuning (default is 50).
#' @return dataframe with predictions column indicating predictions made by the model.
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
#' df <- calculate_vi(dtm, dsm, my_points, "omnidir", 500)
#'
#' # generate predictors
#' df_preds <- gen_preds(dtm, dsm, df, "omnidir", save = FALSE)
#'
#' # model
#' mod <- model_vi(df_preds, TRUE, FALSE, 80)
#'
#' # print the aggregate model results
#' print(mod[[2]])


model_vi <- function(df, cross_validate = FALSE, tune = FALSE, percent_cores = 50){
  # first thing check for column names
  contains_vi <- FALSE
  num_vi_cols <- 0
  contains_xy <- FALSE
  for (colname in colnames(df)){
    if(grepl("vi_", colname)){
      contains_vi <- TRUE
      num_vi_cols <- num_vi_cols + 1
    }
  }
  if ("x" %in% colnames(df) & "y" %in% colnames(df)){
    contains_xy <- TRUE
  }

  if(contains_vi==TRUE & contains_xy==TRUE){
    if (num_vi_cols>1){
      stop("Input dataframe must only contain one 'vi_x' column, please remove additional or format accordingly (where x = vi radius).")
    }
    print("A column name contains vi, using that as target variable")
  }
  if (contains_vi ==FALSE){
    stop("Input dataframe must contain a column with 'vi_x' to serve as the target for modeling. Please format accordinly (where x = vi radius).")
  }

  if(contains_xy ==FALSE & cross_validate==TRUE){
    stop("Input dataframe must contain separate 'x' and 'y' columns with coordinates for spatial cross-validation.")
  }

  # get your vi column name & rename for use in ranger formula
  vi_col_num <- grep("vi_", colnames(df), ignore.case = TRUE)
  colnames(df)[vi_col_num] <- "vi"

  # remove any nas
  df <- df %>%
    na.omit()

  # if azimuth is a column remove that too
  if ("azimuth" %in% colnames(df)){
    df <- df %>%
      dplyr::select(-c(azimuth))
  }

  # get your cores
  useCores <- floor(paralell::detectCores() * (percent_cores/100))
  if (useCores<1){
    stop("Too few cores selected, please increase the percentage of cores used (values between 0-100) must result in at least one core.")
  }

  if (cross_validate == FALSE){
    if(contains_xy==TRUE){
      df <- df %>%
        dplyr::select(-c(x,y))
    }
    rf_all <- ranger::ranger(formula = vi~., data= df, importance = "permutation", num.threads = useCores)
    if(tune ==TRUE){
      #tune by making task, tuning, and getting tune vars
      rf_task <- makeRegrTask(data = df, target = "vi")
      tuned <- tuneRanger(rf_task, num.threads = useCores)
      mtry_val <- tuned$recommended.pars$mtry
      min_node_val <- tuned$recommended.pars$min.node.size
      sample_fraction_val <- tuned$recommended.pars$sample.fraction

      rf <- ranger::ranger(formula = vi~., data= df, mtry = mtry_val, importance = "permutation", min.node.size = min_node_val, sample.fraction = sample_fraction_val, num.threads = useCores)
      return(rf)
    } else {
      #final rf is that one
      rf <- rf_all
      return(rf)
    }
  }
  if (cross_validate == TRUE){
    df_cv <- df # df_cv will retain fold info, df will not
    # get midpoints for spatial cv
    mid_x <- ((max(df$x) + min(df$x))/2)
    mid_y <- ((max(df$y) + min(df$y))/2)
    ## assign folds
    df_cv$fold <- 0
    ## fold 1: < x & > y (top left)
    df_cv$fold[df_cv$x < mid_x & df_cv$y > mid_y] <- "1"
    ## fold 2: > x & > y (top right)
    df_cv$fold[df_cv$x > mid_x & df_cv$y > mid_y] <- "2"
    ## fold 3: < x & < y (bottom left)
    df_cv$fold[df_cv$x < mid_x & df_cv$y < mid_y] <- "3"
    ## fold 4: > x & < y (bottom right)
    df_cv$fold[df_cv$x > mid_x & df_cv$y < mid_y] <- "4"

    # plot this to show whats going on
    p1 <- ggplot(data = df_cv, aes(x = x, y = y, color = fold)) +
      geom_point() + geom_hline(yintercept = mid_y) + geom_vline(xintercept = mid_x) +
      ggtitle("Folds used for spatial cross-validation")+
      theme(panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth=0.5),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(color = NA, fill = NA),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
    print(p1)

    if(tune==TRUE){
      # remove xy to start the tune
      df_noxy <- df %>%
        dplyr::select(-c(x,y))
      # model, tune, get hyperparams
      rf_all <- ranger::ranger(formula = vi~., data= df_noxy, importance = "permutation", num.threads = useCores)
      rf_task <- makeRegrTask(data = df, target = "vi")
      tuned <- tuneRanger(rf_task, num.threads = useCores)
      mtry_val <- tuned$recommended.pars$mtry
      min_node_val <- tuned$recommended.pars$min.node.size
      sample_fraction_val <- tuned$recommended.pars$sample.fraction

      # then do everything in the else statment but using hyperparams
      df_cv <- df_cv %>%
        dplyr::select(-c(x,y))
      fold1 <- df_cv[df_cv$fold=="1",]
      fold2 <- df_cv[df_cv$fold=="2",]
      fold3 <- df_cv[df_cv$fold=="3",]
      fold4 <- df_cv[df_cv$fold=="4",]

      # fold 1
      # 1: need to remove fold for modeling
      test1_df <- fold1 %>%
        dplyr::select(-c(fold)) %>%
        na.omit()
      train <- rbind(fold2, fold3, fold4)
      train_df <- train %>%
        dplyr::select(-c(fold)) %>%
        na.omit()

      rf1 <- ranger::ranger(formula = vi~., data= train_df, mtry = mtry_val, importance = "permutation", min.node.size = min_node_val, sample.fraction = sample_fraction_val, num.threads = useCores)
      predictions1 <- predict(rf1, test1_df, num.threads = useCores)
      test1_df$predictions <- predictions1$predictions

      # fold 2
      test2_df <- fold2 %>%
        dplyr::select(-c(fold)) %>%
        na.omit()
      train <- rbind(fold1, fold3, fold4)
      train_df <- train %>%
        dplyr::select(-c(fold)) %>%
        na.omit()

      rf2 <- ranger::ranger(formula = vi~., data= train_df, mtry = mtry_val, importance = "permutation", min.node.size = min_node_val, sample.fraction = sample_fraction_val, num.threads = useCores)
      predictions2 <- predict(rf2, test2_df, num.threads = useCores)
      test2_df$predictions <- predictions2$predictions

      # fold 3
      test3_df <- fold3 %>%
        dplyr::select(-c(fold)) %>%
        na.omit()
      train <- rbind(fold2, fold1, fold4)
      train_df <- train %>%
        dplyr::select(-c(fold)) %>%
        na.omit()

      rf3 <- ranger::ranger(formula = vi~., data= train_df, mtry = mtry_val, importance = "permutation", min.node.size = min_node_val, sample.fraction = sample_fraction_val, num.threads = useCores)
      predictions3 <- predict(rf3, test3_df, num.threads = useCores)
      test3_df$predictions <- predictions3$predictions

      # fold 4
      test4_df <- fold3 %>%
        dplyr::select(-c(fold)) %>%
        na.omit()
      train <- rbind(fold2, fold1, fold3)
      train_df <- train %>%
        dplyr::select(-c(fold)) %>%
        na.omit()

      rf4 <- ranger::ranger(formula = vi~., data= train_df, mtry = mtry_val, importance = "permutation", min.node.size = min_node_val, sample.fraction = sample_fraction_val, num.threads = useCores)
      predictions4 <- predict(rf4, test4_df, num.threads = useCores)
      test4_df$predictions <- predictions4$predictions

      # combine all your tests
      all_test <- rbind(test1_df,test2_df, test3_df, test4_df)

      # get your test stats
      lm <- lm(predictions ~ vi, all_test)
      rmse <- rmse(all_test[,"vi"], all_test[,"predictions"])
      nrmse <- rmse/(max(all_test[,"vi"])-min(all_test[,"vi"]))
      r2 <- summary(lm)$r.squared

      # run the full model too
      df_noxy <- df %>%
        dplyr::select(-c(x,y))
      rf_all <- ranger::ranger(formula = vi~., data= df, mtry = mtry_val, importance = "permutation", min.node.size = min_node_val, sample.fraction = sample_fraction_val, num.threads = useCores)

      # return df with predictions, R2 and RMSE of test preds, and full rf
      prnt_stmt <- paste0("The cross-validated, tuned model had an R squared of ", as.character(round(r2,2)),
                          "and an nRMSE of ", as.character(round(nrmse,2)), ".",
                          " Based on the test predictions from each quadrant.")
      print(prnt_stmt)
      return(list(all_test, rf_all))


    } else {
      # just do your 4 folds
      # remove xy to start the modeling
      df_cv <- df_cv %>%
        dplyr::select(-c(x,y))
      fold1 <- df_cv[df_cv$fold=="1",]
      fold2 <- df_cv[df_cv$fold=="2",]
      fold3 <- df_cv[df_cv$fold=="3",]
      fold4 <- df_cv[df_cv$fold=="4",]

      # fold 1
      # 1: need to remove fold for modeling
      test1_df <- fold1 %>%
        dplyr::select(-c(fold)) %>%
        na.omit()
      train <- rbind(fold2, fold3, fold4)
      train_df <- train %>%
        dplyr::select(-c(fold)) %>%
        na.omit()

      rf1 <- ranger::ranger(formula = vi~., data= train_df, importance = "permutation", num.threads = useCores)
      predictions1 <- predict(rf1, test1_df, num.threads = useCores)
      test1_df$predictions <- predictions1$predictions

      # fold 2
      test2_df <- fold2 %>%
        dplyr::select(-c(fold)) %>%
        na.omit()
      train <- rbind(fold1, fold3, fold4)
      train_df <- train %>%
        dplyr::select(-c(fold)) %>%
        na.omit()

      rf2 <- ranger::ranger(formula = vi~., data= train_df, importance = "permutation", num.threads = useCores)
      predictions2 <- predict(rf2, test2_df, num.threads = useCores)
      test2_df$predictions <- predictions2$predictions

      # fold 3
      test3_df <- fold3 %>%
        dplyr::select(-c(fold)) %>%
        na.omit()
      train <- rbind(fold2, fold1, fold4)
      train_df <- train %>%
        dplyr::select(-c(fold)) %>%
        na.omit()

      rf3 <- ranger::ranger(formula = vi~., data= train_df, importance = "permutation", num.threads = useCores)
      predictions3 <- predict(rf3, test3_df, num.threads = useCores)
      test3_df$predictions <- predictions3$predictions

      # fold 4
      test4_df <- fold3 %>%
        dplyr::select(-c(fold)) %>%
        na.omit()
      train <- rbind(fold2, fold1, fold3)
      train_df <- train %>%
        dplyr::select(-c(fold)) %>%
        na.omit()

      rf4 <- ranger::ranger(formula = vi~., data= train_df, importance = "permutation", num.threads = useCores)
      predictions4 <- predict(rf4, test4_df, num.threads = useCores)
      test4_df$predictions <- predictions4$predictions

      # combine all your tests
      all_test <- rbind(test1_df,test2_df, test3_df, test4_df)

      # get your test stats
      lm <- lm(predictions ~ vi, all_test)
      rmse <- rmse(all_test[,"vi"], all_test[,"predictions"])
      nrmse <- rmse/(max(all_test[,"vi"])-min(all_test[,"vi"]))
      r2 <- summary(lm)$r.squared

      # run the full model too
      df_noxy <- df %>%
        dplyr::select(-c(x,y, fold))
      rf_all <- ranger::ranger(formula = vi~., data= df_noxy, importance = "permutation", num.threads = useCores)

      # return df with predictions, R2 and RMSE of test preds, and full rf
      prnt_stmt <- paste0("The cross-validated model had an R squared of ", as.character(round(r2,2)),
                          " and an nRMSE of ", as.character(round(nrmse,2)), ". ",
                          " Based on the test predictions from each quadrant.")
      print(prnt_stmt)
      return(list(all_test, rf_all))
    }

  }

}
