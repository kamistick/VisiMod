#' Generate predictors for modeling visibility index
#'
#' Extracts predictor values at points. Predictors are generated at a resolution 10x the input raster resolution. If using 1 m input resolution, predictors will be generated at 10 m resolution. Predictors include: canopy cover, canopy height, elevation, slope, slope_derivative, curvature (plan and profile), aspect cosine, aspect sine, slope aspect cosine, slope aspect sine, as well as mean and standard deviation predictors  - for each of the the previous predictors - calculated within a focal area with radius x pixels for x in (2, 4, 6, 8, 16, 32)
#' @param dtm SpatRaster digital terrain model at finest resolution available
#' @param dsm SpatRaster digital surface model at finest resolution available
#' @param pts dataframe of your points, must include x, y columns, may include vi columns
#' @param vi_type which type of VI calculation did you do? One of: 'omnidir' (omnidirectional, 360 degree), 'directional_single' (a single specified view direction), 'directional_random' (multiple random view directions will be sampled)
#' @param vi_fov optional argument, with a default value of 180, if vi_type "directional_single" or "directional_random" are selected. Specifies the field of view for the vi calculation
#' @param vi_azi optional argument, with a default value of 0, if vi_type "directional_single" specify the direction of vi in degrees values between 0 and 359 (e.g. North = 0, South = 180, East = 90, West = 270)
#' @param save  Default is TRUE. If TRUE and if 'omnidir' or 'single_directional' are used then the multiband raster of predictors will also be saved.
#' @param save_dir if save = TRUE specify the directory where the multiband raster of predictors (ONLY if using 'omnidir' or 'directional_single') should be saved, default is your working directory.
#' @return A dataframe with predictor columns. If save = TRUE a .tif will also be written to file.
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
#' # generate predictors
#' preds <- gen_preds(dtm, dsm, my_points, "omnidir", save = FALSE)


gen_preds <- function(dtm, dsm, pts, vi_type, vi_fov=180, vi_azi=0, save=TRUE, save_dir=getwd()){
  print(paste0("Starting predictor generation at ", Sys.time()))

  # next get your chm and aggregate rasters
  chm <- dsm-dtm
  chm[chm < 0] <- 0

  ## start by aggregating dtm and chm to 10
  dtm10 <- terra::aggregate(dtm, 10, "mean")
  chm10 <- terra::aggregate(chm, 10, "mean")

  # generate predictor rasters
  # Slope, slope derivative, elevation, canopy cover, canopy height, curvature, aspect
  elevation <- dtm10
  slope <- terra::terrain(dtm10, "slope")
  slope_derivative <- terra::terrain(slope, "slope")
  aspect <- terra::terrain(dtm10, "aspect")
  ch <- chm10
  curvature <- spatialEco::curvature(dtm10, "total")
  curvature_plan <- spatialEco::curvature(dtm10, "planform")
  curvature_prof <- spatialEco::curvature(dtm10, "profile")

  # canopy cover requires multiple steps
  chm1 <- chm
  chm1[chm1 <= 2] <- 0
  chm1[chm1 > 2] <- 1
  cc <- terra::aggregate(chm1, 10, "mean")

  # start gathering all your rasters
  stack <- c(cc, ch, elevation, slope, slope_derivative, curvature, curvature_plan, curvature_prof)
  names(stack) <- c("cc", "ch", "elevation", "slope", "slope_derivative", "curvature", "curvature_plan", "curvature_prof")

  print(paste0("Local predictors generated at  ", Sys.time()))

  if(vi_type =="omnidir"){
    # deal with aspect first b/c it is unique to vi type
    aspect_sin <- sin(aspect*pi/180)
    aspect_cos <- cos(aspect*pi/180)
    slope_aspect_sin <- slope * aspect_sin
    slope_aspect_cos <- slope * aspect_cos
    aspect_stack <- c(aspect_sin, aspect_cos, slope_aspect_sin, slope_aspect_cos)
    names(aspect_stack) <- c("aspect_sin", "aspect_cos", "slope_aspect_sin", "slope_aspect_cos")
    stack <- c(stack, aspect_stack)
  }
  if(vi_type=="directional_single"){
    # deal with aspect first b/c it is unique to vi type
    aspect <- aspect - vi_azi
    aspect_sin <- sin(aspect*pi/180)
    aspect_cos <- cos(aspect*pi/180)
    slope_aspect_sin <- slope * aspect_sin
    slope_aspect_cos <- slope * aspect_cos
    aspect_stack <- c(aspect_sin, aspect_cos, slope_aspect_sin, slope_aspect_cos)
    names(aspect_stack) <- c("aspect_sin", "aspect_cos", "slope_aspect_sin", "slope_aspect_cos")
    stack <- c(stack, aspect_stack)
  }

  if (vi_type == "omnidir" | vi_type =="directional_single"){
    for(pred_rad in c(20,40,60,80,160,320)){
      print(paste0("Running focal predictors at ", as.character(pred_rad), " radius... ", Sys.time()))

      # for omnidirectional your focal matrix is a circle
      if (vi_type == "omnidir"){

        # define the focal region as a CIRCULAR BUFFER
        pred_rad_foc <- pred_rad/10
        buf <- terra::buffer(vect(cbind(0,0), crs=crs(dtm)), pred_rad_foc)
        r <- rast(xmin=-pred_rad_foc, xmax = pred_rad_foc, ymax = pred_rad_foc, ymin= -pred_rad_foc, ncols=pred_rad_foc*2, nrows=pred_rad_foc*2, vals = 1)
        buf_ras <- terra::rasterize(buf, r, touches = TRUE)
        mat <- as.matrix(buf_ras, wide=TRUE)
        mat_odd <- rbind(mat, NaN)
        mat_odd <- cbind(mat_odd, NaN)

        # define a second focal region for TPI
        inner_buf <- terra::buffer(vect(cbind(0,0), crs=crs(dtm)), pred_rad_foc/2)
        annulus <- terra::erase(buf, inner_buf)
        ar <- rast(xmin=-pred_rad_foc, xmax = pred_rad_foc, ymax = pred_rad_foc, ymin= -pred_rad_foc, ncols=pred_rad_foc*2, nrows=pred_rad_foc*2, vals = 1)
        annul_ras <- terra::rasterize(annulus, ar, touches = TRUE)
        Amat <- as.matrix(annul_ras, wide=TRUE)
        Amat_odd <- rbind(Amat, NaN)
        Amat_odd <- cbind(Amat_odd, NaN)


      } else {

        # for single directional  the focal region is a wedge
        pred_rad_foc <- (pred_rad/10)
        wdg <- wedge(0,0, pred_rad_foc, vi_azi, vi_fov)
        r <- rast(xmin=-pred_rad_foc, xmax = pred_rad_foc, ymax = pred_rad_foc, ymin= -pred_rad_foc, ncols=(2*pred_rad_foc), nrows=(2*pred_rad_foc), vals = 1)
        wdg_ras <- terra::rasterize(wdg, r, touches = TRUE)
        mat <- as.matrix(wdg_ras, wide=TRUE)
        mat_odd <- rbind(mat, NaN)
        mat_odd <- cbind(mat_odd, NaN)

        # define second focal region for TPI
        inner_wdg <- wedge(0,0, pred_rad_foc/2, vi_azi, vi_fov)
        annulus <- terra::erase(wdg, inner_wdg)
        ar <- rast(xmin=-pred_rad_foc, xmax = pred_rad_foc, ymax = pred_rad_foc, ymin= -pred_rad_foc, ncols=(2*pred_rad_foc), nrows=(2*pred_rad_foc))
        annul_ras <- terra::rasterize(annulus, ar, touches = TRUE)
        Amat <- as.matrix(annul_ras, wide=TRUE)
        Amat_odd <- rbind(Amat, NaN)
        Amat_odd <-cbind(Amat_odd, NaN)

      }

      # now use that matrix in your focal statistics
      # start with MEAN b/c its a built in fxn
      fslope_mean <- terra::focal(slope, mat_odd, "mean")
      fslope_derivative_mean <-  terra::focal(slope_derivative, mat_odd, "mean")
      felevation_mean <- terra::focal(elevation, mat_odd, "mean")
      fcc_mean <- terra::focal(cc, mat_odd, "mean")
      fch_mean <- terra::focal(ch, mat_odd, "mean")
      fcurvature_mean <-terra::focal(curvature, mat_odd, "mean")
      faspect_sin <- terra::focal(aspect_sin, mat_odd, "mean")
      faspect_cos <- terra::focal(aspect_cos, mat_odd, "mean")
      fslope_aspect_sin <- terra::focal(slope_aspect_sin, mat_odd, "mean")
      fslope_aspect_cos <- terra::focal(slope_aspect_cos, mat_odd, "mean")

      focal_stack <- c(fslope_mean, fslope_derivative_mean, felevation_mean, fcc_mean, fch_mean, fcurvature_mean, faspect_sin, faspect_cos, fslope_aspect_sin, fslope_aspect_cos)
      names(focal_stack) <-  c(paste0("slope_mean_", as.character(pred_rad)), paste0("slope_derivative_mean_", as.character(pred_rad)), paste0("elevation_mean_", as.character(pred_rad)), paste0("cc_mean_", as.character(pred_rad)), paste0("ch_mean_", as.character(pred_rad)), paste0("curvature_mean_", as.character(pred_rad)), paste0("aspect_sin_mean_", as.character(pred_rad)), paste0("aspect_cos_mean_", as.character(pred_rad)), paste0("slope_aspect_sin_mean_", as.character(pred_rad)), paste0("slope_aspect_cos_mean_", as.character(pred_rad)))

      # make that stack even bigggggerrrrrrr
      stack <- c(stack, focal_stack)

      # next do your standard devs -- need the custom fxn
      slope_sd <- terra::focal(x = slope, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      slope_derivative_sd <- terra::focal(x = slope_derivative, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      elevation_sd <- terra::focal(x = elevation, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      cc_sd <- terra::focal(x = cc, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      ch_sd <- terra::focal(x = ch, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      curvature_sd <- terra::focal(x = curvature, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      aspect_sin_sd <- terra::focal(x = aspect_sin, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      aspect_cos_sd <- terra::focal(x = aspect_cos, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      slope_aspect_sin_sd <- terra::focal(x = slope_aspect_sin, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })
      slope_aspect_cos_sd <- terra::focal(x = slope_aspect_cos, w = mat_odd, fun = function(x) {
        n <- length(x)
        if (n < 2) {
          stop("At least two observations are required.")
        }
        mean_x <- mean(x)
        sum_squares <- sum((x - mean_x)^2)
        sd_value <- sqrt(sum_squares / (n - 1))
        return(sd_value)
      })

      sd_stack <- c(slope_sd, slope_derivative_sd,elevation_sd, cc_sd, ch_sd, curvature_sd, aspect_sin_sd, aspect_cos_sd, slope_aspect_sin_sd, slope_aspect_cos_sd)
      names(sd_stack) <- c(paste0("slope_sd_", as.character(pred_rad)), paste0("slope_derivative_sd_", as.character(pred_rad)), paste0("elevation_sd_", as.character(pred_rad)), paste0("cc_sd_", as.character(pred_rad)), paste0("ch_sd_", as.character(pred_rad)), paste0("curvature_sd_", as.character(pred_rad)), paste0("aspect_sin_sd_", as.character(pred_rad)), paste0("aspect_cos_sd_", as.character(pred_rad)), paste0("slope_aspect_sin_sd_", as.character(pred_rad)), paste0("slope_aspect_cos_sd_", as.character(pred_rad)))

      # then EVEN BIGGER stack nowwww
      stack <- c(stack, sd_stack)

      # need TPI
      annul_elev <- terra::focal(x=elevation, w = Amat_odd, fun = "mean")
      tpi <- elevation - annul_elev

      stack <- c(stack, tpi)
      names(stack)[length(names(stack))] <- paste0("tpi_", as.character(pred_rad))

    }
    if (save == TRUE){
      writeRaster(stack, paste0(save_dir, "\\predictor_raster_stack.tif"), overwrite = TRUE)
    }

    # now extract all values to the dataframe
    xy_only <- pts %>%
      dplyr::select(c(x,y))
    vals <- terra::extract(stack, xy_only, method = "bilinear", bind=TRUE)

    df <- as.data.frame(vals, geom="XY")
    final_df <- merge(df, pts)

  }  else if (vi_type == "directional_random"){
    # start by extracting all local variables the easy way:
    xy_only <- pts %>%
      dplyr::select(c(x,y))
    vals <- terra::extract(stack, xy_only, method = "bilinear", bind=TRUE)
    df <- as.data.frame(vals, geom="XY")
    final_df <- merge(df, pts)
    pts <- final_df

    # next make sure you have an azimuth column before jumping into the loop
    if ("azimuth" %in% colnames(pts) == FALSE){
      stop("Input dataframe must have an 'azimuth' column when using vi_type 'directional_random'.")
    }

    # for this one you have to do pt by pt for FOCAL
    for (i in 1:nrow(pts)){
      # print an update every 50 points
      if(i %% 25 == 0 | i == 1){
        print(paste0("Starting point ", as.character(i), " at ", Sys.time(), " ..."))
      }

      # isolate point and azimuth
      x <- pts[i, "x"]
      y <- pts[i, "y"]
      az <- pts[i, "azimuth"]

      # make it a point
      pt <- sf::st_point(c(x,y))
      pt_sv <- terra::vect(pt)
      terra::crs(pt_sv) <- crs(dtm)

      # start with local aspect vars
      rel_asp <- aspect - az
      asp_cos <-  cos(rel_asp * pi / 180)
      asp_sin <- sin(rel_asp * pi / 180)
      sac <- slope * asp_cos
      sas <- slope * asp_sin

      # extract those values
      asp_cos_val <- terra::extract(asp_cos, pt_sv)
      asp_sin_val <- terra::extract(asp_sin, pt_sv)
      sac_val <- terra::extract(sac, pt_sv)
      sas_val <- terra::extract(sas, pt_sv)

      # add to df
      pts[i, "aspect_cos"] <- asp_cos_val[1,2]
      pts[i, "aspect_sin"] <- asp_sin_val[1,2]
      pts[i, "slope_aspect_cos"] <- sac_val[1,2]
      pts[i, "slope_aspect_sin"] <- sas_val[1,2]

      # avoid repetitive code with a list of lists
      predlist <- list(names = c("slope", "slope_derivative", "elevation",
                                 "cc", "ch", "curvature", "aspect_cos",
                                 "aspect_sin", "slope_aspect_cos", "slope_aspect_sin"),
                       rasts = list(slope, slope_derivative, elevation, cc, ch,
                                    curvature, asp_cos, asp_sin, sac, sas))
      for(pred_rad in c(20,40,60,80,160,320)){
        # draw wedge
        poly_rad <- wedge(x,y, pred_rad, az, vi_fov)
        terra::crs(poly_rad) <- terra::crs(ch)

        # loop through all predictor categories
        for (k in 1:length(predlist$names)){
          # isolate raster
          pred_ras <- predlist$rasts[[k]]

          # extract values from wedge
          ext_mean <- terra::extract(pred_ras, poly_rad, "mean", exact = TRUE)
          ext_sd <- terra::extract(pred_ras, poly_rad, "sd", touches = TRUE)

          # put value in correct column
          name <- predlist$names[k]
          colname_mean <- paste0(name, "_mean_", as.character(pred_rad))
          colname_sd <- paste0(name, "_sd_", as.character(pred_rad))
          pts[i, colname_mean] <- ext_mean[1,2]
          pts[i, colname_sd] <- ext_sd[1,2]

        }

        # deal with TPI which is special...
        # define radii for annulus
        out_rad <- pred_rad
        in_rad <- pred_rad/2

        # draw inner radius wedge
        poly_in_rad <- wedge(x,y,in_rad, az, vi_fov)

        # crop outer radius wedge to inner radius wedge to get annulus wedge
        wedge_annulus <- terra::erase(poly_rad, poly_in_rad)

        # get wedge annulus mean elevation
        wa_elev <- terra::extract(elevation, wedge_annulus, "mean", exact = TRUE)

        # get elevation value at point
        elev_val <- pts[i, "elevation"]

        # get TPI value and add to datafrme
        tpi <- elev_val - wa_elev[1,2]
        pred_col <- paste0("TPI_", as.character(pred_rad))
        pts[i, pred_col] <- tpi

      }

    }
    # rename your df to return it
    final_df <- pts

  } else {
    stop("vi_type must be one of: 'omnidir', 'directional_single', or 'directional_random'.")
  }

  print(paste0("Process completed at ", Sys.time()))
  ## at the end return and save the dataframe with all the extracted values.
  return(final_df)
}
