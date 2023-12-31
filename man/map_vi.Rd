% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_vi.R
\name{map_vi}
\alias{map_vi}
\title{Map VI across a landscape}
\usage{
map_vi(
  mod,
  predictors,
  cores = floor(parallel::detectCores()/2),
  save = T,
  out_file = file.path(getwd(), "vi.tif")
)
}
\arguments{
\item{mod}{ranger. A random forest model generated using \code{\link[ranger:ranger]{ranger::ranger()}} aimed at predicting VI.}

\item{predictors}{SpatRaster. A multiband SpatRaster of all predictors used for VI modeling.}

\item{cores}{Numeric. Defines the number of cores you would like to use for parallel processing. Defaults to half of the cores on your machine.}

\item{save}{Boolean. Defines whether or not you would like to save the resulting VI prediction map to file. Default is TRUE.}

\item{out_file}{Character. The full file path (including the ".tif" extension) of your output VI map.}
}
\value{
A SpatRaster containing pixel-level VI predictions.
}
\description{
This function takes a random forest model (generated using \code{\link[=mod_vi]{mod_vi()}}), representing the predictive relationship between visibility index (VI) (calculated at a series of point locations using \code{\link[=calc_vi]{calc_vi()}}), and a suite of visibility predictor layers representing terrain and vegetation structure (derived using \code{\link[=gen_preds]{gen_preds()}}), and applies this model to the prediction of VI across an entire study area of interest.
}
\details{
Successful execution of this function will result in a spatially exhaustive ("wall-to-wall") map, where each pixel represents an estimate of VI based on a spatial scope of interest (viewing distance, and direction and field of view, if mapping directional VI).
\itemize{
\item This is the sixth and final suggested function in the VisiMod workflow, following (1) \code{\link[=prep_dems]{prep_dems()}} (2) \code{\link[=gen_pts]{gen_pts()}}, (3) \code{\link[=calc_vi]{calc_vi()}}, (4) \code{\link[=gen_preds]{gen_preds()}}, and (5) \code{\link[=mod_vi]{mod_vi()}}.
\item The spatial resolution, coordinate system, origin, and extent of the resulting map will be based entirely on those of your multiband predictors raster.
\item It is recommended to be somewhat conservative with your definition of \code{cores}. Using a large proportion of your computer's cores can cause a memory overload, as a large amount of information has to be stored in memory for each core.
\item Note that, by default, if \code{save == TRUE} and you supply an existing filename to \code{out_file}, it will be overwritten.
}
}
\examples{

# get your dtm and dsm
dsm <- rast("dsm.tif")
dtm <- rast("dtm.tif")

# check dtm and dsm
pd <- prep_dems(dtm, dsm, "C:/temp/dtm_filled.tif", "C:/temp/dsm_filled.tif")

# get your points
my_points <- gen_pts(dtm, dsm, 100, 1000)

# calculate vi
my_points <- calc_vi(dtm, dsm, my_points, "directional_single", c(500, 1000), 90, 90, 4L, 5L)

# generate predictors
preds <- gen_preds(dtm, dsm, my_points, "directional_single", 90, 90, T, "C:/temp")

# model vi
mod <- mod_vi(preds, 500, T, T, 5L)

# map vi
rf <- mod$ranger_mod
predictors <- preds$pred_rast
mymap <- map_vi(rf, predictors, 5L, T, "C:/temp/vi_map.tif")
plot(mymap)
}
