% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_dems.R
\name{prep_dems}
\alias{prep_dems}
\title{Check and prepare DEMs for use in VisiMod workflow}
\usage{
prep_dems(in_dtm, in_dsm, out_dtm, out_dsm)
}
\arguments{
\item{in_dtm}{SpatRaster. Digital terrain model at finest resolution available.}

\item{in_dsm}{SpatRaster. Digital surface model at finest resolution available.}

\item{out_dtm}{Character. File name for saving your output DTM. Only used if NAs are found and filled.}

\item{out_dsm}{Character. File name for saving your output DSM. Only used if NAs are found and filled.}
}
\value{
If the \code{dtm} and \code{dsm} are spatially aligned, it returns a list with two SpatRasters: \code{dtm} and \code{dsm}. If interior NAs were found, the returned SpatRasters will be filled. If no interior NAs are found, the original SpatRasters will be returned.
}
\description{
This function performs three checks to determine if data are ready for the VisiMod workflow:
\itemize{
\item (1) Are the two input SpatRasters are spatially aligned (same CRS, extent, resolution, origin, and number of rows/columns)?
\item (2) Are the two input SpatRasters in appropriate projected coordinate systems, where linear units are measured in meters (e.g., UTM)?
\item (3) Are there any NA cells within the interior of either SpatRaster?
}

Interior NA values (i.e., within a polygon representing the outer dimensions of the data), will be problematic for VisiMod. If present, the function will fill these NA values using focal means from surrounding, non-NA pixels and output new DEMs.
}
\details{
\itemize{
\item This is the first function in the suggested VisiMod workflow
\item This workflow relies on two main input digital elevation models (DEMs): (1) a digital terrain model (DTM), a raster dataset where each pixel value represents the elevation of the ground surface; and (2) a digital surface model (DSM), a raster dataset where each pixel represents the elevation of the ground plus any above-ground features (e.g., trees).
\item If any of the three checks fail, the function terminates and the user is asked to remedy noted issue with the data.
\item The resulting NA-free DEMs should be used throughout the remainder of the VisiMod workflow.
}
}
\examples{
# read in your dtm and dsm
dsm <- rast("dsm.tif")
dtm <- rast("dtm.tif")

# run prep_dems()
pd <- prep_dems(dtm, dsm, "C:/temp/dtm_filled.tif", "C:/temp/dsm_filled.tif")
}
