#' Check and prepare DEMs for use in VisiMod workflow
#'
#' @description
#' The first suggested function in the VisiMod workflow. This workflow relies on two main input digital elevation models (DEMs): (1) a digital terrain model (DTM), a raster dataset where each pixel value represents the elevation of the ground surface; and (2) a digital surface model (DSM), a raster dataset where each pixel represents the elevation of the ground plus any above-ground features (e.g., trees). It is recommended that these DEMs are generated at a high spatial resolution (e.g., 1m) from lidar data. This function ensures that both of these inputs are suitable for use in the VisiMod workflow. It begins by checking whether the two datasets are spatially aligned (same CRS, extent, resolution, origin, and number of rows/columns). If not, the function terminates and the user is asked to ensure proper alignment before proceeding with other VisiMod functions. If so, a secondary check is performed to determine if they are in appropriate projected coordinate systems, where linear units are measured in meters (e.g., UTM). If not, the function terminates and the user is asked to reproject the data before proceeding with other VisiMod functions. If so, a tertiary check is performed to determine if there are any NA cells within the interior of either dataset. Interior NA values (i.e., within a polygon representing the outer dimensions of the data), will be problematic for VisiMod. If present, the function will fill these NA values using focal means from surrounding, non-NA pixels and output new DEMs. The resulting NA-free DEMs should be used throughout the remainder of the VisiMod workflow.
#' 
#' @param in_dtm SpatRaster. Digital terrain model at finest resolution available.
#' @param in_dsm SpatRaster. Digital surface model at finest resolution available.
#' @param out_dtm Character. File name for saving your output DTM. Only used if NAs are found and filled.
#' @param out_dsm Character. File name for saving your output DSM. Only used if NAs are found and filled.
#' @return If `dtm` and `dsm` are spatially aligned, a list with two SpatRasters: `dtm` and `dsm`. If interior NAs were found, the returned SpatRasters will be filled. If not, the original SpatRasters will be returned.
#'
#' @export
#' @examples
#' # read in your dtm and dsm
#' dsm <- rast("dsm.tif")
#' dtm <- rast("dtm.tif")
#'
#' # run prep_dems()
#' pd <- prep_dems(dtm, dsm, "C:/temp/dtm_filled.tif", "C:/temp/dsm_filled.tif")

prep_dems <- function(in_dtm, in_dsm, out_dtm, out_dsm){
  
  #-------------------alignment
  
  # print message
  message(paste0(Sys.time(), " prep_dems() has begun"))
  message(paste0(Sys.time(), "   Checking for spatial alignment between dtm and dsm..."))
  
  # compare crs
  crs_dtm <- terra::crs(in_dtm)
  crs_dsm <- terra::crs(in_dsm)
  crs_comp <- crs_dtm == crs_dsm
  if (crs_comp) message(paste0(Sys.time(), "     crs: match"))
  if (!crs_comp) message(paste0(Sys.time(), "     crs: do not match"))
  
  # compare ext
  ext_dtm <- terra::ext(in_dtm)
  ext_dsm <- terra::ext(in_dsm)
  ext_comp <- ext_dtm == ext_dsm
  if (ext_comp) message(paste0(Sys.time(), "     ext: match"))
  if (!ext_comp) message(paste0(Sys.time(), "     ext: do not match"))
  
  # compare nrow
  nrow_dtm <- terra::nrow(in_dtm)
  nrow_dsm <- terra::nrow(in_dsm)
  nrow_comp <- nrow_dtm == nrow_dsm
  if (nrow_comp) message(paste0(Sys.time(), "     nrow: match"))
  if (!nrow_comp) message(paste0(Sys.time(), "     nrow: do not match"))
  
  # compare ncol
  ncol_dtm <- terra::ncol(in_dtm)
  ncol_dsm <- terra::ncol(in_dsm)
  ncol_comp <- ncol_dtm == ncol_dsm
  if (ncol_comp) message(paste0(Sys.time(), "     ncol: match"))
  if (!ncol_comp) message(paste0(Sys.time(), "     ncol: do not match"))
  
  # compare res
  res_dtm <- terra::res(in_dtm)
  res_dsm <- terra::res(in_dsm)
  res_comp <- all(res_dtm == res_dsm)
  if (res_comp) message(paste0(Sys.time(), "     res: match"))
  if (!res_comp) message(paste0(Sys.time(), "     res: do not match"))
  
  # compare origin
  origin_dtm <- terra::origin(in_dtm)
  origin_dsm <- terra::origin(in_dsm)
  origin_comp <- all(origin_dtm == origin_dsm)
  if (origin_comp) message(paste0(Sys.time(), "     origin: match"))
  if (!origin_comp) message(paste0(Sys.time(), "     origin: do not match"))
  
  # compare all
  if (!(crs_comp & ext_comp & nrow_comp & ncol_comp & res_comp & origin_comp)){
    stop("\nYou must ensure spatial alignment of dtm and dsm before proceeding.\n")
  } else {
    message(paste0(Sys.time(), "   dtm and dsm are spatially aligned"))
  }
  
  #-------------------linear units
  
  # print message
  message(paste0(Sys.time(), "   Checking to make sure dtm/dsm coordinate system has linear units in meters..."))
  
  # check linear units
  units <- terra::linearUnits(in_dtm)
  if (!units == 1){
    stop("\nYou must ensure that your input dtm and dsm are in a projected coordinate system with linear units in meters (e.g., UTM).\n")  
  } else {
    message(paste0(Sys.time(), "   dtm and dsm have suitable coordinate systems"))
  }
  
  #-------------------interior NAs
  
  # print message
  message(paste0(Sys.time(), "   Checking for interior NA values..."))
  
  # create study area boundary
  sab <- terra::ifel(!is.na(in_dtm), 1, NA) |>
    terra::as.polygons() |>
    terra::fillHoles()
  
  # create output flags
  dtm_fill_flag <- F
  dsm_fill_flag <- F
  
  # fill nas -- dtm
  na_count <- terra::freq(in_dtm, value = NA, zones = sab)$count
  while(na_count != 0){
    dtm_fill_flag <- T
    message(paste0(Sys.time(), "    dtm has ", na_count, " NAs. Removing..."))
    in_dtm <- terra::focal(in_dtm, 9, mean, na.policy = "only", na.rm = T)
    na_count <- terra::freq(in_dtm, value = NA, zones = sab)$count
  }
  message(paste0(Sys.time(), "    dtm has no interior NA values"))
  if (dtm_fill_flag) in_dtm <- terra::crop(in_dtm, sab, mask = T)
  
  # fill nas -- dsm
  na_count <- terra::freq(in_dsm, value = NA, zones = sab)$count
  while(na_count != 0){
    dsm_fill_flag <- T
    message(paste0(Sys.time(), "    dsm has ", na_count, " NAs. Removing..."))
    in_dsm <- terra::focal(in_dsm, 9, mean, na.policy = "only", na.rm = T)
    na_count <- terra::freq(in_dsm, value = NA, zones = sab)$count
  }
  message(paste0(Sys.time(), "    dsm has no interior NA values"))
  if (dsm_fill_flag) in_dsm <- terra::crop(in_dsm, sab, mask = T)
  
  # write to file
  if (dtm_fill_flag) terra::writeRaster(in_dtm, out_dtm, overwrite = T)
  if (dsm_fill_flag) terra::writeRaster(in_dsm, out_dsm, overwrite = T)
  
  # read back in and return as SpatRasters
  if (dtm_fill_flag) in_dtm <- terra::rast(out_dtm)
  if (dsm_fill_flag) in_dsm <- terra::rast(out_dsm)
  message(paste0(Sys.time(), " prep_dems() is complete"))
  return(list(dtm = in_dtm, dsm = in_dsm))
  
}