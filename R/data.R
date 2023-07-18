#' Digital Terrain Model 
#'
#' A 1 m digital terrain model from: 
#'
#' @format ## `dtm`
#' A SpatRaster with one band, 'dtm', representing elevation in meters, exclusive of vegetation. 
#' 
#' @source USGS LPC AZ Yavapai 2021 b 21 lidar. Tiles 6714-6712, 6812-6814, 6912-6914. Processed into 1 m raster using las2dem from lastools. <https://https://apps.nationalmap.gov/downloader/>
"dtm"

#' Digital Surface Model 
#'
#' A 1 m digital terrain model from: 
#'
#' @format ## `dsm`
#' A SpatRaster with one band, 'dsm', representing elevation in meters, inclusive of vegetation. 
#' 
#' @source USGS LPC AZ Yavapai 2021 b 21 lidar. Tiles 6714-6712, 6812-6814, 6912-6914. Processed into 1 m raster using las2dem from lastools. <https://https://apps.nationalmap.gov/downloader/>
"dsm"