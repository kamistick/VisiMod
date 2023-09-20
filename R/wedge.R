#' Generate wedge-shaped polygon
#'
#' Function that draws a wedge (angular slice of a circle from a central origin point) with user-defined radius, azimuth (central viewing direction), and theta (angular field of view). This is a utility function within the VisiMod library that is used in several functions to aid in directional visibility estimation.
#' @param x Numeric. The x coordinate of the wedge origin point.
#' @param y Numeric. The y coordinate of the wedge origin point.
#' @param rad Numeric. The radius of the wedge.
#' @param az Numeric. The azimuth of the wedge (direction in which the wedge will be drawn) in degrees. Must be between 0 and 359.
#' @param theta Numeric. The angular field of view of the wedge in degrees, must be between 1 and 360.
#' @return SpatVector representation of a wedge.
#' @export
#' @examples
#' # draw a east-facing wedge
#' w <- wedge(0,0, 500, 90, 90)
#'
#' # plot it
#' plot(w)

wedge <- function(x, y, rad, az, theta){
  
  # define the angles of the wedges two radial lines
  theta.1 <- (360 + az - theta/2) %% 360
  theta.2 <- (360 + az + theta/2) %% 360
  
  # define the angles of a series of radii between those two lines
  if (theta.1 > theta.2){
    prop.1 <- (360-theta.1)/theta
    prop.2 <- 1-prop.1
    theta.is <- c(seq(theta.1,360,length.out=round(30*prop.1)),
                  seq(0,theta.2,length.out=round(30*prop.2)))
  } else {
    theta.is <- seq(theta.1,theta.2,length.out=30)
  }
  
  # create the points that define the wedge's geometry
  pts <- c(x,y)
  for (theta.i in theta.is){
    theta.i <- theta.i * pi / 180
    x.i <- x + rad * sin(theta.i)
    y.i <- y + rad * cos(theta.i)
    pts <- rbind(pts, c(x.i, y.i))
  }
  
  # convert to SpatVector
  colnames(pts) <- c("x", "y")
  poly <- terra::vect(pts, "polygons")
  
  # return the result
  return(poly)
  
}
