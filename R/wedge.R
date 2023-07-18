#' Generate wedge-shaped polygon
#'
#' Function that draws a wedge with user-defined radius, azimuth, and theta.
#' @param x the x coordinate of the wedge origin point
#' @param y the y coordinate of the wedge origin point
#' @param rad the radius of the wedge
#' @param az the azimuth of the wedge (direction in which the wedge will be drawn) in degrees. Must be between 0 and 359.
#' @param theta the opening angle of the wedge in degrees, must be between 1 and 360
#' @return SpatVector
#' @export
#' @examples
#' # draw a east-facing wedge
#' w <- wedge(0,0, 500, 90, 90)
#'
#' # plot it
#' plot(w)



wedge <- function(x, y, rad, az, theta){
  theta.1 <- (360 + az - theta/2) %% 360
  theta.2 <- (360 + az + theta/2) %% 360
  if (theta.1 > theta.2){
    prop.1 <- (360-theta.1)/theta
    prop.2 <- 1-prop.1
    theta.is <- c(seq(theta.1,360,length.out=round(30*prop.1)),
                  seq(0,theta.2,length.out=round(30*prop.2)))
  } else {
    theta.is <- seq(theta.1,theta.2,length.out=30)
  }
  pts <- c(x,y)
  for (theta.i in theta.is){
    theta.i <- theta.i * pi / 180
    x.i <- x + rad * sin(theta.i)
    y.i <- y + rad * cos(theta.i)
    pts <- rbind(pts, c(x.i, y.i))
  }
  colnames(pts) <- c("x", "y")
  poly <- terra::vect(pts, "polygons")
  return(poly)
}
