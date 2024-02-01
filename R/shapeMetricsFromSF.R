#' Calculate the length of feature axes of an sf polygon
#'
#' @param sfPoly `POLYGON ` of class `sf`
#'
#' @return list;  with the major and minor axis lengths
#' @export
#'
#' @examples

st_feature_axes <- function(sfPoly) {
  minRect <- sf::st_minimum_rotated_rectangle(sfPoly)

  coords <- sf::st_coordinates(minRect)[, c("X", "Y")]

  # Calculate the lengths of the rectangle sides
  side1Length <- sqrt((coords[1, "X"] - coords[2, "X"])^2 + (coords[1, "Y"] - coords[2, "Y"])^2)
  side2Length <- sqrt((coords[2, "X"] - coords[3, "X"])^2 + (coords[2, "Y"] - coords[3, "Y"])^2)

  # Determine which side is the major axis
  majorAxisLength <- max(side1Length, side2Length)
  minorAxisLength <- min(side1Length, side2Length)

  # Return the mean and sum of the curvature values
  return(list(
    majorAxisLength = majorAxisLength,
    minorAxisLength = minorAxisLength))
}

#' Title
#'
#' @param sfPoly `POLYGON ` of class `sf`
#' @param smoothness list; curvature measures
#'
#' @return
#' @export
#'
#' @references https://stackoverflow.com/questions/62250151/calculate-curvature-of-a-closed-object-in-r
#'
#' @examples
st_calculateCurvature <- function(sfPoly, smoothness = 5) {
  # Smooth the data using concave hull and ksmooth method
  smooth_poly <- smoothr::smooth(sf::st_boundary(sfPoly),
                                 method = 'ksmooth',
                                 smoothness = smoothness)

  smooth_poly <- sfPoly

  # Convert the smoothed data to a data frame
  smooth_df <- as.data.frame(sf::st_coordinates(smooth_poly))

  # Calculate the change in x and y per unit length
  dx <- diff(c(smooth_df$X, smooth_df$X[1]))
  dy <- diff(c(smooth_df$Y, smooth_df$Y[1]))
  ds <- sqrt(dx^2 + dy^2)
  ddx <- dx/ds
  ddy <- dy/ds
  ds2 <- (ds + c(ds[-1], ds[1]))/2
  smooth_df$Cx <- diff(c(ddx, ddx[1]))/ds2
  smooth_df$Cy <- diff(c(ddy, ddy[1]))/ds2

  # Calculate the change in curvature (K) per unit length
  smooth_df$K <- (ddy * smooth_df$Cx - ddx * smooth_df$Cy) /
    ((ddx^2 + ddy^2)^(3/2))
  # Replace NA's
  smooth_df$K[length(smooth_df$K)] <- smooth_df$K[1]
  smooth_df$K[length(smooth_df$K) - 1] <- smooth_df$K[length(smooth_df$K) - 2]

  smooth_df$K[smooth_df$K > 10] <- 5
  smooth_df$K[smooth_df$K < -10] <- -5

  # Return the mean and sum of the curvature values
  return(list(
    meanAbsCurv = mean(abs(smooth_df$K[!is.na(smooth_df$K)])),
    totalAbsCurv = sum(abs(smooth_df$K[!is.na(smooth_df$K)])),
    minCurv = min(smooth_df$K[!is.na(smooth_df$K)]),
    maxCurv = max(smooth_df$K[!is.na(smooth_df$K)]))
  )
}


#' Calculate curls of a polygon
#'
#' @param sfPoly `POLYGON ` of class `sf`
#'
#' @return numeric; the curl of the polygon
#' @export
#'
#' @examples
st_calculateShapeCurl <- function(sfPoly){
  # Major axis length
  length <- st_feature_axes(sfPoly)$majorAxisLength
  # Calculate perimeter
  perimeter <- sf::st_length(sf::st_boundary(sfPoly))
  # Calculate radicand
  radicand <- perimeter^2 - 16 * sf::st_area(sfPoly)
  # Check if sqrt of radicand is real
  radicand <- ifelse(radicand > 0, radicand, 0)
  # Calculate fibre length
  fibreLength <- (perimeter - sqrt(radicand)) / 4
  # Return curl
  return(length/fibreLength)
}
