#' Find the center of a set of coordinates
#'
#' @param sfe SpatialFeatureExperiment
#' @param m vector; vector of TRUE/FALSE indicating the indeces of the cells of interest
#' @param mark character; name of a column in `colData` that corresponds to the marks of interest
#' @param targetMark character; mark of interest, must be present in mark
#'
#' @return
#' @export
#' @import SpatialFeatureExperiment
#'
#' @examples
findCenterCoord <- function(sfe, m, mark, targetMark) {

  # select centroids that fit the target mark
  sel <- centroids(sfe[,m])[colData(sfe[,m])[[mark]] == targetMark,]
  # TODO: add function to calculate based on medians
  centers <- st_coordinates(sel) |> colMeans()

  return(c(centers[1], centers[2]))
}


#' Function to generate gradient / line based on a set of points
#'
#' @param coords matrix; matrix of coordinates
#' @param geom 	object of class sf, sfc or sfg
#'
#' @return
#' @export
#'
#' @examples
calculateLine <- function(coords, geom) {

  # fit a linear model
  mod <- lm(y ~ x, coords)

  # predict new data
  new_data <- data.frame(x = c(st_bbox(geom)$xmin, st_bbox(geom)$xmax))
  ynew <- predict(mod, new_data)

  # create a new line
  new_data$y <- ynew
  newline <- st_linestring(as.matrix(new_data))

  return(newline)
}
