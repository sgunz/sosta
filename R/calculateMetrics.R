#' Calculate a set of shape metrics of a polygon
#'
#' @param sfPoly POLYGON of class sf
#'
#' @return list; list of shape metrics
#'
#' @import sf
#' @export
#'
#' @examples

shapeMetrics <- function(sfPoly){
  # Area
  shapeArea <- st_area(sfPoly)
  # Perimeter
  shapePerimeter <- st_length(st_boundary(sfPoly))
  # Convex Hull
  shapeConvexHull <- st_convex_hull(sfPoly)
  # Perimeter of convex hull
  shapeConvexPerimeter <- st_length(st_boundary(shapeConvexHull))
  # Feature axes
  FeatureAxes <- st_feature_axes(sfPoly)
  # Compactness: 0 and 1 (circle)
  shapeCompactness <- (4 * pi * shapeArea) /
    (shapePerimeter)^2
  # Eccentricity: between 0 and 1
  shapeEccentricity <- FeatureAxes$minorAxis /
    FeatureAxes$majorAxis
  # Circularity / roundness: 0 and 1 for round object
  shapeCircularity <- (4 * pi * shapeArea) / shapeConvexPerimeter^2
  #  Curl
  shapeCurl <- st_calculateShapeCurl(sfPoly)
  # Solidity: measures the density of an object
  shapeSolidity <- shapeArea / st_area(shapeConvexHull)
  # Max. inscribed circle
  areaInscribedCricle <- max(st_area(st_inscribed_circle((sfPoly))))
  # Curvature
  curvature <- st_calculateCurvature(sfPoly)
  # Area of holes
  AreaHoles <- abs(st_area(st_multipolygon(
    lapply(st_union(sfPoly), function(x) x[1]))) -
      st_area(sfPoly))
  # Return all metrics in a list
  return(list(
    Area = shapeArea,
    Perimeter = shapePerimeter,
    ConvexPerimeter = shapeConvexPerimeter,
    MajorAxis = FeatureAxes$majorAxisLength,
    MinorAxis = FeatureAxes$minorAxisLength,
    Compactness = shapeCompactness,
    Eccentricity = shapeEccentricity,
    Circularity = shapeCircularity,
    Solidity = shapeSolidity,
    areaMaxInscribedCircle = areaInscribedCricle,
    Curl = shapeCurl,
    AreaHoles = AreaHoles,
    meanAbsCurvature = curvature$meanAbsCurv,
    totalAbsCurvature = curvature$totalAbsCurv,
    minCurvature = curvature$minCurv,
    maxCurvature = curvature$maxCurv))
}



#' Calculate a set of shape metrics of a set of polygons
#'
#' @details
#' Calculate a set of shape metrics of a set of polygons.
#' The function calculates all metrics that are implemented in the function
#' `shapeMetrics()`
#'
#' @param sfInput `MULTIPOLYGON` of class sf
#'
#' @return matrix; matrix of shape metrics
#' @import sf
#' @export
#'
#' @examples
totalShapeMetrics <- function(sfInput){
  # cast into different objects, i.e the substructures
  cast_sf <- suppressWarnings(st_cast(sfInput,  "POLYGON"))
  # calculate tissue metrics on all substructures
  if (length(cast_sf) > 1) {
    shapeStruct <- sapply(st_geometry(cast_sf),
                          function(x) shapeMetrics(st_sfc(x)))
  }
  else {shapeStruct <- t(data.frame(shapeMetrics(cast_sf)))}
  # matrix of metrics of all substructures
  shapeMat <- matrix(as.numeric(shapeStruct), nrow = dim(shapeStruct)[1])
  rownames(shapeMat) <- rownames(shapeStruct)
  colnames(shapeMat) <- paste0(deparse(substitute(sfInput)),
                               rep(c(1:dim(shapeMat)[2]), each = 1))
  return(shapeMat)
}


#' Calculate mean shape metrics of a set of polygons
#'
#' @param totalShapeMetricMatrix matrix of shape metrics
#'
#' @return matrix; matrix of mean shape metrics
#' @import stats
#' @export
#'
#' @examples
meanShapeMetrics <- function(totalShapeMetricMatrix){

  meanShapeMat <- rowMeans(totalShapeMetricMatrix)
  meanShapeMat['numberStructures'] <- dim(totalShapeMetricMatrix)[2]
  meanShapeMat['totalArea'] <- sum(totalShapeMetricMatrix['Area',])
  meanShapeMat['sdMeanCurvature'] <- sd(totalShapeMetricMatrix['meanAbsCurvature',])
  meanShapeMat['sdEccentricity'] <- sd(totalShapeMetricMatrix['Eccentricity',])

  meanShapeMat <- as.matrix(meanShapeMat)
  meanShapeMat[is.na(meanShapeMat)] <- 0
  return(meanShapeMat)
}
