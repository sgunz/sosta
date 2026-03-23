#' Calculate a set of shape metrics of a polygon
#'
#' @param sfPoly POLYGON of class sfc
#'
#' @return list; list of shape metrics
#'
#' @importFrom sf st_area st_length st_convex_hull st_boundary st_geometry_type
#' @export
#'
#' @examples
#' matrix_R <- matrix(c(
#'     0, 0, 0, 0, 0, 0, 0, 0, 0,
#'     0, 1, 1, 1, 1, 1, 0, 0, 0,
#'     0, 1, 1, 0, 0, 1, 1, 0, 0,
#'     0, 1, 1, 0, 0, 1, 1, 0, 0,
#'     0, 1, 1, 1, 1, 1, 0, 0, 0,
#'     0, 1, 1, 0, 1, 1, 0, 0, 0,
#'     0, 1, 1, 0, 0, 1, 1, 0, 0,
#'     0, 1, 1, 0, 0, 1, 1, 0, 0,
#'     0, 0, 0, 0, 0, 0, 0, 0, 0
#' ), nrow = 9, byrow = TRUE)
#' polyR <- binaryImageToSF(matrix_R, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#' shapeMetrics(polyR)
shapeMetrics <- function(sfPoly) {
    # Input checks
    stopifnot("'sfPoly' must be a valid sfc object" = inherits(sfPoly, "sfc"))
    stopifnot("'sfPoly' must be of type POLYGON" = st_geometry_type(sfPoly) == "POLYGON")
    stopifnot("'sfPoly' must be a single POLYGON,
    use `sosta::totalShapeMetrics` for a set of POLYGONs" = length(st_geometry_type(sfPoly)) == 1)
    # Area
    shapeArea <- st_area(sfPoly)
    # Perimeter
    shapePerimeter <- st_length(st_boundary(sfPoly))
    # Convex Hull
    shapeConvexHull <- st_convex_hull(sfPoly)
    # Perimeter of convex hull
    shapeConvexPerimeter <- st_length(st_boundary(shapeConvexHull))
    # Feature axes
    FeatureAxes <- stFeatureAxes(sfPoly)
    # Compactness: 0 and 1 (circle)
    shapeCompactness <- (4 * pi * shapeArea) / (shapePerimeter)^2
    # Eccentricity: between 0 and 1
    shapeEccentricity <- FeatureAxes$minorAxisLength / FeatureAxes$majorAxisLength
    # Circularity / roundness: 0 and 1 for round object
    shapeCircularity <- (4 * pi * shapeArea) / shapeConvexPerimeter^2
    #  Curl
    shapeCurl <- stCalculateShapeCurl(sfPoly)
    # Solidity: measures the density of an object
    shapeSolidity <- shapeArea / st_area(shapeConvexHull)

    # Return all metrics in a list
    return(list(
        Area = shapeArea,
        Compactness = shapeCompactness,
        Eccentricity = shapeEccentricity,
        Circularity = shapeCircularity,
        Solidity = shapeSolidity,
        Curl = shapeCurl$Curl,
        fibreLength = shapeCurl$fibreLength,
        fibreWidth = shapeCurl$fibreWidth
    ))
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
#' @importFrom sf st_cast st_geometry st_sfc st_geometry_type
#' @export
#'
#' @examples
#' data(sostaSPE)
#' struct <- reconstructShapeDensityImage(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     imageId = "image1", markSelect = "A", dim = 500
#' )
#' totalShapeMetrics(struct)
totalShapeMetrics <- function(sfInput) {
    # Input checks
    stopifnot("'sfInput' must be a valid sf object" = inherits(sfInput, "sf"))
    stopifnot("'sfInput' must be of type MULTIPOLYGON" = all(st_geometry_type(sfInput) == "POLYGON"))

    # cast into different objects, i.e the substructures
    cast_sf <- st_cast(st_geometry(sfInput), "POLYGON")
    # calculate tissue metrics on all substructures
    if (length(cast_sf) > 1) {
        shapeStruct <- vapply(
            st_geometry(cast_sf),
            function(x) unlist(shapeMetrics(st_sfc(x))),
            numeric(8) # length of metrics
        )
    } else {
        shapeStruct <- t(data.frame(shapeMetrics(cast_sf)))
    }
    # matrix of metrics of all substructures
    shapeMat <- matrix(as.numeric(shapeStruct), nrow = dim(shapeStruct)[1])
    rownames(shapeMat) <- rownames(shapeStruct)

    if (!is.null(sfInput[["structID"]])) {
        colNames <- sfInput[["structID"]]
    } else {
        colNames <- paste0(
            deparse(substitute(sfInput)),
            seq_len(dim(shapeMat)[2])
        )
    }

    colnames(shapeMat) <- colNames
    return(shapeMat)
}


#' Calculate mean shape metrics of a set of polygons
#'
#' @param totalShapeMetricMatrix matrix of shape metrics
#'
#' @return matrix; matrix of mean shape metrics
#' @export
#' @examples
#' data(sostaSPE)
#' struct <- reconstructShapeDensityImage(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     imageId = "image1", markSelect = "A", dim = 500
#' )
#' shapeMetrics <- totalShapeMetrics(struct)
#' meanShapeMetrics(shapeMetrics)
meanShapeMetrics <- function(totalShapeMetricMatrix) {
    # Check Input
    stopifnot("'totalShapeMetricMatrix' must be a matrix" = is.matrix(totalShapeMetricMatrix))

    meanShapeMat <- rowMeans(totalShapeMetricMatrix)
    meanShapeMat["numberStructures"] <- dim(totalShapeMetricMatrix)[2]
    meanShapeMat["totalArea"] <- sum(totalShapeMetricMatrix["Area", ])

    meanShapeMat <- as.matrix(meanShapeMat)
    meanShapeMat[is.na(meanShapeMat)] <- 0
    return(meanShapeMat)
}
