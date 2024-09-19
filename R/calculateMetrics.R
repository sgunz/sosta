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
#' poly_R <- binaryImageToSF(matrix_R, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#' shapeMetrics(poly_R)
shapeMetrics <- function(sfPoly) {
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
    shapeCompactness <- (shapePerimeter)^2 / (4 * pi * shapeArea)
    # Eccentricity: between 0 and 1
    shapeEccentricity <- FeatureAxes$minorAxisLength / FeatureAxes$majorAxisLength
    # Circularity / roundness: 0 and 1 for round object
    shapeCircularity <- (4 * pi * shapeArea) / shapeConvexPerimeter^2
    #  Curl
    shapeCurl <- st_calculateShapeCurl(sfPoly)
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
#' @import sf
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' islet_poly <- reconstructShapeDensityImage(spe,
#'     marks = "cell_category",
#'     image_col = "image_name", image_id = "E04", mark_select = "islet", dim = 500
#' )
#' totalShapeMetrics(islet_poly)
totalShapeMetrics <- function(sfInput) {
    # cast into different objects, i.e the substructures
    cast_sf <- st_cast(st_geometry(sfInput), "POLYGON")
    # calculate tissue metrics on all substructures
    if (length(cast_sf) > 1) {
        shapeStruct <- vapply(
            st_geometry(cast_sf),
            function(x) unlist(shapeMetrics(st_sfc(x))),
            numeric(8) #length of metrics
        )
    } else {
        shapeStruct <- t(data.frame(shapeMetrics(cast_sf)))
    }
    # matrix of metrics of all substructures
    shapeMat <- matrix(as.numeric(shapeStruct), nrow = dim(shapeStruct)[1])
    rownames(shapeMat) <- rownames(shapeStruct)
    colnames(shapeMat) <- paste0(
        deparse(substitute(sfInput)),
        seq_len(dim(shapeMat)[2])
    )
    return(shapeMat)
}


#' Calculate mean shape metrics of a set of polygons
#'
#' @param totalShapeMetricMatrix matrix of shape metrics
#'
#' @return matrix; matrix of mean shape metrics
#' @export
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' islet_poly <- reconstructShapeDensityImage(spe,
#'     marks = "cell_category",
#'     image_col = "image_name", image_id = "E04", mark_select = "islet", dim = 500
#' )
#' shape_metrics <- totalShapeMetrics(islet_poly)
#' meanShapeMetrics(shape_metrics)
meanShapeMetrics <- function(totalShapeMetricMatrix) {
    meanShapeMat <- rowMeans(totalShapeMetricMatrix)
    meanShapeMat["numberStructures"] <- dim(totalShapeMetricMatrix)[2]
    meanShapeMat["totalArea"] <- sum(totalShapeMetricMatrix["Area", ])

    meanShapeMat <- as.matrix(meanShapeMat)
    meanShapeMat[is.na(meanShapeMat)] <- 0
    return(meanShapeMat)
}
