#' Calculate the length of feature axes of a single sf polygon
#'
#' @param sfPoly `POLYGON ` of class `sf`
#' @importFrom sf st_minimum_rotated_rectangle st_coordinates
#'
#' @return list; list containing the major and minor axis lengths
#' @export
#'
#' @examples
#' matrixR <- matrix(c(
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
#' polyR <- binaryImageToSF(matrixR, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#' stFeatureAxes(polyR)
stFeatureAxes <- function(sfPoly) {
    # Input checks
    stopifnot("'sfPoly' must be a valid sfc object" = inherits(sfPoly, "sfc"))
    stopifnot("'sfPoly' must be of type POLYGON" = st_geometry_type(sfPoly) == "POLYGON")

    minRect <- st_minimum_rotated_rectangle(sfPoly)
    coords <- st_coordinates(minRect)[, c("X", "Y")]

    # Calculate the lengths of the rectangle sides
    side1Length <- sqrt((coords[1, "X"] - coords[2, "X"])^2 + (coords[1, "Y"] - coords[2, "Y"])^2)
    side2Length <- sqrt((coords[2, "X"] - coords[3, "X"])^2 + (coords[2, "Y"] - coords[3, "Y"])^2)

    # Return the mean and sum of the curvature values
    return(list(
        majorAxisLength = max(side1Length, side2Length),
        minorAxisLength = min(side1Length, side2Length)
    ))
}

#' Calculate curvature of sf object
#'
#' @param sfPoly `POLYGON ` of class `sf`
#' @param smoothness list; curvature measures
#'
#' @return list; list of curvatures values
#' @importFrom sf st_boundary st_coordinates
#' @importFrom smoothr smooth
#' @export
#'
#' @references https://stackoverflow.com/questions/62250151/calculate-curvature-of-a-closed-object-in-r
#' @examples
#' matrixR <- matrix(c(
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
#' polyR <- binaryImageToSF(matrixR, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#' stCalculateCurvature(polyR)
stCalculateCurvature <- function(sfPoly, smoothness = 5) {
    # Input checks
    stopifnot("'sfPoly' must be a valid sfc object" = inherits(sfPoly, "sfc"))
    stopifnot("'sfPoly' must be of type POLYGON" = st_geometry_type(sfPoly) == "POLYGON")

    # Smooth the data using concave hull and ksmooth method
    smooth_poly <- smooth(st_boundary(sfPoly),
        method = "ksmooth",
        smoothness = smoothness
    )

    smooth_poly <- sfPoly
    # Convert the smoothed data to a data frame
    smooth_df <- as.data.frame(st_coordinates(smooth_poly))

    # Calculate the change in x and y per unit length
    dx <- diff(c(smooth_df$X, smooth_df$X[1]))
    dy <- diff(c(smooth_df$Y, smooth_df$Y[1]))
    ds <- sqrt(dx^2 + dy^2)
    ddx <- dx / ds
    ddy <- dy / ds
    ds2 <- (ds + c(ds[-1], ds[1])) / 2
    smooth_df$Cx <- diff(c(ddx, ddx[1])) / ds2
    smooth_df$Cy <- diff(c(ddy, ddy[1])) / ds2

    # Calculate the change in curvature (K) per unit length
    smooth_df$K <- (ddy * smooth_df$Cx - ddx * smooth_df$Cy) /
        ((ddx^2 + ddy^2)^(3 / 2))
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
        maxCurv = max(smooth_df$K[!is.na(smooth_df$K)])
    ))
}


#' Calculate curl of a polygon
#'
#' @param sfPoly `POLYGON ` of class `sf`
#'
#' @return numeric; the curl of the polygon
#' @importFrom sf st_area st_length st_boundary
#' @export
#'
#' @examples
#' matrixR <- matrix(c(
#'     1, 1, 1, 1, 1, 0,
#'     1, 1, 0, 0, 1, 1,
#'     1, 1, 0, 0, 1, 1,
#'     1, 1, 1, 1, 1, 0,
#'     1, 1, 0, 1, 1, 0,
#'     1, 1, 0, 0, 1, 1,
#'     1, 1, 0, 0, 1, 1
#' ), nrow = 7, byrow = TRUE)
#' polyR <- binaryImageToSF(matrixR, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
#' stCalculateShapeCurl(polyR)
stCalculateShapeCurl <- function(sfPoly) {
    # Input checks
    stopifnot("'sfPoly' must be a valid sfc object" = inherits(sfPoly, "sfc"))
    stopifnot("'sfPoly' must be of type POLYGON" = st_geometry_type(sfPoly) == "POLYGON")
    # Major axis length
    length <- stFeatureAxes(sfPoly)$majorAxisLength
    # Calculate perimeter
    perimeter <- st_length(st_boundary(sfPoly))
    # Calculate radicand
    radicand <- perimeter^2 - 16 * st_area(sfPoly)
    # Check if sqrt of radicand is real
    radicand <- ifelse(radicand > 0, radicand, 0)
    # Calculate fibre length
    fibreLength <- 4 * st_area(sfPoly) / (perimeter - sqrt(radicand))
    # Return metrics
    return(list(
        Curl = (1 - (length / fibreLength)),
        fibreLength = fibreLength,
        fibreWidth = st_area(sfPoly) / fibreLength
    ))
}
