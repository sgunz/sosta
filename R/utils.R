#' Converts a binary matrix to an sf polygon
#'
#' @param binaryMatrix matrix; binary matrix
#' @param xmin integer; minimum x coordinate of the window
#' @param xmax integer; maximum x coordinate of the window
#' @param ymin integer; minimum y coordinate of the window
#' @param ymax integer; maximum y coordinate of the window
#'
#' @return sf object
#' @import sf
#' @importFrom terra rast as.polygons ext
#' @export
#'
#' @examples
binaryImageToSF <- function(binaryMatrix, xmin, xmax,
                            ymin, ymax) {
    # turn 90 degrees anti clockwise for correspondance with spatstat
    binaryMatrix <- apply(t(binaryMatrix), 2, rev)
    # get raster
    r <- rast(binaryMatrix)
    # rescale to correct windwow
    ext(r) <- c(xmin, xmax,
                ymin, ymax)
    # convert to polygons
    poly <- as.polygons(r)
    # polygons is a SpatVector. Convert it to an sf object
    polygonsSF <- st_as_sf(poly)
    # Merge polygons to a single multipolygon
    return(st_union(polygonsSF[polygonsSF$lyr.1 == 1, ]))
}



#' Function to extract x y coordinates from binary image
#'
#' @param inputMatrix a binary matrix
#'
#' @return matrix; matrix with x,y coordinates of the cell of the input matrix
#' @export
#'
#' @examples
xyCoordinates <- function(inputMatrix) {
    indices <- which(inputMatrix == 1, arr.ind = TRUE)
    colnames(indices) <- c("x", "y")
    return(as.matrix(indices))
}

#' Function to normalize coodinates between zero and one while keep scaling
#'
#' @param coords matrix; matrix with coordinates
#'
#' @return matrix; coordinates scaled between 0 and 1
#' @export
#'
#' @examples
normalizeCoordinates <- function(coords) {
    # Calculate the range of x and y coordinates
    xRange <- max(coords[, 1]) - min(coords[, 1])
    yRange <- max(coords[, 2]) - min(coords[, 2])

    # Determine which axis is longer
    if (xRange >= yRange) {
        # Normalize x while maintaining the aspect ratio
        coords[, 1] <- (coords[, 1] - min(coords[, 1])) / xRange
        coords[, 2] <- (coords[, 2] - min(coords[, 2])) / xRange
    } else {
        # Normalize y while maintaining the aspect ratio
        coords[, 1] <- (coords[, 1] - min(coords[, 1])) / yRange
        coords[, 2] <- (coords[, 2] - min(coords[, 2])) / yRange
    }
    return(coords)
}



#' Function to get the dimension based on dim of y axis
#'
#' @param ppp point pattern object of class `ppp`
#' @param ydim dimension of y axis
#'
#' @return vector; vector with x and y dimension
#' @export
#'
#' @examples
getDimXY <- function(ppp, ydim) {
    xratio <- abs(diff(ppp$window$xrange)) / abs(diff(ppp$window$yrange))
    dimyx <- c(ydim, round(xratio * ydim))
    return(dimyx)
}

#' Function to convert spatial coordinates of a `SpatialExperiment` object to a `ppp` object
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds to the image
#' @param image_id character; image id, must be present in image_col
#'
#' @return ppp; object of type `ppp`
#' @export
#' @import SpatialExperiment
#'
#' @examples
SPE2ppp <- function(spe, marks,
    image_col = NULL,
    image_id = NULL) {
    if (!is.null(image_col) & !is.null(image_id)) {
        spe <- spe[, colData(spe)[[image_col]] == image_id]
    }

    ppp <- as.ppp(
        spatialCoords(spe),
        c(
            min(spatialCoords(spe)[, 1]),
            max(spatialCoords(spe)[, 1]),
            min(spatialCoords(spe)[, 2]),
            max(spatialCoords(spe)[, 2])
        )
    )
    marks(ppp) <- colData(spe)[[marks]]
    return(ppp)
}


#' Estimate the intensity threshold for the reconstruction of spatial strucutres
#'
#' @param ppp point pattern object of class `ppp`
#' @param bndw numeric; bandwith of the sigma parameter in the density estimation,
#' if no value is given the bandwith is estimated using cross validation with the `bw.diggle` function.
#' @param dimyx vector; pixel dimensions of the kernel density image
#' @param breaks numeric; number of bins of the intensity estimates
#'
#' @return numeric; estimated intensity threshold
#' @import spatstat.explore
#' @export
#'
#' @examples
findIntensityThreshold <- function(ppp, bndw = NULL, dimyx,
                                   breaks = 1000) {
  # define default of the sigma threshold
  if(is.null(bndw)) bndw <- bw.diggle(ppp)
  # create data frame
  den_df <- as.data.frame(density(ppp, sigma = bndw, dimyx = dimyx))
  # bin the value to get better estimates of the peaks
  cuts <- cut(den_df$value, breaks = breaks)
  # take all densities bigger than first bin
  # TODO: this has to be improved
  new_den <- density(den_df$value[cuts != levels(cuts)[1]])
  # define the peaks
  peaks <- new_den$x[which(diff(sign(diff(new_den$y))) == -2)]
  # for simplicity the threshold is the mean between the first and second peak,
  # we could in the future find a better estimation
  if(length(peaks) == 1){thres <- peaks}
  else {thres <- (peaks[2] - peaks[1]) / 2}
  return(thres)

}

