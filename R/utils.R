#' Converts a binary matrix to an sf polygon
#'
#' @param binaryMatrix matrix; binary matrix
#' @param xmin integer; minimum x coordinate of the coordinate system
#' @param xmax integer; maximum x coordinate of the coordinate system
#' @param ymin integer; minimum y coordinate of the coordinate system
#' @param ymax integer; maximum y coordinate of the coordinate system
#'
#' @return sf object
#' @importFrom sf st_as_sf st_union
#' @importFrom terra rast as.polygons set.ext
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
#' plot(poly_R)
binaryImageToSF <- function(
        binaryMatrix,
        xmin, xmax,
        ymin, ymax) {
    # Input checking
    stopifnot("'binaryMatrix' must be a matrix" = is.matrix(binaryMatrix))
    stopifnot(
        "'xmin', 'xmax', 'ymin', and 'ymax' must be numeric" =
            is.numeric(c(xmin, xmax, ymin, ymax))
    )
    stopifnot("'xmin' must be less than 'xmax'" = xmin < xmax)
    stopifnot("'ymin' must be less than 'ymax'" = ymin < ymax)
    # turn 90 degrees anti clockwise for correspondance with spatstat
    binaryMatrix <- apply(t(binaryMatrix), 2, rev)
    # get raster
    r <- rast(binaryMatrix)
    # rescale to correct windwow
    set.ext(r, c(xmin, xmax, ymin, ymax))
    # convert to polygons
    poly <- as.polygons(r) # Why does it print "hardcopy" here?
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
#' xyCoordinates(matrix_R)
xyCoordinates <- function(inputMatrix) {
    # Input checking
    stopifnot("'inputMatrix' must be a matrix" = is.matrix(inputMatrix))
    # Code
    indices <- which(inputMatrix == 1, arr.ind = TRUE)
    colnames(indices) <- c("x", "y")
    return(as.matrix(indices))
}

#' Function to normalize coordinates between zero and one while keep scaling
#'
#' @param coords matrix; matrix with coordinates
#'
#' @return matrix; coordinates scaled between 0 and 1
#' @export
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
#' coords <- xyCoordinates(matrix_R)
#' normalizeCoordinates(coords)
normalizeCoordinates <- function(coords) {
    stopifnot("'coords' must be a matrix" = is.matrix(coords))
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
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' ppp <- SPE2ppp(spe,
#'     marks = "cell_category", image_col = "image_name",
#'     image_id = "E04"
#' )
#' getDimXY(ppp, 500)
getDimXY <- function(ppp, ydim) {
    # Input checking
    stopifnot("'ppp' must be an object of class 'ppp'" = inherits(ppp, "ppp"))
    stopifnot("'ydim' must be a single numeric value" = is.numeric(ydim) && length(ydim) == 1)

    xratio <- abs(diff(ppp$window$xrange)) / abs(diff(ppp$window$yrange))
    dimyx <- c(ydim, round(xratio * ydim))
    return(dimyx)
}

#' Function to convert spatial coordinates of a `SpatialExperiment` object to a `ppp` object
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to
#' the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds to
#' the image
#' @param image_id character; image id, must be present in image_col
#'
#' @return ppp; object of type `ppp`
#' @export
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @importFrom spatstat.geom as.ppp setmarks
#'
#' @examples
#' data(sostaSPE)
#' SPE2ppp(sostaSPE, marks = "cell_type", image_col = "image_name",
#' image_id = "image1")
SPE2ppp <- function(spe,
    marks,
    image_col = NULL,
    image_id = NULL) {
    # Input checking
    stopifnot(
        "'spe' must be an object of class 'SpatialExperiment'" =
            inherits(spe, "SpatialExperiment")
    )
    stopifnot(
        "'marks' must exist in colData(spe)" =
            (marks %in% colnames(colData(spe)))
    )

    # Subset the SPE object
    if (!is.null(image_col) & !is.null(image_id)) {
        stopifnot(
            "'image_col' must exist in colData(spe)" =
                image_col %in% colnames(colData(spe))
        )
        stopifnot(
            "'image_id' must exist in colData(spe)[['image_col']]" =
                image_id %in% colData(spe)[[image_col]]
        )
        spe <- spe[, colData(spe)[[image_col]] %in% image_id]
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
    ppp <- setmarks(ppp, colData(spe)[[marks]])
    return(ppp)
}


#' Estimate the intensity threshold for the reconstruction of spatial structures
#'
#' @param ppp point pattern object of class `ppp`
#' @param mark_select character; name of mark that is to be selected for the
#'  reconstruction
#' @param bndw numeric; bandwith of the sigma parameter in the density estimation,
#' if no value is given the bandwith is estimated using cross validation with
#' the `bw.diggle` function.
#' @param dim numeric; x dimension of the final reconstruction.
#' @param steps numeric; value used to filter the density estimates, where only
#' densities greater than the maximum value divided by \code{threshold} are considered.
#' Default is 250.
#'
#' @return numeric; estimated intensity threshold
#' @importFrom spatstat.explore bw.diggle density.ppp
#' @importFrom spatstat.geom subset.ppp
#' @importFrom stats density
#' @export
#'
#' @examples
#' data(sostaSPE)
#' ppp <- SPE2ppp(sostaSPE, marks = "cell_type", image_col = "image_name", image_id = "image1")
#' findIntensityThreshold(ppp, mark_select = "A", dim = 250)
findIntensityThreshold <- function(ppp, mark_select = NULL,
    bndw = NULL, dim,
    steps = 250) {
    stopifnot("'steps' must be a single numeric value" = is.numeric(dim) && length(dim) == 1)
    # get density image
    density_image <- .intensityImage(ppp, mark_select, bndw, dim)$den_im
    # calculate threshold
    thres <- .intensityThreshold(density_image, steps)
    return(thres)
}


#' Function to estimate the intensity image of a point pattern
#' @param ppp point pattern object of class `ppp`
#' @param mark_select character; name of mark that is to be selected for the
#'  reconstruction
#' @param bndw bandwidth of kernel density estimator
#' @param dim numeric; x dimension of the final reconstruction.
#'
#' @return list; list with the intensity image and the bandwidth and dimension parameters
#' @importFrom spatstat.explore bw.diggle density.ppp
#' @importFrom spatstat.geom subset.ppp
.intensityImage <- function(ppp,
    mark_select = NULL,
    bndw = NULL,
    dim) {
    # Input checking
    stopifnot("'ppp' must be an object of class 'ppp'" = inherits(ppp, "ppp"))
    stopifnot("'dim' must be a single, positive, numeric value" = is.numeric(dim) &&
        length(dim) == 1 & dim > 0)

    if (!is.null(bndw)) {
        stopifnot("'bndw' must be a single numeric value" = is.numeric(bndw) && length(bndw) == 1)
    }

    # Extract the islet cells
    if (!is.null(mark_select)) {
        stopifnot(
            "All values in 'mark_select' must exist in 'marks' of 'ppp'" =
                all(mark_select %in% marks(ppp))
        )
        pp_sel <- subset.ppp(ppp, marks %in% mark_select)
    } else {
        pp_sel <- ppp
    }


    # Set the dimensions of the resulting reconstruction
    dimyx <- getDimXY(pp_sel, dim)

    # Set default of sigma bandwith
    if (is.null(bndw)) bndw <- bw.diggle(pp_sel)

    # plot the density of the image
    den <- density.ppp(pp_sel,
        sigma = bndw,
        dimyx = dimyx,
        positive = TRUE
    )

    return(list(den_im = den, bndw = bndw, dimyx = dimyx))
}

#' Function to estimate the intensity threshold for the reconstruction of spatial structures
#'
#' @param density_image real-valued pixel image; output from the function `.intensityImage`
#' @param steps numeric; value used to filter the density estimates, where only
#' densities greater than the maximum value divided by \code{threshold} are considered.
#' Default is 250.
#'
#' @return numeric; estimated threshold
#' @importFrom stats density
.intensityThreshold <- function(density_image, steps = 250) {
    # take all densities greater than certain threshold due to numerical properties
    # of the density estimation
    den_df <- density_image |> as.data.frame()
    new_den <- density(den_df$value[den_df$value > max(den_df$value) / steps])
    # define the peaks x values
    peaks <- new_den$x[which(diff(sign(diff(new_den$y))) == -2)]
    # define peak values
    peak_vals <- new_den$y[which(diff(sign(diff(new_den$y))) == -2)]
    # the threshold is the mean between the two main modes of the distribution
    if (length(peaks) == 1) {
        thres <- peaks
    } else {
        thres <- (peaks[order(peaks, decreasing = FALSE)[2]] -
            peaks[order(peaks, decreasing = FALSE)[1]]) / 2 +
            peaks[order(peaks, decreasing = FALSE)[1]]
    }
    return(thres)
}

#' Function to convert spatialCoords to an sf object
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#'
#' @returns sf; Simple feature collection of geometry type POINT
#' @importFrom sf  st_as_sf
#' @importFrom SpatialExperiment  spatialCoords
#'
#' @examples
#' data(sostaSPE)
#' spe_sel <- sostaSPE[, sostaSPE[["image_name"]] == "image1"]
#' spatialCoords2SF(spe_sel)
#' @export
spatialCoords2SF <- function(spe){
    # Input checking
    stopifnot(
        "'spe' must be an object of class 'SpatialExperiment'" =
            inherits(spe, "SpatialExperiment")
    )
    # creates sf points object from SPE
    spatial_coords_sf <- st_as_sf(data.frame(spatialCoords(spe)),
                                  coords = c(colnames(spatialCoords(spe))[1],
                                             colnames(spatialCoords(spe))[2]))
    return(spatial_coords_sf)
}





