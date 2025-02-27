#' Simulate Tissue Blobs
#'
#' This function generates a simulated tissue-like structure using a Gaussian blur technique.
#'
#' @param size Integer; The size (width and height) of the simulated tissue image.
#' @param seedNumber Integer; The number of random seed points used to generate tissue blobs.
#' @param clumpSize Numeric; The standard deviation (sigma) of the Gaussian blur applied to generate tissue clumps.
#'
#' @return A binary matrix representing the simulated tissue structure.
#'
#' @importFrom EBImage gblur
#' @importFrom stats runif
#'
#' @examples
#' tissueImage <- simulateTissueBlobs(128, 100, 7)
#' image(tissueImage)
#'
#' @export
simulateTissueBlobs <- function(size, seedNumber, clumpSize) {
    # Create an initial random image
    image <- matrix(rep(0, size * size), nrow = size)

    # Generate random indices
    indices <- sample(size * size, seedNumber)

    # Set these indices to 1
    image[indices] <- 1

    # Use a Gaussian filter to create the tissue clumps
    image <- gblur(image, sigma = clumpSize)

    # Threshold the final image to get a binary image
    image <- image > mean(image)

    return(image)
}

#' Create a Point Pattern on a Simulated Tissue Image
#'
#' This function generates a spatial point pattern with different types of points (`A`, `B`, `C`) distributed over the simulated tissue structure.
#'
#' @param tissueImage Matrix; A binary matrix representing the simulated tissue.
#' @param intA Numeric; Intensity of type "A" points (points per unit area) on tissue regions.
#' @param intB Numeric; Intensity of type "B" points (points per unit area) on non-tissue regions.
#' @param intCInA Numeric; Intensity of type "C" points placed in extended regions around tissue.
#' @param intCInB Numeric; Intensity of type "C" points placed within tissue.
#'
#' @return A `ppp` object representing the spatial point pattern.
#'
#' @importFrom spatstat.geom owin as.rectangle grow.rectangle marks superimpose
#' @importFrom spatstat.random rpoispp
#'
#' @examples
#' tissueImage <- simulateTissueBlobs(128, 100, 7)
#' createPointPatternTissue(tissueImage, 0.1, 0.1, 0.005, 0.005)
#'
#' @export
createPointPatternTissue <- function(tissueImage, intA, intB, intCInA, intCInB) {
    # Create a binary image of non-tissue
    nonTissue <- (tissueImage == 0)

    # Convert binary image to window
    tissueWindow <- owin(mask = t(tissueImage))
    extendedWindow <- grow.rectangle(as.rectangle(tissueWindow))
    nonTissueWindow <- owin(mask = t(nonTissue))

    # Create point pattern with noise
    pointsA <- rpoispp(intA, win = tissueWindow)
    aNoise <- rpoispp(intA / 20, win = extendedWindow)
    pointsB <- rpoispp(intB, win = nonTissueWindow)
    pointsC <- rpoispp(intCInA, win = extendedWindow)
    pointsC2 <- rpoispp(intCInB, win = tissueWindow)

    # Assign marks (labels) to points
    spatstat.geom::marks(pointsA) <- "A"
    spatstat.geom::marks(aNoise) <- "A"
    spatstat.geom::marks(pointsB) <- "B"
    spatstat.geom::marks(pointsC) <- "C"
    spatstat.geom::marks(pointsC2) <- "C"

    # Combine all points into a single pattern
    pointPattern <- superimpose(pointsA, aNoise, pointsB, pointsC, pointsC2,
        W = extendedWindow
    )

    return(pointPattern)
}
