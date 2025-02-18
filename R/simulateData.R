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
#' tissue_image <- simulateTissueBlobs(128, 100, 7)
#' image(tissue_image)
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
#' @param tissue_image Matrix; A binary matrix representing the simulated tissue.
#' @param int_a Numeric; Intensity of type "A" points (points per unit area) on tissue regions.
#' @param int_b Numeric; Intensity of type "B" points (points per unit area) on non-tissue regions.
#' @param int_c_in_a Numeric; Intensity of type "C" points placed in extended regions around tissue.
#' @param int_c_in_b Numeric; Intensity of type "C" points placed within tissue.
#'
#' @return A `ppp` object representing the spatial point pattern.
#'
#' @importFrom spatstat.geom owin as.rectangle grow.rectangle marks superimpose
#' @importFrom spatstat.random rpoispp
#'
#' @examples
#' tissue_image <- simulateTissueBlobs(128, 100, 7)
#' point_pattern <- createPointPatternTissue(tissue_image, 0.1, 0.1, 0.005, 0.005)
#'
#' @export
createPointPatternTissue <- function(tissue_image, int_a, int_b, int_c_in_a, int_c_in_b) {
    # Create a binary image of non-tissue
    non_tissue <- (tissue_image == 0)

    # Convert binary image to window
    tissue_window <- owin(mask = t(tissue_image))
    extended_window <- grow.rectangle(as.rectangle(tissue_window))
    non_tissue_window <- owin(mask = t(non_tissue))

    # Create point pattern with noise
    points_a <- rpoispp(int_a, win = tissue_window)
    a_noise <- rpoispp(int_a / 20, win = extended_window)
    points_b <- rpoispp(int_b, win = non_tissue_window)
    points_c <- rpoispp(int_c_in_a, win = extended_window)
    points_c2 <- rpoispp(int_c_in_b, win = tissue_window)

    # Assign marks (labels) to points
    spatstat.geom::marks(points_a) <- "A"
    spatstat.geom::marks(a_noise) <- "A"
    spatstat.geom::marks(points_b) <- "B"
    spatstat.geom::marks(points_c) <- "C"
    spatstat.geom::marks(points_c2) <- "C"

    # Combine all points into a single pattern
    point_pattern <- superimpose(points_a, a_noise, points_b, points_c, points_c2,
                                 W = extended_window)

    return(point_pattern)
}
