#' Simulate Tissue Blobs
#'
#' @param size the size of the binary image
#' @param seedNumber number of seed points
#' @param clumpSize the size of the Gaussian Kernel
#'
#' @return A binary image with simulated tissue blobs
#' @export
#'
#' @examples
#'
simulateTissueBlobs <- function(size, seedNumber, clumpSize){
  # Create an initial random image
  image <- matrix(rep(0, size * size), nrow = size)

  # Generate random indices
  indices <- sample(size * size, seedNumber)

  # Set these indices to 1
  image[indices] <- 1

  # Use a Gaussian filter to create the tissue clumps
  image <- EBImage::gblur(image, sigma = clumpSize)

  # Threshold the final image to get a binary image
  image <- image > mean(image)

  return(image)
}
