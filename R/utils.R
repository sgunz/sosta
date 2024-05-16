

# Custom color palette ----------------------------------------------------

tramCol <- c('#e3000d', '#4b3b90', '#ffd500', '#119fe4',
             '#019e3b', '#a01c65', '#89d0e3', '#e92e89',
             '#945a24', '#afc903', '#dd9e42')


# -------------------------------------------------------------------------

#' Plot binary matrix
#'
#' @param inputMatrix binary matrix
#' @param main title of the plot
#'
#' @return ggplot object
#'
#' @export
#'
#' @examples
#
plotSim <- function(inputMatrix, main = ''){
  labeledMat <- EBImage::bwlabel(inputMatrix)
  # create data frame to store index
  xyIndexed <- data.frame(row = integer(),
                          col = integer(),
                          index = factor())
  # loop to get x y coordinates
  for (i in c(1:max(labeledMat))) {

    dfLabel <- data.frame(which(labeledMat == i, arr.ind = TRUE))
    dfLabel$index <- i

    xyIndexed <- rbind(xyIndexed, dfLabel)

  }

  ratio <- dim(inputMatrix)[2] / dim(inputMatrix)[1]

  #plot using ggplot
  p <- xyIndexed |>
    dplyr::mutate(index = as.factor(index)) |>
    ggplot(aes(x = row, y = col)) +
      geom_tile(aes(fill = index), show.legend = FALSE) +
      theme(aspect.ratio = ratio) +
      theme_void() +
      labs(title = main) +
      scale_fill_manual(values = c(rep(tramCol,
                                       ceiling(max(labeledMat)/10))))
  return(p)
}


#' Plot sf object
#'
#' @param sfImage sf object
#'
#' @return ggplot object
#' @export
#'
#' @examples
st_plotSim <- function(sfImage){
  # cast substructures
  sfImage_cast <- sf::st_as_sf(sf::st_cast(sfImage,  "POLYGON"))
  # number of strucutres
  nstruct <- dim(sfImage_cast)[1]
  # id for each substructure
  sfImage_cast$id <- factor(c(1:nstruct))
  # window of sf
  win <- sf::st_bbox(sfImage_cast)
  ratio <- abs(win$xmin - win$xmax) / abs(win$ymin - win$ymax)
  # plot
  sfImage_cast |>
    ggplot() +
      geom_sf(aes(fill = id, color = id), show.legend = FALSE) +
      theme_void() +
      theme(aspect.ratio = ratio) +
      scale_fill_manual(values = c(rep(tramCol, ceiling(max(nstruct)/10)))) +
      scale_color_manual(values = c(rep(tramCol, ceiling(max(nstruct)/10))))
}



# -------------------------------------------------------------------------


#' Converts a binary matrix to an sf polygon
#'
#' @param binaryMatrix binary matrix
#'
#' @return sf object
#' @import terra
#' @export
#'
#' @examples
binaryImageToSF <- function(binaryMatrix){
  # turn 90 degrees anti clockwise for correspondance with spatstat
  binaryMatrix <- apply(t(binaryMatrix),2,rev)
  # get raster
  r <- terra::rast(binaryMatrix)
  # convert to polygons
  poly <- terra::as.polygons(r)
  # polygons is a SpatVector. Convert it to an sf object
  polygonsSF <- sf::st_as_sf(poly)
  # Merge polygons to a single multipolygon
  return(sf::st_union(polygonsSF[polygonsSF$lyr.1 == 1,]))
}

# function to extract x y coordinates from binary image
xyCoordinates <- function(inputMatrix) {
  indices <- which(inputMatrix == 1, arr.ind = TRUE)
  colnames(indices) <- c('x', 'y')
  return(as.matrix(indices))
}

# function to normalize coodinates to [0,1] while keep scaling
normalizeCoordinates <- function(coords) {
  # Calculate the range of x and y coordinates
  xRange <- max(coords[,1]) - min(coords[,1])
  yRange <- max(coords[,2]) - min(coords[,2])

  # Determine which axis is longer
  if (xRange >= yRange) {
    # Normalize x while maintaining the aspect ratio
    coords$x <- (coords$x - min(coords$x)) / xRange
    coords$y <- (coords$y - min(coords$y)) / xRange
  } else {
    # Normalize y while maintaining the aspect ratio
    coords$x <- (coords$x - min(coords$x)) / yRange
    coords$y <- (coords$y - min(coords$y)) / yRange
  }
  return(coords)
}

binaryToLogical <- function(binaryMatrix) {
  # Convert binary matrix to logical
  logicalMatrix <- binaryMatrix >= 1
  return(logicalMatrix)
}

logicalToBinary <- function(Matrix) {
  # Convert binary matrix to logical
  binaryMatrix <- ifelse(Matrix == TRUE, 1, 0)
  return(binaryMatrix)
}


# -------------------------------------------------------------------------


edgeDetectionConvolution <- function(inputMatrix, filter =
                                       matrix(c(1,1,1,1,-8,1,1,1,1),
                                              nrow = 3, ncol = 3)) {

  # get the edges by convolution using high pass filter
  convMatrix <- OpenImageR::convolution(inputMatrix, filter)
  # binarize matrix
  convMatrix[convMatrix > 1] <- 1
  convMatrix[convMatrix < 0] <- 0
  # perimeter as sum of boundary pixels
  return(convMatrix)

}
