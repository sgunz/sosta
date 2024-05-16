#' Reconstruct shape from point pattern density
#'
#' @param ppp point pattern object of class `ppp`
#' @param bndw banddwith of kernel density estimator
#' @param thres intensity threshold
#' @param dimyx pixel dimension of the output image
#'
#' @return sf object of class `POLYGON`
#' @import spatstat.explore spatstat.geom sf
#' @export
#'
#' @examples

reconstructShapeDensity <- function(ppp,
                                    bndw,
                                    thres,
                                    dimyx){
  # estimate density
  density_image <- density.ppp(ppp,
                               bndw,
                               dimyx = c(dimyx))

  # TODO:
  # Error when threshold is too small

  # construct spatstat window from matrix with true false entries
  mat <- ifelse(t(as.matrix(density_image)) > thres, TRUE, FALSE)
  shapeWin <- owin(mask = mat)
  # rescale based on original ppp
  rescaleShapeWin <- rescale.owin(shapeWin,
                                  diff(shapeWin$xrange)/diff(ppp$window$xrange))

  #TODO: rescaling directly using sf objects
  # using custom function
  stCast <- st_cast(st_make_valid(binaryImageToSF(as.matrix(rescaleShapeWin))), "POLYGON")  # make valid is important!!
  stCast <- stCast[!st_is_empty(stCast), drop = FALSE]
  # return sf object
  return(st_sf(st_cast(stCast, "POLYGON")))
}

#' Reconstruction of shape from point pattern using likelihood ratio test
#'
#' @param ppp point pattern object of class `ppp`
#' @param r radius of the circle to use in the function `scanLRTS`
#' @param dimyx pixel dimension of the output image
#'
#' @return
#' @import spatstat.explore spatstat.geom sf
#' @importFrom EBImage resize
#' @export
#'
#'
#' @examples
reconstructShapeLRT <- function(ppp,
                                 r = bw.diggle(ppp),
                                 dimyx){
  # likelihood ratio test
  LR <- scanLRTS(ppp,
                 r,
                 dimyx = c(dimyx))

  # construct spatstat window from matrix with true false entries
  diff <- ceiling(c(dim(LR) - dimyx)[1]/2)

  LRsub <- LR[diff:(dim(LR)[1] - diff),
              diff:(dim(LR)[2] - diff)]

  matResc <- resize(t(as.matrix(LRsub)), w = diff(ppp$window$xrange),
                    h = diff(ppp$window$yrange))

  #TODO: Problem if matrix is too big, add some scaling factor
  #TODO: determine scaling factor
  # matResc <- resize(t(as.matrix(LRsub)), w = diff(tonsilPPP$window$xrange)/100,
  #                   h = diff(tonsilPPP$window$yrange)/100)
  #
  # shapeWinResc <- owin(mask = t(ifelse(matResc > 0, TRUE, FALSE)))
  #
  # rescaleShapeWin <- rescale.owin(shapeWinResc,
  #                                 1/100)



  shapeWinResc <- owin(mask = t(ifelse(matResc > 0, TRUE, FALSE)))
  # cast to sf object
  stWinResc <- st_cast(st_as_sf(shapeWinResc), "POLYGON")
  # return remove 0.5 * r, to account for radius
  stWinResc <- st_buffer(stWinResc, -0.5*r)
  # remove empty polygons
  stWinResc <- stWinResc[!st_is_empty(stWinResc),drop=FALSE]
  # this ensures that we have single polygons
  return(st_cast(st_union(stWinResc), "POLYGON"))
}
