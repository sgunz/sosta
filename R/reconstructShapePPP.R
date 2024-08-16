#' Reconstruct structure from point pattern density
#'
#' @param ppp point pattern object of class `ppp`
#' @param bndw banddwith of kernel density estimator
#' @param thres intensity threshold for the reconstruction
#' @param dimyx pixel dimensions of the output image
#'
#' @return sf object of class `POLYGON`
#' @import spatstat.explore spatstat.geom sf
#' @export
#'
#' @examples
reconstructShapeDensity <- function(ppp,
    bndw,
    thres,
    dimyx) {
    # estimate density
    density_image <- density.ppp(ppp,
        bndw,
        dimyx = c(dimyx)
    )

    # TODO:
    # Error when threshold is too small

    # construct spatstat window from matrix with true false entries
    mat <- ifelse(t(as.matrix(density_image)) > thres, TRUE, FALSE)

    # using custom function
    stCast <- st_cast(st_make_valid(
      binaryImageToSF(
        mat,
        xmin = ppp$window$xrange[1],
        xmax = ppp$window$xrange[2],
        ymin = ppp$window$yrange[1],
        ymax = ppp$window$yrange[2]
      )
    )
    , "POLYGON") # make valid is important!!
    stCast <- stCast[!st_is_empty(stCast), drop = FALSE]

    # return sf object
    return(st_sf(st_cast(stCast, "POLYGON")))
}

#' Reconstruct structure from point pattern using likelihood ratio test
#'
#' @param ppp point pattern object of class `ppp`
#' @param r radius of the circle to use in the function `scanLRTS`
#' @param dimyx pixel dimension of the output image
#'
#' @return sf object of class `POLYGON`
#' @import spatstat.explore spatstat.geom sf
#' @importFrom EBImage resize
#' @export
#'
#'
#' @examples
reconstructShapeLRT <- function(ppp,
    r = bw.diggle(ppp),
    dimyx) {
    # likelihood ratio test
    LR <- scanLRTS(ppp,
        r,
        dimyx = c(dimyx)
    )

    # construct spatstat window from matrix with true false entries
    diff <- ceiling(c(dim(LR) - dimyx)[1] / 2)

    LRsub <- LR[
        diff:(dim(LR)[1] - diff),
        diff:(dim(LR)[2] - diff)
    ]

    matResc <- resize(t(as.matrix(LRsub)),
        w = diff(ppp$window$xrange),
        h = diff(ppp$window$yrange)
    )

    # TODO: Problem if matrix is too big, add some scaling factor
    # TODO: determine scaling factor
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
    stWinResc <- st_buffer(stWinResc, -0.5 * r)
    # remove empty polygons
    stWinResc <- stWinResc[!st_is_empty(stWinResc), drop = FALSE]
    # this ensures that we have single polygons
    return(st_cast(st_union(stWinResc), "POLYGON"))
}




#' Get intensity plot and intensity from spe object with given image id
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds to the image
#' @param image_id character; image id, must be present in image_col
#' @param mark_select character; name of mark that is to be selected for the reconstruction
#' @param bndw numeric; bandwith of the sigma parameter in the density estimation,
#' if no value is given the bandwith is estimated using cross validation with the `bw.diggle` function.
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction. Default = 500
#' @return ggplot object with intensity image and histogram
#' @importFrom ggplot2 ggplot aes geom_histogram theme_light geom_tile
#' labs coord_equal theme_classic scale_color_viridis_c
#' @export
#'
#' @examples
shapeIntensityImage <- function(spe, marks,
                              image_col,
                              image_id,
                              mark_select,
                              bndw = NULL,
                              dim = 500){

  # Convert the spe object to a point pattern object
  pp <- SPE2ppp(spe, marks = marks, image_col = image_col, image_id = image_id)

  # Extract the islet cells
  pp.islet <- subset(pp, marks %in% mark_select)

  # Set the dimensions of the resulting reconstruction
  dimyx <- getDimXY(pp.islet, dim)

  # Set default of sigma bandwith
  if (is.null(bndw)) bndw <- bw.diggle(pp.islet)

  # plot the density of the image
  im_df <- density(pp.islet, sigma = bndw, dimyx = dimyx) |> as.data.frame()

  # plot density image
  den_im <- im_df |>
    ggplot(aes(x = x, y = y, color = value)) +
    geom_tile() +
    scale_color_viridis_c(option = "C") +
    coord_equal() +
    labs(color = "intensity") +
    theme_classic()

  # plot histogram
  den_hist <- im_df |>
    ggplot(aes(x = value)) +
    geom_histogram(bins = 50) +
    labs(x = "pixel intensity") +
    theme_light()


  p <- patchwork::wrap_plots(den_im, den_hist, ncol = 2)

  return(p)
}


#' Reconstruct structure from spe object with given image id
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds to the image
#' @param image_id character; image id, must be present in image_col
#' @param mark_select character; name of mark that is to be selected for the reconstruction
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction. Default = 500
#' @param bndw numeric; bandwith of the sigma parameter in the density estimation,
#' if no value is given the bandwith is estimated using cross validation with the `bw.diggle` function.
#' @param thres numeric; intensity threshold for the reconstruction
#'
#' @return sf object of class `POLYGON`
#' @export
#'
#' @examples
reconstructShapeDensityImage <- function(spe, marks,
                                         image_col, image_id, mark_select,
                                         dim = 500, bndw = NULL,
                                         thres = NULL) {

  # Convert the spe object to a point pattern object
  pp <- SPE2ppp(spe,
                marks,
                image_col,
                image_id)

  # Extract the islet cells
  pp_sel <- subset(pp, marks %in% mark_select)

  # Get dimension of reconstruction
  dimyx <- getDimXY(pp_sel, dim)

  # Set default of sigma bandwith
  if (is.null(bndw)) bndw <- bw.diggle(pp_sel)

  # Find the threshold for the structure if not set
  if (is.null(thres)) thres <- findIntensityThreshold(pp_sel,
                                                      bndw = bndw,
                                                      dimyx = dimyx)
  # Get the structure
  struct <- reconstructShapeDensity(pp_sel,
                                    bndw = bndw,
                                    thres = thres,
                                    dimyx = dimyx)

  return(struct)
}


#' Reconstruct structure from spatial experiment object per image id
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds to the image
#' @param mark_select character; name of mark that is to be selected for the reconstruction
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction. Default = 500
#' @param bndw numeric; bandwith of the sigma parameter in the density estimation,
#' if no value is given the bandwith is estimated using cross validation with the `bw.diggle` function.
#' @param thres numeric; intensity threshold for the reconstruction
#' @param ncores numeric; number of cores for parallel processing using `mclapply`. Default = 1
#'
#' @importFrom parallel mclapply
#' @importFrom dplyr bind_rows
#'
#' @return simple feature collection
#' @export
#'
#' @examples
reconstructShapeDensitySPE <- function(spe, marks,
                                       image_col, mark_select,
                                       dim = 500, bndw = NULL, thres,
                                       ncores = 1){
  # Get all unique image ids
  all_images <- spe[[image_col]] |> unique()
  # Calculate polygon for each id using multiple cores
  res_all <- mclapply(all_images, function(x){
    res <- reconstructShapeDensityImage(spe, marks, image_col,
                                        x, mark_select,
                                        dim = 500, bndw = NULL,
                                        thres)
    # assign image_id
    res[[image_col]] <- x
    return(res)
  }, mc.cores = ncores
  )
  # return data frame with all structures
  # TODO: do.call(rbind, a)
  return(bind_rows(res_all))
}

#TODO: for one image

#' Title
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds to the image
#' @param mark_select character; name of mark that is to be selected for the reconstruction
#' @param nimages integer; number of images for the estimation. Will be randomly sampled
#' @param fun character; function to estimate the kernel density. Default bw.diggle.
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction. Default = 500
#' @param ncores numeric; number of cores for parallel processing using `mclapply`. Default = 1
#' @param plot_hist logical; if histogram of estimated densities and thresholds should be plotted. Default = TRUE
#'
#' @importFrom parallel mclapply
#' @importFrom dplyr bind_rows
#' @import spatstat.explore
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggplot aes geom_histogram theme_light
#' @return list; list of estimated intensities
#' @export
#'
#' @examples
estimateReconstructionParametersSPE <- function(spe,
                                     marks,
                                     image_col,
                                     mark_select = NULL,
                                     nimages = NULL,
                                     fun = "bw.diggle",
                                     dim = 500,
                                     ncores = 1,
                                     plot_hist = TRUE) {

  # get the id's of all images
  all_images <- colData(spe)[[image_col]] |> unique()
  # default is to take all values
  if(is.null(nimages)) nimages <- length(all_images)
  # alternitavely we sample some images
  sample_images <- sample(all_images, nimages)
  # we calculate the bandwiths and thresholds
  res <- parallel::mclapply(sample_images, function(x) {
    pp <- SPE2ppp(spe,
                  marks = marks,
                  image_col = image_col,
                  image_id = x)
    # If selection of mark
    if(!is.null(mark_select)) pp <- subset(pp, marks %in% mark_select)
    # Estimate the bandwidth for the kernel estimation of point process intensity
    bndw <- do.call(fun, args = list(X = pp, warn = FALSE))

    # Estimate the theshold for the recinstruction
    # Get dimension of reconstruction
    dimyx <- getDimXY(pp, dim)
    thres <- findIntensityThreshold(pp, bndw = bndw, dimyx = dimyx)

    return(c("img" = x, "bndw" = as.numeric(bndw), "thres" = as.numeric(thres)))
  }, mc.cores = ncores)

  #TODO: do.call(rbind, a) ?
  res <- dplyr::bind_rows(res)
  res$thres <- as.numeric(res$thres)
  res$bndw <- as.numeric(res$bndw)

  if (plot_hist == TRUE & nimages > 1){

    p1 <- res |>
      ggplot(aes(x = bndw)) +
      geom_histogram(bins = nimages/2) +
      theme_light()

    p2 <- res |>
      ggplot(aes(x = thres)) +
      geom_histogram(bins = nimages/2) +
      theme_light()

    print(patchwork::wrap_plots(p1, p2, ncol = 2))
  }

  return(res)
}
