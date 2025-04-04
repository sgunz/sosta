#' Reconstruct polygon from point pattern density
#'
#' This function estimates the density of a spatial point pattern (`ppp`),
#' thresholds the density to create a binary image, and then converts it
#' to a valid `sf` object (polygons).
#'
#' @param ppp point pattern object of class `ppp`
#' @param markSelect character; name of mark that is to be selected for the
#'  reconstruction
#' @param bndw bandwidth of kernel density estimator
#' @param thres intensity threshold for the reconstruction
#' @param dim numeric; x dimension of the final reconstruction.
#'
#' @return sf object of class `POLYGON`
#' @importFrom sf st_cast st_make_valid st_sf st_is_empty st_geometry
#' @export
#'
#' @examples
#' data("sostaSPE")
#' ppp <- SPE2ppp(sostaSPE, marks = "cellType", imageCol = "imageName", imageId = "image1")
#' thres <- findIntensityThreshold(ppp, markSelect = "A", dim = 500)
#' struct <- reconstructShapeDensity(ppp, markSelect = "A", thres = thres, dim = 500)
#' plot(struct)
reconstructShapeDensity <- function(ppp, markSelect = NULL,
    bndw = NULL, thres = NULL, dim) {
    # estimate density
    res <- .intensityImage(ppp, markSelect, bndw, dim)

    if (!is.null(thres)) {
        stopifnot("'thres' must be a single numeric value" = is.numeric(thres) &&
            length(thres) == 1)
    }

    # Check if intensity threshold exists
    if (is.null(thres)) {
        thres <- .intensityThreshold(res$denIm)
    }

    # construct spatstat window from matrix with true false entries
    mat <- ifelse(t(as.matrix(res$denIm)) > thres, TRUE, FALSE)

    # Check if we get empty or full polygon
    if (all(mat == 1)) {
        warning("Full image converted to polygon; threshold might be too low")
    }

    if (all(mat == 0)) {
        warning("No structure found; threshold might be too high")
        return()
    }

    # using custom function
    stCast <- st_cast(
        st_make_valid(
            binaryImageToSF(
                mat,
                xmin = ppp$window$xrange[1], xmax = ppp$window$xrange[2],
                ymin = ppp$window$yrange[1], ymax = ppp$window$yrange[2]
            )
        ),
        "POLYGON"
    ) # make valid is important
    stCast <- stCast[!st_is_empty(stCast), drop = FALSE]

    obj <- st_sf(st_cast(stCast, "POLYGON"))
    sf::st_geometry(obj) <- "sostaPolygon"
    return(obj)
}


#' Intensity plot
#'
#' This function plots the intensity of a point pattern image and displays
#' a histogram of the intensity values. Note that intensities less than
#' largest intensity value divided by 250 are not displayed in the histogram.
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to
#' the `ppp` marks
#' @param imageCol character; name of a column in `colData` that corresponds to
#' the image
#' @param imageId character; image id, must be present in imageCol
#' @param markSelect character; name of mark that is to be selected for the
#' reconstruction
#' @param bndw numeric; smoothing bandwidth in the density estimation,
#' corresponds to the `sigma` parameter in the `density.ppp` function,
#' if no value is given the bandwidth is estimated using cross validation with
#' the `bw.diggle` function.
#' @param dim numeric; x dimension of the final reconstruction. A lower resolution
#' speeds up computation but lead to less exact reconstruction. Default = 500
#' @return ggplot object with intensity image and histogram
#' @importFrom ggplot2 ggplot aes_string geom_histogram theme_light geom_tile
#' labs coord_equal theme_classic scale_color_viridis_c geom_vline theme element_text
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom spatstat.geom subset.ppp
#' @export
#'
#' @examples
#' data("sostaSPE")
#' shapeIntensityImage(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     imageId = "image1", markSelect = "A"
#' )
shapeIntensityImage <- function(
        spe, marks,
        imageCol,
        imageId,
        markSelect,
        bndw = NULL,
        dim = 500) {
    # Convert the spe object to a point pattern object
    ppp <- SPE2ppp(spe, marks = marks, imageCol = imageCol, imageId = imageId)

    # plot the density of the image
    res <- .intensityImage(ppp, markSelect, bndw, dim)
    im_df <- res$denIm |> as.data.frame()
    thres <- findIntensityThreshold(ppp, markSelect, res$bndw, dim)

    # plot density image
    denIm <- im_df |>
        ggplot(aes(x = .data$x, y = .data$y, color = .data$value)) +
        geom_tile() +
        coord_equal() +
        labs(color = "intensity") +
        scale_color_viridis_c(option = "C") +
        theme_classic()

    # plot histogram
    den_hist <- im_df |>
        filter(.data$value > max(.data$value) / 250) |>
        ggplot(aes(x = abs(.data$value))) + # Use .data pronoun
        geom_histogram(bins = 50) +
        labs(x = "pixel intensity") +
        theme_light() +
        geom_vline(xintercept = thres, color = "seagreen")


    p <- wrap_plots(denIm, den_hist, ncol = 2) +
        plot_annotation(
            title = paste0(imageCol, ": ", imageId),
            subtitle = paste0(
                "bndw: ", round(res$bndw, 4), "; estimated thres: ",
                round(thres, 4)
            ),
            caption = paste0(
                "Dimension of the density image (pixels): ", res$dimyx[1],
                "x", res$dimyx[2]
            ),
            theme = theme(plot.caption = element_text(hjust = 0))
        )

    return(p)
}


#' Reconstruct structure from spe object with given image id
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond
#' to the `ppp` marks
#' @param imageCol character; name of a column in `colData` that corresponds
#' to the image
#' @param imageId character; image id, must be present in imageCol
#' @param markSelect character; name of mark that is to be selected for the
#'  reconstruction
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction.
#'  Default = 500
#' @param bndw numeric; smoothing bandwidth in the density estimation,
#' corresponds to the `sigma` parameter in the `density.ppp` function,
#' if no value is given the bandwidth is estimated using cross validation with
#' the `bw.diggle` function.
#' @param thres numeric; intensity threshold for the reconstruction;
#' if NULL the threshold is set as the mean between the mode of the pixel intensity
#' distributions
#' @return sf object of class `POLYGON`
#' @importFrom spatstat.geom subset.ppp
#' @export
#'
#' @examples
#' data("sostaSPE")
#' struct <- reconstructShapeDensityImage(sostaSPE,
#'     marks = "cellType", imageCol = "imageName", imageId = "image1",
#'     markSelect = "A", dim = 500
#' )
#' plot(struct)
reconstructShapeDensityImage <- function(spe, marks,
    imageCol, imageId, markSelect, dim = 500, bndw = NULL, thres = NULL) {
    # Convert the spe object to a point pattern object
    ppp <- SPE2ppp(spe, marks, imageCol, imageId)

    # Get the structure
    struct <- reconstructShapeDensity(ppp, markSelect, bndw, thres, dim)
    return(struct)
}


#' Reconstruct structure from spatial experiment object per image id
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond
#' to the `ppp` marks
#' @param imageCol character; name of a column in `colData` that corresponds
#' to the image
#' @param markSelect character; name of mark that is to be selected for the
#' reconstruction
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction.
#' Default = 500
#' @param bndw numeric; bandwidth of the sigma parameter in the density estimation,
#' if no value is given the bandwidth is estimated using cross validation with
#' the `bw.diggle` function for each image individually.
#' @param thres numeric; intensity threshold for the reconstruction;
#' if NULL the threshold is set as the mean between the mode of the pixel intensity
#' distributions estimated for each image individual
#' @param nCores numeric; number of cores for parallel processing using
#' `mclapply`. Default = 1
#'
#' @importFrom parallel mclapply
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom spatstat.geom as.ppp setmarks
#'
#' @return simple feature collection
#' @export
#'
#' @examples
#' data("sostaSPE")
#' allStructs <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     markSelect = "A", bndw = 3.5, thres = 0.005
#' )
#' allStructs
reconstructShapeDensitySPE <- function(spe, marks,
    imageCol, markSelect,
    dim = 500, bndw = NULL, thres = NULL,
    nCores = 1) {
    # Create a data frame with all necessary variables for computational (memory) reasons
    df <- .SPE2df(spe, imageCol, marks)
    # Remove SPE to free up memory
    rm(spe)
    gc()
    # Split up by image name
    ls <- split(df, as.factor(df[, 3]))
    # Calculate polygon for each id using multiple cores
    res_all <- mclapply(ls, function(x) {
        ppp <- .df2ppp(x)
        # Reconstruct the structure
        res <- reconstructShapeDensity(ppp, markSelect, bndw, thres, dim)
        if (is.null(res)) {
            message(paste0(
                "No structure found in: ",
                unique(x[, 3])
            ), " (see Warning below)")
            return()
        }
        # Assign imageId
        res[["structID"]] <- paste0(unique(x[, 3]), "_", c(1:dim(res)[1]))
        res[[imageCol]] <- unique(x[, 3])
        return(res)
    }, mc.cores = nCores)
    # Return data frame with all structures
    return(do.call(rbind, res_all))
}

#' Estimate reconstruction parameters from a set of images
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to
#' the `ppp` marks
#' @param imageCol character; name of a column in `colData` that corresponds
#' to the image
#' @param markSelect character; name of mark that is to be selected for the
#' reconstruction
#' @param nImages integer; number of images for the estimation. Will be randomly
#' sampled
#' @param fun character; function to estimate the kernel density. Default
#' bw.diggle.
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction.
#' Default = 500
#' @param nCores numeric; number of cores for parallel processing using `mclapply`.
#' Default = 1
#' @param plotHist logical; if histogram of estimated densities and thresholds
#' should be plotted. Default = TRUE
#'
#' @importFrom spatstat.geom subset.ppp
#' @importFrom SummarizedExperiment colData assays
#' @importFrom parallel mclapply
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggplot aes_string geom_histogram theme_light
#' @importFrom rlang .data
#'
#' @return tibble; tibble with estimated intensities and thresholds
#' @export
#'
#' @examples
#' data("sostaSPE")
#' estimateReconstructionParametersSPE(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     markSelect = "A", plotHist = TRUE
#' )
estimateReconstructionParametersSPE <- function(
        spe,
        marks,
        imageCol,
        markSelect = NULL,
        nImages = NULL,
        fun = "bw.diggle",
        dim = 500,
        nCores = 1,
        plotHist = TRUE) {
    # Input checks
    if (!is.null(nImages)) {
        stopifnot("'nImages' must be numeric" = is.numeric(nImages))
        stopifnot(
            "'nImages' must be smaller or equal to the number of images in the `SpatialExperiment`" =
                (nImages < length(unique(colData(spe)[[imageCol]])))
        )
    }

    # get the id's of all images
    allImages <- colData(spe)[[imageCol]] |> unique()
    # default is to take all values
    if (is.null(nImages)) nImages <- length(allImages)
    # alternatively we sample some images
    sampleImages <- sample(allImages, nImages)
    # Select sampled images
    spe <- spe[, colData(spe)[[imageCol]] %in% sampleImages]

    # Create a data frame with all necessary variables for computational (memory) reasons
    df <- .SPE2df(spe, imageCol, marks)
    # Remove SPE to free up memory
    rm(spe)
    gc()
    # Split up by image name
    ls <- split(df, as.factor(df[, 3]))

    # we calculate the bandwidths and thresholds
    res <- mclapply(ls, function(x) {
        ppp <- .df2ppp(x)
        res_x <- .intensityImage(ppp, markSelect, dim = dim)
        thres <- .intensityThreshold(res_x$denIm)
        return(list(img = x, bndw = as.numeric(res_x$bndw), thres = as.numeric(thres)))
    }, mc.cores = nCores)

    # collect in one data frame
    res <- as.data.frame(do.call(rbind, res))
    res$bndw <- as.numeric(res$bndw)
    res$thres <- as.numeric(res$thres)

    if (plotHist == TRUE & nImages > 1) {
        p1 <- res |>
            ggplot(aes(x = .data$bndw)) +
            geom_histogram(bins = round(nImages / 2)) +
            theme_light()

        p2 <- res |>
            ggplot(aes(x = .data$thres)) +
            geom_histogram(bins = round(nImages / 2)) +
            theme_light()

        plot(wrap_plots(p1, p2, ncol = 2))
    }

    return(res)
}
