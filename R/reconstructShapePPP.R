#' Reconstruct polygon from point pattern density
#'
#' This function estimates the density of a spatial point pattern (`ppp`),
#' thresholds the density to create a binary image, and then converts it
#' to a valid `sf` object (polygons).
#'
#' @param ppp point pattern object of class `ppp`
#' @param mark_select character; name of mark that is to be selected for the
#'  reconstruction
#' @param bndw bandwidth of kernel density estimator
#' @param thres intensity threshold for the reconstruction
#' @param dim numeric; x dimension of the final reconstruction.
#'
#' @return sf object of class `POLYGON`
#' @importFrom sf st_cast st_make_valid st_sf st_is_empty
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' ppp <- SPE2ppp(spe, marks = "cell_category", image_col = "image_name", image_id = "E04")
#' thres <- findIntensityThreshold(ppp, mark_select = "islet", dim = 500)
#' islet_poly <- reconstructShapeDensity(ppp, mark_select = "islet", thres = thres, dim = 500)
#' plot(islet_poly)
reconstructShapeDensity <- function(ppp, mark_select = NULL,
    bndw = NULL, thres = NULL, dim) {
    # estimate density
    res <- .intensityImage(ppp, mark_select, bndw, dim)

    if (!is.null(thres)) {
        stopifnot("'thres' must be a single numeric value" = is.numeric(thres) &&
            length(thres) == 1)
    }

    # Check if intensity threshold exists
    if (is.null(thres)) {
        thres <- .intensityThreshold(res$den_im)
    }

    # construct spatstat window from matrix with true false entries
    mat <- ifelse(t(as.matrix(res$den_im)) > thres, TRUE, FALSE)

    # Check if we get empty or full polygon
    stopifnot("Threshold too low" = (!all(mat == 1)))
    stopifnot("Threshold too high" = (!all(mat == 0)))

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

    # return sf object
    return(st_sf(st_cast(stCast, "POLYGON")))
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
#' @param image_col character; name of a column in `colData` that corresponds to
#' the image
#' @param image_id character; image id, must be present in image_col
#' @param mark_select character; name of mark that is to be selected for the
#' reconstruction
#' @param bndw numeric; bandwith of the sigma parameter in the density estimation,
#' if no value is given the bandwith is estimated using cross validation with
#' the `bw.diggle` function.
#' @param dim numeric; x dimension of the final reconstruction. A lower resolution
#' speeds up computation but lead to less exact reconstruction. Default = 500
#' @return ggplot object with intensity image and histogram
#' @importFrom ggplot2 ggplot aes_string geom_histogram theme_light geom_tile
#' labs coord_equal theme_classic scale_color_viridis_c geom_vline
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom spatstat.geom subset.ppp
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' shapeIntensityImage(spe,
#'     marks = "cell_category", image_col = "image_name",
#'     image_id = "E04", mark_select = "islet"
#' )
shapeIntensityImage <- function(spe, marks,
    image_col,
    image_id,
    mark_select,
    bndw = NULL,
    dim = 500) {
    # Convert the spe object to a point pattern object
    ppp <- SPE2ppp(spe, marks = marks, image_col = image_col, image_id = image_id)

    # plot the density of the image
    res <- .intensityImage(ppp, mark_select, bndw, dim)
    im_df <- res$den_im |> as.data.frame()
    thres <- findIntensityThreshold(ppp, mark_select, res$bndw, dim)

    # plot density image
    den_im <- im_df |>
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


    p <- wrap_plots(den_im, den_hist, ncol = 2) +
        plot_annotation(
            title = paste0(image_col, ": ", image_id),
            subtitle = paste0(
                "bndw: ", round(res$bndw, 4), "; estimated thres: ",
                round(thres, 4)
            ),
            caption = paste0(
                "Pixel image dimensions: ", res$dimyx[1],
                "x", res$dimyx[2]
            )
        )

    return(p)
}


#' Reconstruct structure from spe object with given image id
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond
#' to the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds
#' to the image
#' @param image_id character; image id, must be present in image_col
#' @param mark_select character; name of mark that is to be selected for the
#'  reconstruction
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction.
#'  Default = 500
#' @param bndw numeric; bandwith of the sigma parameter in the density estimation,
#' if no value is given the bandwith is estimated using cross validation with
#' the `bw.diggle` function.
#' @param thres numeric; intensity threshold for the reconstruction
#' @return sf object of class `POLYGON`
#' @importFrom spatstat.geom subset.ppp
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' islet_poly <- reconstructShapeDensityImage(spe,
#'     marks = "cell_category",
#'     image_col = "image_name", image_id = "E04", mark_select = "islet", dim = 500
#' )
#' plot(islet_poly)
reconstructShapeDensityImage <- function(
        spe, marks,
        image_col, image_id, mark_select, dim = 500, bndw = NULL, thres = NULL) {
    # Convert the spe object to a point pattern object
    ppp <- SPE2ppp(spe, marks, image_col, image_id)

    # Get the structure
    struct <- reconstructShapeDensity(ppp, mark_select, bndw, thres, dim)

    return(struct)
}


#' Reconstruct structure from spatial experiment object per image id
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond
#' to the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds
#' to the image
#' @param mark_select character; name of mark that is to be selected for the
#' reconstruction
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction.
#' Default = 500
#' @param bndw numeric; bandwith of the sigma parameter in the density estimation,
#' if no value is given the bandwith is estimated using cross validation with
#' the `bw.diggle` function.
#' @param thres numeric; intensity threshold for the reconstruction
#' @param ncores numeric; number of cores for parallel processing using
#' `mclapply`. Default = 1
#'
#' @importFrom parallel mclapply
#' @importFrom SummarizedExperiment assays
#'
#' @return simple feature collection
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' spe_sel <- spe[, spe[["image_name"]] %in% c("E02", "E03", "E04")]
#' all_islets <- reconstructShapeDensitySPE(spe_sel,
#'     marks = "cell_category",
#'     image_col = "image_name", mark_select = "islet", bndw = sigma, thres = 0.0025
#' )
#' all_islets
reconstructShapeDensitySPE <- function(
        spe, marks,
        image_col, mark_select,
        dim = 500, bndw = NULL, thres,
        ncores = 1) {
    # For computational reasonos delete all assays in SPE
    SummarizedExperiment::assays(spe) <- list()
    # Get all unique image ids
    all_images <- spe[[image_col]] |> unique()
    # Calculate polygon for each id using multiple cores
    res_all <- mclapply(all_images, function(x) {
        res <- reconstructShapeDensityImage(spe, marks, image_col,
            x, mark_select,
            dim = 500, bndw = NULL,
            thres
        )
        # assign image_id
        res[["structID"]] <- paste0(x, "_", c(1:dim(res)[1]))
        res[[image_col]] <- x
        return(res)
    }, mc.cores = ncores)
    # return data frame with all structures
    return(do.call(rbind, res_all))
}

#' Estimate reconstruction parameters from a set of images
#'
#' @param spe SpatialExperiment; a object of class `SpatialExperiment`
#' @param marks character; name of column in `colData` that will correspond to
#' the `ppp` marks
#' @param image_col character; name of a column in `colData` that corresponds
#' to the image
#' @param mark_select character; name of mark that is to be selected for the
#' reconstruction
#' @param nimages integer; number of images for the estimation. Will be randomly
#' sampled
#' @param fun character; function to estimate the kernel density. Default
#' bw.diggle.
#' @param dim numeric; x dimension of the final reconstruction.
#' A lower resolution speed up computation but lead to less exact reconstruction.
#' Default = 500
#' @param ncores numeric; number of cores for parallel processing using `mclapply`.
#' Default = 1
#' @param plot_hist logical; if histogram of estimated densities and thresholds
#' should be plotted. Default = TRUE
#'
#' @importFrom spatstat.geom subset.ppp
#' @importFrom SummarizedExperiment colData
#' @importFrom parallel mclapply
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggplot aes_string geom_histogram theme_light
#' @importFrom rlang .data
#'
#' @return tibble; tibble with estimated intensities and thresholds
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' spe_sel <- spe[, spe[["image_name"]] %in% c("E02", "E03", "E04")]
#' estimateReconstructionParametersSPE(spe_sel,
#'     marks = "cell_category",
#'     image_col = "image_name", mark_select = "islet", plot_hist = TRUE
#' )
estimateReconstructionParametersSPE <- function(spe,
    marks,
    image_col,
    mark_select = NULL,
    nimages = NULL,
    fun = "bw.diggle",
    dim = 500,
    ncores = 1,
    plot_hist = TRUE) {
    # Input checks
    if (!is.null(nimages)) {
        stopifnot("'nimages' must be numeric" = is.numeric(nimages))
        stopifnot(
            "'nimages' must be smaller or equal to the number of images in the `SpatialExperiment`" =
                (nimages < length(unique(colData(spe)[[image_col]])))
        )
    }

    # get the id's of all images
    all_images <- colData(spe)[[image_col]] |> unique()
    # default is to take all values
    if (is.null(nimages)) nimages <- length(all_images)
    # alternatively we sample some images
    sample_images <- sample(all_images, nimages)
    # we calculate the bandwidths and thresholds
    res <- mclapply(sample_images, function(x) {
        ppp <- SPE2ppp(spe, marks = marks, image_col = image_col, image_id = x)
        res_x <- .intensityImage(ppp, mark_select, dim = dim)
        thres <- .intensityThreshold(res_x$den_im)
        return(list(img = x, bndw = as.numeric(res_x$bndw), thres = as.numeric(thres)))
    }, mc.cores = ncores)

    # collect in one data frame
    res <- as.data.frame(do.call(rbind, res))
    res$bndw <- as.numeric(res$bndw)
    res$thres <- as.numeric(res$thres)

    if (plot_hist == TRUE & nimages > 1) {
        p1 <- res |>
            ggplot(aes(x = .data$bndw)) +
            geom_histogram(bins = round(nimages / 2)) +
            theme_light()

        p2 <- res |>
            ggplot(aes(x = .data$thres)) +
            geom_histogram(bins = round(nimages / 2)) +
            theme_light()

        plot(wrap_plots(p1, p2, ncol = 2))
    }

    return(res)
}
