#' Reconstruct polygon from point pattern density
#'
#' This function estimates the density of a spatial point pattern (`ppp`),
#' thresholds the density to create a binary image, and then converts it
#' to a valid `sf` object (polygons).
#'
#' @param ppp point pattern object of class `ppp`
#' @param bndw bandwidth of kernel density estimator
#' @param thres intensity threshold for the reconstruction
#' @param dimyx pixel dimensions of the output image
#'
#' @return sf object of class `POLYGON`
#' @import spatstat.explore spatstat.geom sf
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' pp <- SPE2ppp(spe, marks = "cell_category", image_col = "image_name", image_id = "E04")
#' pp_sel <- subset(pp, marks == "islet")
#' bndw <- spatstat.explore::bw.diggle(pp_sel)
#' dimyx <- getDimXY(pp_sel, 500)
#' thres <- findIntensityThreshold(pp_sel, bndw = bndw, dimyx = dimyx)
#' islet_poly <- reconstructShapeDensity(pp_sel, bndw = bndw, thres = thres, dimyx = dimyx)
#' plot(islet_poly)
reconstructShapeDensity <- function(ppp, bndw, thres, dimyx) {
    # estimate density
    density_image <- density.ppp(ppp, bndw, dimyx = c(dimyx), positive = TRUE)

    # TODO: Error when threshold is too small
    # construct spatstat window from matrix with true false entries
    mat <- ifelse(t(as.matrix(density_image)) > thres, TRUE, FALSE)

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
    ) # make valid is important!!
    stCast <- stCast[!st_is_empty(stCast), drop = FALSE]

    # return sf object
    return(st_sf(st_cast(stCast, "POLYGON")))
}


#' Get intensity plot and intensity from spe object with given image id
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
#' labs coord_equal theme_classic scale_color_viridis_c
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom dplyr filter
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' shapeIntensityImage(spe,
#'     marks = "cell_category", image_col = "image_name",
#'     image_id = "E04", mark_select = "islet"
#' )
shapeIntensityImage <- function(
        spe, marks,
        image_col,
        image_id,
        mark_select,
        bndw = NULL,
        dim = 500) {
    # Convert the spe object to a point pattern object
    pp <- SPE2ppp(spe, marks = marks, image_col = image_col, image_id = image_id)

    # Extract the islet cells
    pp.islet <- subset(pp, marks %in% mark_select)

    # Set the dimensions of the resulting reconstruction
    dimyx <- getDimXY(pp.islet, dim)

    # Set default of sigma bandwith
    if (is.null(bndw)) bndw <- bw.diggle(pp.islet)

    # plot the density of the image
    im_df <- density.ppp(pp.islet,
        sigma = bndw,
        dimyx = dimyx,
        positive = TRUE
    ) |> as.data.frame()

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
        dplyr::filter(.data$value > max(.data$value) / 250) |>
        ggplot(aes(x = abs(.data$value))) + # Use .data pronoun
        geom_histogram(bins = 50) +
        labs(x = "pixel intensity") +
        theme_light()


    p <- patchwork::wrap_plots(den_im, den_hist, ncol = 2) +
        patchwork::plot_annotation(
            title = paste0(image_col, ": ", image_id),
            subtitle = paste0("bndw: ", round(bndw, 4)),
            caption = paste0("Pixel image dimensions: ", dimyx[1], "x", dimyx[2])
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
#'
#' @return sf object of class `POLYGON`
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' islet_poly <- reconstructShapeDensityImage(spe,
#'     marks = "cell_category",
#'     image_col = "image_name", image_id = "E04", mark_select = "islet", dim = 500
#' )
#' plot(islet_poly)
reconstructShapeDensityImage <- function(spe, marks,
    image_col, image_id, mark_select,
    dim = 500, bndw = NULL,
    thres = NULL) {
    # Convert the spe object to a point pattern object
    pp <- SPE2ppp(
        spe,
        marks,
        image_col,
        image_id
    )

    # Extract the islet cells
    pp_sel <- subset(pp, marks %in% mark_select)

    # Get dimension of reconstruction
    dimyx <- getDimXY(pp_sel, dim)

    # Set default of sigma bandwith
    if (is.null(bndw)) bndw <- bw.diggle(pp_sel)

    # Find the threshold for the structure if not set
    if (is.null(thres)) {
        thres <- findIntensityThreshold(pp_sel,
            bndw = bndw,
            dimyx = dimyx
        )
    }
    # Get the structure
    struct <- reconstructShapeDensity(pp_sel,
        bndw = bndw,
        thres = thres,
        dimyx = dimyx
    )

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
#' @importFrom dplyr bind_rows
#'
#' @return simple feature collection
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' spe <- spe[, spe[["image_name"]] %in% c("E02", "E03", "E04")]
#' all_islets <- reconstructShapeDensitySPE(spe,
#'     marks = "cell_category",
#'     image_col = "image_name", mark_select = "islet", bndw = sigma, thres = 0.0025
#' )
#' all_islets
reconstructShapeDensitySPE <- function(spe, marks,
    image_col, mark_select,
    dim = 500, bndw = NULL, thres,
    ncores = 1) {
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
        res[[image_col]] <- x
        return(res)
    }, mc.cores = ncores)
    # return data frame with all structures
    # TODO: do.call(rbind, a)
    return(bind_rows(res_all))
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
#' @importFrom parallel mclapply
#' @importFrom dplyr bind_rows
#' @import spatstat.explore
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggplot aes_string geom_histogram theme_light
#' @return tibble; tibble with estimated intensities and thresholds
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' spe <- spe[, spe[["image_name"]] %in% c("E02", "E03", "E04")]
#' estimateReconstructionParametersSPE(spe,
#'     marks = "cell_category",
#'     image_col = "image_name", mark_select = "islet", plot_hist = TRUE
#' )
estimateReconstructionParametersSPE <- function(
        spe, marks, image_col,
        mark_select = NULL, nimages = NULL, fun = "bw.diggle", dim = 500,
        ncores = 1, plot_hist = TRUE) {
    # get the id's of all images
    all_images <- colData(spe)[[image_col]] |> unique()
    # default is to take all values
    if (is.null(nimages)) nimages <- length(all_images)
    # alternitavely we sample some images
    sample_images <- sample(all_images, nimages)
    # we calculate the bandwiths and thresholds
    res <- parallel::mclapply(sample_images, function(x) {
        pp <- SPE2ppp(spe, marks = marks, image_col = image_col, image_id = x)
        # If selection of mark
        if (!is.null(mark_select)) pp <- subset(pp, marks %in% mark_select)
        # Estimate the bandwidth for the kernel estimation of point process intensity
        bndw <- do.call(fun, args = list(X = pp, warn = FALSE))
        # Estimate the threshold for the reconstruction
        # Get dimension of reconstruction
        dimyx <- getDimXY(pp, dim)
        thres <- findIntensityThreshold(pp, bndw = bndw, dimyx = dimyx)
        return(c("img" = x, "bndw" = as.numeric(bndw), "thres" = as.numeric(thres)))
    }, mc.cores = ncores)

    res <- dplyr::bind_rows(res)
    res$thres <- as.numeric(res$thres)
    res$bndw <- as.numeric(res$bndw)

    if (plot_hist == TRUE & nimages > 1) {
        p1 <- res |>
            ggplot(aes(x = .data$bndw)) +
            geom_histogram(bins = round(nimages / 2)) +
            theme_light()

        p2 <- res |>
            ggplot(aes(x = .data$thres)) +
            geom_histogram(bins = round(nimages / 2)) +
            theme_light()

        plot(patchwork::wrap_plots(p1, p2, ncol = 2))
    }

    return(res)
}
