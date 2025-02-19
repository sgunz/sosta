#' Function to assign spatial points to structures
#'
#' This function assigns each spatial point in a `SpatialExperiment` object (`spe`) to the first intersecting structure from a given set of spatial structures.
#'
#' @param spe SpatialExperiment; An object of class `SpatialExperiment` containing spatial point data.
#' @param all_structs sf; A simple feature collection (sf object) representing spatial structures. Must contain a column which contains a unique identifier for each structure. Default = `structID`.
#' @param image_col character; The column name in `spe` and `all_structs` that identifies the corresponding image.
#' @param n_cores integer; The number of cores to use for parallel processing (default is 1).
#'
#' @returns A vector with structure assignments for each spatial point in `spe`. Points that do not overlap with any structure are assigned `NA`.
#'
#' @importFrom sf st_intersects
#' @importFrom SummarizedExperiment assays
#' @importFrom parallel mclapply
#'
#' @examples
#' data(sostaSPE)
#' all_structures <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cell_type", image_col = "image_name",
#'     mark_select = "A", bndw = 3.5, thres = 0.045)
#' colData(sostaSPE)$struct_assign <- assingCellsToStructures(sostaSPE,
#'     all_structures, "image_name")
#' ggspavis::plotSpots(sostaSPE[, sostaSPE[["image_name"]] == "image1"],
#'     annotate = "struct_assign", sample_id = "sample_id",
#'     in_tissue = NULL, y_reverse = FALSE) + facet_wrap(~image_name)
#'
#' @export
assingCellsToStructures <- function(spe, all_structs, image_col, unique_id = "structID", n_cores = 1) {
    # Input checking
    stopifnot(
        "'spe' must be an object of class 'SpatialExperiment'" =
            inherits(spe, "SpatialExperiment")
    )
    stopifnot(
        "'all_structs' must be an object of class 'sf'" =
            inherits(all_structs, "sf")
    )
    stopifnot(
        "'image_col' must exist in colnames(all_structs)" =
            image_col %in% colnames(all_structs)
    )
    stopifnot(
        "'unique_id' must exist in colnames(all_structs)" =
            unique_id %in% colnames(all_structs)
    )
    # Extract unique image identifiers
    all_images <- unique(all_structs[[image_col]])
    # In order no to create memory problems we remove non relevant SPE entries
    SummarizedExperiment::assays(spe) <- list()

    # Using lapply to process each image separately
    res <- mclapply(all_images, function(sel) {
        # Create results vector with NA values
        res_vect <- rep(NA, ncol(spe))

        # Subset SPE object for the current image
        spe_sel <- spe[, spe[[image_col]] == sel]

        # Subset structure object for the current image
        structs_sel <- all_structs[all_structs[[image_col]] == sel, ]

        # Convert spatial coordinates to sf points object
        spatial_coords_sf <- spatialCoords2SF(spe_sel)

        # Compute intersections between spatial points and structures
        n <- sf::st_intersects(spatial_coords_sf, structs_sel, sparse = FALSE)

        # Extract the first structure ID for each point (if multiple, take the first)
        n_list <- apply(n, 1, function(x) which(x == TRUE)[1])

        # Assign structure ID or NA if no intersection
        res <- ifelse(n_list == 0, NA, structs_sel[[unique_id]][n_list])

        # Store results in the vector
        res_vect[spe[[image_col]] == sel] <- res
        return(res_vect)
    }, mc.cores = n_cores)

    # Create a data frame from the results
    df <- data.frame(do.call(cbind, res))

    # Extract the first non-NA value per row
    overlap_vect <- apply(df, 1, function(x) x[which(!is.na(x))[1]])

    return(overlap_vect)
}

#' Calculate the proportion of each cell type within spatial structures
#'
#' @param spe SpatialExperiment object
#' @param struct_column character; name of the `colData` column specifying the structure assignments
#' @param cell_type_column character; name of the `colData` column specifying cell types
#'
#' @return A data frame where rows correspond to unique structures and columns correspond to cell types,
#' containing the proportion of each cell type within each structure.
#'
#' @importFrom SingleCellExperiment colData
#'
#' @export
#'
#' @examples
#' data(sostaSPE)
#' all_structures <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cell_type", image_col = "image_name",
#'     mark_select = "A", bndw = 3.5, thres = 0.045)
#' colData(sostaSPE)$struct_assign <- assingCellsToStructures(sostaSPE,
#'     all_structures, "image_name")
#' cellTypeProportions(sostaSPE, "struct_assign", "cell_type")
cellTypeProportions <- function(spe, struct_column, cell_type_column) {
    # Extract structure assignments from column in SPE
    structs <- unique(spe[[struct_column]])
    # Remove NA values
    structs <- structs[!is.na(structs)]
    # Unique cell types from the specified column
    all_types <- unique(spe[[cell_type_column]])
    # Compute the proportion of each cell type within each structure
    res <- lapply(structs, function(sel) {
        sub_df <- colData(spe[, spe[[struct_column]] %in% sel])
        # Compute the frequency and normalize
        return(table(factor(sub_df[[cell_type_column]], levels = all_types)) /
                   length(sub_df[[cell_type_column]]))
    })
    # Combine into single df and name with structs id
    res_mat <- do.call(rbind, res) |> as.data.frame()
    rownames(res_mat) <- structs
    return(res_mat)
}


