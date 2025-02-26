#' Function to assign spatial points to structures
#'
#' This function assigns each spatial point in a `SpatialExperiment` object (`spe`) to the first intersecting structure from a given set of spatial structures.
#'
#' @param spe SpatialExperiment; An object of class `SpatialExperiment` containing spatial point data.
#' @param allStructs sf; A simple feature collection (sf object) representing spatial structures. Must contain a column which contains a unique identifier for each structure. Default = `structID`.
#' @param imageCol character; The column name in `spe` and `allStructs` that identifies the corresponding image.
#' @param uniqueId character; The column name in the simple feature collection for which to compute the assignment.
#' @param nCores integer; The number of cores to use for parallel processing (default is 1).
#'
#' @returns A vector with structure assignments for each spatial point in `spe`. Points that do not overlap with any structure are assigned `NA`.
#'
#' @importFrom sf st_intersects
#' @importFrom SummarizedExperiment assays
#' @importFrom parallel mclapply
#'
#' @examples
#' library(ggspavis)
#' library(SpatialExperiment)
#' data(sostaSPE)
#' allStructs <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     markSelect = "A", bndw = 3.5, thres = 0.045)
#' colData(sostaSPE)$structAssign <- assingCellsToStructures(sostaSPE,
#'     allStructs, "imageName")
#' plotSpots(sostaSPE[, sostaSPE[["imageName"]] == "image1"],
#'     annotate = "structAssign", sample_id = "sample_id",
#'     in_tissue = NULL, y_reverse = FALSE) + facet_wrap(~imageName)
#'
#' @export
assingCellsToStructures <- function(spe, allStructs, imageCol, uniqueId = "structID", nCores = 1) {
    # Input checking
    stopifnot(
        "'spe' must be an object of class 'SpatialExperiment'" =
            inherits(spe, "SpatialExperiment")
    )
    stopifnot(
        "'allStructs' must be an object of class 'sf'" =
            inherits(allStructs, "sf")
    )
    stopifnot(
        "'imageCol' must exist in colnames(allStructs)" =
            imageCol %in% colnames(allStructs)
    )
    stopifnot(
        "'uniqueId' must exist in colnames(allStructs)" =
            uniqueId %in% colnames(allStructs)
    )
    # Extract unique image identifiers
    allImages <- unique(allStructs[[imageCol]])
    # In order no to create memory problems we remove non relevant SPE entries
    SummarizedExperiment::assays(spe) <- list()

    # Using lapply to process each image separately
    res <- mclapply(allImages, function(sel) {
        # Create results vector with NA values
        resVect <- rep(NA, ncol(spe))

        # Subset SPE object for the current image
        speSel <- spe[, spe[[imageCol]] == sel]

        # Subset structure object for the current image
        structsSel <- allStructs[allStructs[[imageCol]] == sel, ]

        # Convert spatial coordinates to sf points object
        spatialCoordsSf <- spatialCoords2SF(speSel)

        # Compute intersections between spatial points and structures
        n <- st_intersects(spatialCoordsSf, structsSel, sparse = FALSE)

        # Extract the first structure ID for each point (if multiple, take the first)
        n_list <- apply(n, 1, function(x) which(x == TRUE)[1])

        # Assign structure ID or NA if no intersection
        res <- ifelse(n_list == 0, NA, structsSel[[uniqueId]][n_list])

        # Store results in the vector
        resVect[spe[[imageCol]] == sel] <- res
        return(resVect)
    }, mc.cores = nCores)

    # Create a data frame from the results
    df <- data.frame(do.call(cbind, res))

    # Extract the first non-NA value per row
    overlapVect <- apply(df, 1, function(x) x[which(!is.na(x))[1]])

    return(overlapVect)
}

#' Calculate the proportion of each cell type within spatial structures
#'
#' @param spe SpatialExperiment object
#' @param structColumn character; name of the `colData` column specifying the structure assignments
#' @param cellTypeColumn character; name of the `colData` column specifying cell types
#'
#' @return A data frame where rows correspond to unique structures and columns correspond to cell types,
#' containing the proportion of each cell type within each structure.
#'
#' @importFrom SingleCellExperiment colData
#'
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' data(sostaSPE)
#' allStructs <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     markSelect = "A", bndw = 3.5, thres = 0.045)
#' colData(sostaSPE)$structAssign <- assingCellsToStructures(sostaSPE,
#'     allStructs, "imageName")
#' cellTypeProportions(sostaSPE, "structAssign", "cellType")
cellTypeProportions <- function(spe, structColumn, cellTypeColumn) {
    # Extract structure assignments from column in SPE
    structs <- unique(spe[[structColumn]])
    # Remove NA values
    structs <- structs[!is.na(structs)]
    # Unique cell types from the specified column
    allTypes <- unique(spe[[cellTypeColumn]])
    # Compute the proportion of each cell type within each structure
    res <- lapply(structs, function(sel) {
        sub_df <- colData(spe[, spe[[structColumn]] %in% sel])
        # Compute the frequency and normalize
        return(table(factor(sub_df[[cellTypeColumn]], levels = allTypes)) /
                   length(sub_df[[cellTypeColumn]]))
    })
    # Combine into single df and name with structs id
    res_mat <- do.call(rbind, res) |> as.data.frame()
    rownames(res_mat) <- structs
    return(res_mat)
}


#' Compute minimum boundary distances for each cell within its corresponding image structures
#'
#' @param spe SpatialExperiment object
#' @param imageColumn character; name of the `colData` column specifying the image name
#' @param structColumn character; name of the `colData` column specifying structure assignments
#' @param allStructs sf object; contains spatial structures with corresponding image names
#'
#' @return A numeric vector containing the minimum distances between cells and structure boundaries,
#' values within structures have negative values.
#'
#' @importFrom sf st_distance st_boundary
#'
#' @export
#'
#' @examples
#' library(ggspavis)
#' library(SpatialExperiment)
#' data(sostaSPE)
#' allStructs <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     markSelect = "A", bndw = 3.5, thres = 0.045)
#' colData(sostaSPE)$structAssign <- assingCellsToStructures(sostaSPE,
#'     allStructs, "imageName")
#' colData(sostaSPE)$minDist <- minBoundaryDistances(sostaSPE,
#'     "imageName", "structAssign", allStructs)
#' plotSpots(sostaSPE, annotate = "minDist", in_tissue = NULL, y_reverse = FALSE) +
#'     scale_colour_gradient2() +
#'     geom_sf(data = allStructs, fill = NA, inherit.aes = FALSE) +
#'     facet_wrap(~imageName)
minBoundaryDistances <- function(spe, imageColumn, structColumn, allStructs) {

    # Extract unique image names and remove NAs
    images <- unique(spe[[imageColumn]])
    images <- images[!is.na(images)]

    # Compute the minimum distance to structure boundaries for each cell
    res <- lapply(images, function(sel) {
        resVect <- rep(NA, ncol(spe))

        # Subset SPE and structures for the given image
        speSel <- spe[, spe[[imageColumn]] %in% sel]
        subStruct <- allStructs[allStructs[[imageColumn]] %in% sel, ]

        # If no structures exist for this image, return NA vector
        if (nrow(subStruct) == 0) {
            return(resVect)
        }

        # Compute distances between cell coordinates and structure boundaries
        dist <- sf::st_distance(spatialCoords2SF(speSel), sf::st_boundary(subStruct))

        # Store the minimum distance for each cell
        res <- apply(dist, 1, min)
        resVect[spe[[imageColumn]] %in% sel] <- res
        return(resVect)
    })

    # Combine results into a data frame
    df <- data.frame(do.call(cbind, res))

    # Extract the first non-NA value per row
    overlap_vect <- apply(df, 1, function(x) x[which(!is.na(x))[1]])

    # Negate distances for assigned structures
    overlap_vect[!is.na(spe[[structColumn]])] <- -overlap_vect[!is.na(spe[[structColumn]])]

    return(overlap_vect)
}



