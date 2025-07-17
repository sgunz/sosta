#' Function to assign points / coordinates to structures
#'
#' This function assigns each spatial coordinate in a `SpatialExperiment` object (`spe`) to the first intersecting structure from a given set of spatial structures.
#'
#' @param spe SpatialExperiment; An object of class `SpatialExperiment` containing spatial point data. Must contain `colnames` for correct assignment.
#' @param allStructs sf; A simple feature collection (sf object) representing spatial structures. Must contain a column which contains a unique identifier for each structure. Default = `structID`.
#' @param imageCol character; The column name in `spe` and `allStructs` that identifies the corresponding image.
#' @param uniqueId character; The column name in the simple feature collection for which to compute the assignment.
#' @param nCores integer; The number of cores to use for parallel processing (default is 1).
#'
#' @returns A named list with structure assignments for each spatial point in `spe`. Points that do not overlap with any structure are assigned `NA`. Names correspond to `colnames` of the `SpatialExperiment` input object.
#'
#' @importFrom sf st_intersects
#' @importFrom SummarizedExperiment assays
#' @importFrom parallel mclapply
#'
#' @export
#'
#' @examples
#' library("SpatialExperiment")
#' data("sostaSPE")
#' allStructs <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     markSelect = "A", bndw = 3.5, thres = 0.045
#' )
#' # The function `assingCellsToStructures` needs colnames so we create them here
#' colnames(sostaSPE) <- paste0("cell_", c(1:dim(sostaSPE)[2]))
#'
#' res <- assingCellsToStructures(
#'     spe = sostaSPE, allStructs = allStructs, imageCol = "imageName"
#' )
#' # Assign the structure assignment in the order of the columns in the `SpatialExperiment` object
#' colData(sostaSPE)$structAssign <- res[colnames(sostaSPE)]
#'
#' if (require("ggplot2")) {
#'     cbind(
#'         colData(sostaSPE[, sostaSPE[["imageName"]] == "image1"]),
#'         spatialCoords(sostaSPE[, sostaSPE[["imageName"]] == "image1"])
#'     ) |>
#'         as.data.frame() |>
#'         ggplot(aes(x = x, y = y, color = structAssign)) +
#'         geom_point(size = 0.25) +
#'         coord_equal()
#' }
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
            length(imageCol) == 1 &&
                imageCol %in% colnames(allStructs)
    )
    stopifnot(
        "'uniqueId' must exist in colnames(allStructs)" =
            length(uniqueId) == 1 &&
                uniqueId %in% colnames(allStructs)
    )
    # Convert spe to df
    df <- .SPE2df(spe, imageCol, colNames = TRUE)
    # Split data frame
    ls <- split(df, as.factor(df[, imageCol]))

    # Using lapply to process each image separately
    res <- mclapply(ls, function(dfSel) {
        # Select image name
        sel <- unique(dfSel[, imageCol])

        # Subset structure object for the current image
        structsSel <- allStructs[allStructs[[imageCol]] == sel, ]

        # Convert spatial coordinates to sf points object
        spatialCoordsSf <- st_as_sf(dfSel[, c(1, 2)],
            coords = c(
                colnames(dfSel)[1],
                colnames(dfSel)[2]
            )
        )

        # Compute intersections between spatial points and structures
        n <- st_intersects(spatialCoordsSf, structsSel, sparse = FALSE)

        # Extract the first structure ID for each point (if multiple, take the first)
        n_list <- apply(n, 1, function(x) which(x == TRUE)[1])

        # Assign structure ID or NA if no intersection
        res <- ifelse(n_list == 0, NA, structsSel[[uniqueId]][n_list])
        return(data.frame(colnamesSPE = dfSel$colnamesSPE, structAssign = res))
    }, mc.cores = nCores)
    # bind results, extract assignemnt, name list
    rd <- do.call("rbind", res)
    out <- rd$structAssign
    names(out) <- rd$colnamesSPE
    return(out[colnames(spe)])
}

#' Calculate the proportion of each cell type within spatial structures
#'
#' @param spe SpatialExperiment object
#' @param structColumn character; name of the `colData` column specifying the structure assignments
#' @param cellTypeColumn character; name of the `colData` column specifying cell types
#' @param nCores integer; The number of cores to use for parallel processing (default is 1).
#'
#' @return A data frame where rows correspond to unique structures and columns correspond to cell types,
#' containing the proportion of each cell type within each structure.
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel mclapply
#' @importFrom S4Vectors split
#'
#' @export
#'
#' @examples
#' library("SpatialExperiment")
#' data("sostaSPE")
#' allStructs <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     markSelect = "A", bndw = 3.5, thres = 0.045
#' )
#' # The function `assingCellsToStructures` needs colnames so we create them here
#' colnames(sostaSPE) <- paste0("cell_", c(1:dim(sostaSPE)[2]))
#' # Assign the structure assignment in the order of the columns in the `SpatialExperiment` object
#' colData(sostaSPE)$structAssign <- assingCellsToStructures(
#'     spe = sostaSPE, allStructs = allStructs, imageCol = "imageName"
#' )[colnames(sostaSPE)]
#' cellTypeProportions(sostaSPE, "structAssign", "cellType")
cellTypeProportions <- function(spe, structColumn, cellTypeColumn, nCores = 1) {
    # Input checking
    stopifnot(
        "'spe' must be an object of class 'SpatialExperiment'" =
            inherits(spe, "SpatialExperiment")
    )
    stopifnot(
        "'structColumn' must exist in colData(allStructs)" =
            length(structColumn) == 1 &&
                structColumn %in% colnames(colData(spe))
    )
    stopifnot(
        "'cellTypeColumn' must exist in colData(cellTypeColumn)" =
            length(cellTypeColumn) == 1 &&
                cellTypeColumn %in% colnames(colData(spe))
    )
    # Unique cell types from the specified column
    allTypes <- unique(spe[[cellTypeColumn]])

    df <- colData(spe)[, c(structColumn, cellTypeColumn)]
    df <- df[!is.na(df[[structColumn]]), ]

    # Split data frame
    ls <- split(df, as.factor(df[, 1]))

    # Compute the proportion of each cell type within each structure
    res <- mclapply(ls, function(dfSub) {
        # Compute the frequency and normalize
        return(table(factor(dfSub[[cellTypeColumn]], levels = allTypes)) /
            length(dfSub[[cellTypeColumn]]))
    }, mc.cores = nCores)
    # Combine into single df and name with structs id
    resMat <- do.call(rbind, res) |> as.data.frame()
    return(resMat)
}


#' Compute minimum boundary distances for each cell within its corresponding image structures
#'
#' @param spe SpatialExperiment object
#' @param imageCol character; name of the `colData` column specifying the image name
#' @param structColumn character; name of the `colData` column specifying structure assignments
#' @param allStructs sf object; contains spatial structures with corresponding image names
#' @param nCores integer; The number of cores to use for parallel processing (default is 1).
#'
#' @return A named list containing the minimum distances between cells and structure boundaries,
#' values within structures have negative values. Names correspond to `colnames` of the `SpatialExperiment` input object.
#'
#' @importFrom sf st_distance st_boundary
#' @importFrom parallel mclapply
#'
#' @export
#'
#' @examples
#' library("SpatialExperiment")
#' data("sostaSPE")
#'
#' allStructs <- reconstructShapeDensitySPE(sostaSPE,
#'     marks = "cellType", imageCol = "imageName",
#'     markSelect = "A", bndw = 3.5, thres = 0.045
#' )
#' # The function `assingCellsToStructures` needs colnames so we create them here
#' colnames(sostaSPE) <- paste0("cell_", c(1:dim(sostaSPE)[2]))
#' # Assign the structure assignment in the order of the columns in the `SpatialExperiment` object
#' colData(sostaSPE)$structAssign <- assingCellsToStructures(
#'     spe = sostaSPE, allStructs = allStructs, imageCol = "imageName"
#' )[colnames(sostaSPE)]
#'
#' res <- minBoundaryDistances(
#'     spe = sostaSPE, imageCol = "imageName", structColumn = "structAssign",
#'     allStructs = allStructs
#' )
#'
#' colData(sostaSPE)$minDist <- res[colnames(sostaSPE)]
#'
#' if (require("ggplot2")) {
#'     cbind(colData(sostaSPE), spatialCoords(sostaSPE)) |>
#'         as.data.frame() |>
#'         ggplot(aes(x = x, y = y, color = minDist)) +
#'         geom_point(size = 0.25) +
#'         scale_colour_gradient2() +
#'         geom_sf(data = allStructs, fill = NA, inherit.aes = FALSE) +
#'         facet_wrap(~imageName)
#' }
minBoundaryDistances <- function(spe, imageCol,
    structColumn, allStructs, nCores = 1) {
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
        "'imageCol' must be a character string and exist in colData(spe) and colnames(allStructs)'" =
            is.character(imageCol) && length(imageCol) == 1 &&
                imageCol %in% colnames(colData(spe)) && imageCol %in% colnames(allStructs)
    )
    stopifnot(
        "'structColumn' must be a character string and exist in colData(spe)'" =
            is.character(structColumn) && length(structColumn) == 1 &&
                structColumn %in% colnames(colData(spe))
    )
    stopifnot(
        "'nCores' must be a positive integer'" =
            is.numeric(nCores) && length(nCores) == 1 &&
                nCores >= 1 && round(nCores) == nCores
    )

    # Extract unique image names and remove NAs
    images <- unique(spe[[imageCol]])
    images <- images[!is.na(images)]

    # Convert spe to df
    df <- .SPE2df(spe, imageCol, colNames = TRUE)
    # Split data frame
    ls <- split(df, as.factor(df[, imageCol]))

    # Compute the minimum distance to structure boundaries for each cell
    resdf <- mclapply(ls, function(dfSel) {
        # Select image name
        sel <- unique(dfSel[, imageCol])

        subStruct <- allStructs[allStructs[[imageCol]] %in% sel, ]

        # If no structures exist for this image, return NA vector
        if (nrow(subStruct) == 0) {
            return(rep(NA, nrow(dfSel)))
        }

        # Convert spatial coordinates to sf points object
        spatialCoordsSf <- st_as_sf(dfSel[, c(1, 2)],
            coords = c(
                colnames(dfSel)[1],
                colnames(dfSel)[2]
            )
        )

        # Compute distances between cell coordinates and structure boundaries
        dist <- sf::st_distance(spatialCoordsSf, sf::st_boundary(subStruct))

        # Store the minimum distance for each cell
        res <- apply(dist, 1, min)
        return(data.frame(colnamesSPE = dfSel$colnamesSPE, minDist = res))
    }, mc.cores = nCores)

    # bind results, extract assignment, name list
    rd <- do.call("rbind", resdf)
    out <- rd$minDist
    names(out) <- rd$colnamesSPE

    # Negate distances for assigned structures
    out[colnames(spe)][!is.na(spe[[structColumn]])] <-
        -out[colnames(spe)][!is.na(spe[[structColumn]])]

    return(out[colnames(spe)])
}

