# Create a dummy point pattern
ppp <- spatstat.geom::ppp(
    x = runif(100),
    y = runif(100),
    window = spatstat.geom::owin(c(0, 1), c(0, 1))
)

# load sostaSPE object
data("sostaSPE")

# Test the reconstruction function
polyA <- reconstructShapeDensityImage(
    sostaSPE,
    marks = "cellType",
    imageCol = "imageName",
    imageId = "image1",
    markSelect = "A",
    dim = 500
)

# Test the reconstruction on sostaSPE object
allA <- reconstructShapeDensitySPE(
    sostaSPE,
    marks = "cellType",
    imageCol = "imageName",
    markSelect = "A",
    bndw = 3.5,
    thres = 0.005
)

test_that("reconstructShapeDensity returns valid polygons", {
    # Reconstruct polygons with valid parameters
    result <- reconstructShapeDensity(ppp, markSelect = NULL, bndw = 1, dim = 100)

    expect_s3_class(result, "sf")
    expect_true(any(st_geometry_type(result) == "POLYGON"))
    expect_false(any(st_is_empty(result)))
})

test_that("reconstructShapeDensity handles invalid thresholds", {
    # Test low threshold
    expect_warning(
        reconstructShapeDensity(ppp, thres = 0, bndw = 1, dim = 500),
        "Full image converted to polygon; threshold might be too low"
    )

    # Test high threshold
    expect_warning(
        reconstructShapeDensity(ppp, thres = 1E5, bndw = 1, dim = 500),
        "No structure found; threshold might be too high"
    )
})

test_that("shapeIntensityImage generates valid plots", {
    # Test the intensity image function
    p <- shapeIntensityImage(
        sostaSPE,
        marks = "cellType",
        imageCol = "imageName",
        imageId = "image1",
        markSelect = "A"
    )

    expect_s3_class(p, "ggplot")
})

test_that("reconstructShapeDensityImage returns polygons from SpatialExperiment", {
    expect_s3_class(polyA, "sf")
    expect_true(any(st_geometry_type(polyA) == "POLYGON"))
    expect_false(any(st_is_empty(polyA)))
})

test_that("reconstructShapeDensitySPE handles multiple images", {
    expect_s3_class(allA, "sf")
    expect_true(any(st_geometry_type(allA) == "POLYGON"))
    expect_false(any(st_is_empty(allA)))
})

test_that("estimateReconstructionParametersSPE estimates valid parameters", {
    # Test the estimation function
    res <- estimateReconstructionParametersSPE(
        sostaSPE,
        marks = "cellType",
        imageCol = "imageName",
        markSelect = "A",
        nImages = 2,
        dim = 500,
        plotHist = FALSE
    )

    expect_true(is.data.frame(res))
    expect_named(res, c("img", "bndw", "thres"))
    expect_true(all(res$bndw > 0))
    expect_true(all(res$thres > 0))
})

test_that("reconstructShapeDensity handles invalid input types", {
    # Invalid input (non-ppp object)
    expect_error(
        reconstructShapeDensity("not_a_ppp"),
        "'ppp' must be an object of class 'ppp'"
    )

    expect_error(
        reconstructShapeDensity(ppp, dim = -500),
        "'dim' must be a single, positive, numeric value"
    )
})

test_that("estimateReconstructionParametersSPE handles edge cases", {
    # Test with fewer images than requested
    expect_error(
        estimateReconstructionParametersSPE(
            sostaSPE,
            marks = "cellType",
            imageCol = "imageName",
            markSelect = "A",
            nImages = 50,
            dim = 500,
            plotHist = FALSE
        )
    )
})


test_that("shapeIntensityImage returns output if imageId = NULL", {
    result <- shapeIntensityImage(
        spe = sostaSPE,
        marks = "cellType",
        imageCol = "imageName",
        imageId = NULL,
        markSelect = "A"
    )

    expect_false(is.null(result))
    expect_s3_class(result, "gg")
})

test_that("reconstructShapeDensityImage returns output if imageId = NULL", {
    result <- reconstructShapeDensityImage(
        spe = sostaSPE,
        marks = "cellType",
        imageCol = "imageName",
        imageId = NULL,
        markSelect = "A",
        dim = 50
    )

    expect_false(is.null(result))
    expect_s3_class(result, "sf")
})

test_that("reconstructShapeDensitySPE runs when imageCol is NULL", {
    result <- reconstructShapeDensitySPE(
        spe = sostaSPE,
        marks = "cellType",
        imageCol = NULL,
        markSelect = "A",
        dim = 50
    )

    expect_false(is.null(result))
    expect_s3_class(result, "sf")
})

test_that("estimateReconstructionParametersSPE works if imageCol = NULL", {
    result <- estimateReconstructionParametersSPE(
        spe = sostaSPE,
        marks = "cellType",
        imageCol = NULL,
        markSelect = "A",
        nImages = NULL,
        plotHist = FALSE
    )

    expect_false(is.null(result))
    expect_true(is.list(result) || is.data.frame(result))
})


