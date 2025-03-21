set.seed(42)

test_that("binaryImageToSF works as expected", {
    binaryMatrix <- matrix(c(
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 1, 1, 1, 1, 0, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0,
        0, 1, 1, 1, 1, 1, 0, 0, 0,
        0, 1, 1, 0, 1, 1, 0, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0,
        0, 1, 1, 0, 0, 1, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0
    ), nrow = 9, byrow = TRUE)

    xmin <- 0
    xmax <- 1
    ymin <- 0
    ymax <- 1
    poly_sf <- binaryImageToSF(binaryMatrix, xmin, xmax, ymin, ymax)
    expect_s3_class(poly_sf, "sfc")
})

test_that("xyCoordinates returns correct coordinates", {
    inputMatrix <- matrix(c(
        0, 0, 0, 0, 0,
        0, 1, 0, 1, 0,
        0, 0, 0, 1, 1,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0
    ), nrow = 5, byrow = TRUE)

    coords <- xyCoordinates(inputMatrix)
    expect_equal(ncol(coords), 2)
    expect_true(all(coords >= 1))
    expect_true(all(coords <= 5))
})

test_that("normalizeCoordinates scales values correctly", {
    coords <- matrix(c(
        1, 1,
        1, 3,
        4, 1,
        4, 3
    ), ncol = 2, byrow = TRUE)

    normalized_coords <- normalizeCoordinates(coords)
    expect_equal(min(normalized_coords), 0)
    expect_equal(max(normalized_coords), 1)
    expect_equal(ncol(normalized_coords), 2)
})

test_that("getDimXY calculates dimensions based on y-axis", {
    ppp <- spatstat.geom::ppp(
        x = runif(100, 0, 1), y = runif(100, 0, 1),
        window = spatstat.geom::owin(c(0, 2), c(0, 2))
    )
    ydim <- 500
    dimyx <- getDimXY(ppp, ydim)
    expect_equal(length(dimyx), 2)
    expect_equal(dimyx[1], ydim)
})

test_that("findIntensityThreshold calculates threshold correctly", {
    ppp <- spatstat.geom::ppp(
        x = runif(100, 0, 1), y = runif(100, 0, 1),
        window = spatstat.geom::owin(c(0, 2), c(0, 2))
    )
    threshold <- findIntensityThreshold(ppp, markSelect = NULL, dim = 100)
    expect_true(is.numeric(threshold))
})

test_that(".intensityImage generates density image correctly", {
    # Create a sample ppp object
    ppp <- spatstat.geom::ppp(
        x = runif(100, 0, 1), y = runif(100, 0, 1),
        window = spatstat.geom::owin(c(0, 1), c(0, 1))
    )

    # Test with and without a specified bandwidth
    dim <- 100
    resultWithBandwidth <- .intensityImage(ppp, bndw = 0.1, dim = dim)
    resultWithoutBandwidth <- .intensityImage(ppp, dim = dim)

    # Ensure result is a list with expected elements
    expect_type(resultWithBandwidth, "list")
    expect_named(resultWithBandwidth, c("denIm", "bndw", "dimyx"))

    expect_type(resultWithoutBandwidth, "list")
    expect_named(resultWithoutBandwidth, c("denIm", "bndw", "dimyx"))

    # Check that 'denIm' is of correct class and contains density values
    expect_s3_class(resultWithBandwidth$denIm, "im")
    expect_true(all(resultWithBandwidth$denIm$v >= 0))

    # Validate that dimensions are calculated based on dim input
    expect_equal(resultWithBandwidth$dimyx[1], dim)
})

test_that(".intensityThreshold calculates a threshold based on density image", {
    # Mock a simple density image object using spatstat
    ppp <- spatstat.geom::ppp(
        x = runif(100, 0, 1), y = runif(100, 0, 1),
        window = spatstat.geom::owin(c(0, 1), c(0, 1))
    )

    densityImage <- density.ppp(ppp)

    # Calculate threshold
    threshold <- .intensityThreshold(densityImage, steps = 100)

    # Check threshold is numeric and within density range
    expect_true(is.numeric(threshold))
    expect_true(threshold >= min(densityImage$v))
    expect_true(threshold <= max(densityImage$v))
})



# TODO:  assingCellsToStructures
