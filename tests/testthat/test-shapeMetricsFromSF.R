test_that("stFeatureAxes computes major and minor axis lengths correctly", {
    # Create a simple polygon
    matrixR <- matrix(c(
        0, 1, 1, 0,
        0, 1, 1, 0,
        0, 1, 1, 0,
        0, 1, 1, 0
    ), nrow = 4, byrow = TRUE)

    polyR <- binaryImageToSF(matrixR, xmin = 0, xmax = 1, ymin = 0, ymax = 1)

    # Test function
    axes <- stFeatureAxes(polyR)

    expect_type(axes, "list")
    expect_named(axes, c("majorAxisLength", "minorAxisLength"))
    expect_equal(axes$majorAxisLength, 1)
    expect_equal(axes$minorAxisLength, 0.5)
})

test_that("stCalculateCurvature computes curvature metrics correctly", {
    # Create a simple polygon
    matrixR <- matrix(c(
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

    polyR <- binaryImageToSF(matrixR, xmin = 0, xmax = 1, ymin = 0, ymax = 1)

    # Test function
    curvature <- stCalculateCurvature(polyR)

    expect_type(curvature, "list")
    expect_named(curvature, c("meanAbsCurv", "totalAbsCurv", "minCurv", "maxCurv"))
    expect_true(curvature$meanAbsCurv > 0)
    expect_true(curvature$totalAbsCurv > 0)
})

test_that("stCalculateShapeCurl computes curl correctly", {
    # Create a simple polygon
    matrixR <- matrix(c(
        1, 1, 1, 1, 1, 0,
        1, 1, 0, 0, 1, 1,
        1, 1, 0, 0, 1, 1,
        1, 1, 1, 1, 1, 0,
        1, 1, 0, 1, 1, 0,
        1, 1, 0, 0, 1, 1,
        1, 1, 0, 0, 1, 1
    ), nrow = 7, byrow = TRUE)

    polyR <- binaryImageToSF(matrixR, xmin = 0, xmax = 1, ymin = 0, ymax = 1)

    # Test function
    curl <- stCalculateShapeCurl(polyR)

    expect_type(curl, "list")
    expect_named(curl, c("Curl", "fibreLength", "fibreWidth"))
    expect_true(curl$Curl >= 0)
    expect_true(curl$fibreLength > 0)
    expect_true(curl$fibreWidth > 0)
})

### Invalid inputs

test_that("stFeatureAxes throws an error for invalid input", {
    invalidInput <- list("not a polygon")

    expect_error(stFeatureAxes(invalidInput), "'sfPoly' must be a valid sfc object")
})

test_that("stCalculateCurvature throws an error for invalid input", {
    invalidInput <- list("not a polygon")

    expect_error(stCalculateCurvature(invalidInput), "'sfPoly' must be a valid sfc object")
})

test_that("stCalculateShapeCurl throws an error for invalid input", {
    invalidInput <- list("not a polygon")

    expect_error(stCalculateShapeCurl(invalidInput), "'sfPoly' must be a valid sfc object")
})
