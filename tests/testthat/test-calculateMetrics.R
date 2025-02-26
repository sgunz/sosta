# load data for tests
data(sostaSPE)
allStructs <- reconstructShapeDensitySPE(sostaSPE,
     marks = "cellType", imageCol = "imageName",
     markSelect = "A", bndw = 3.5, thres = 0.005
)
allStructs
metricsMatrix <- totalShapeMetrics(allStructs)


test_that("shapeMetrics computes correct shape metrics", {
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
    metrics <- shapeMetrics(polyR)

    expect_type(metrics, "list")
    expect_named(metrics, c(
        "Area", "Compactness", "Eccentricity", "Circularity",
        "Solidity", "Curl", "fibreLength", "fibreWidth"
    ))
    expect_gt(metrics$Area, 0)
    expect_gt(metrics$Compactness, 0)

    expect_gte(metrics$Eccentricity, 0)
    expect_lte(metrics$Eccentricity, 1)

    expect_gte(metrics$Circularity, 0)
    expect_lte(metrics$Circularity, 1)

    expect_gte(metrics$Solidity, 0)
    expect_lte(metrics$Solidity, 1)

    expect_gt(metrics$Curl, 0)
})

test_that("totalShapeMetrics computes shape metrics matrix correctly", {
    # Test function
    expect_type(metricsMatrix, "double")
    expect_true(ncol(metricsMatrix) > 0)
    expect_true(nrow(metricsMatrix) >= 1)
})

test_that("meanShapeMetrics computes mean shape metrics correctly", {
    # Test function
    meanMetrics <- meanShapeMetrics(metricsMatrix)
    expect_true(ncol(meanMetrics) == 1)
})

### Invalid inputs

test_that("shapeMetrics throws an error for invalid input", {
    invalidInput <- list("not a polygon")

    expect_error(shapeMetrics(invalidInput), "'sfPoly' must be a valid sfc object")
})

test_that("totalShapeMetrics throws an error for invalid input", {
    invalidInput <- list("not a multipolygon")

    expect_error(totalShapeMetrics(invalidInput), "'sfInput' must be a valid sf object")
})

test_that("meanShapeMetrics throws an error for invalid input", {
    invalidInput <- list("not a matrix")

    expect_error(meanShapeMetrics(invalidInput), "'totalShapeMetricMatrix' must be a matrix")
})
