library("SpatialExperiment")
data("sostaSPE")
allStructs <- reconstructShapeDensitySPE(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    markSelect = "A", bndw = 3.5, thres = 0.045
)
colData(sostaSPE)$structAssign <- assingCellsToStructures(
    sostaSPE,
    allStructs, "imageName"
)

# Test assingCellsToStructures
test_that("assingCellsToStructures returns correct output", {
    result <- assingCellsToStructures(sostaSPE, allStructs, "imageName")
    expect_type(result, "character")
    expect_length(result, ncol(sostaSPE))
})

test_that("assingCellsToStructures handles incorrect input", {
    expect_error(assingCellsToStructures(list(), allStructs, "imageName"), "must be an object of class 'SpatialExperiment'")
    expect_error(assingCellsToStructures(sostaSPE, list(), "imageName"), "must be an object of class 'sf'")
})

# Test cellTypeProportions
test_that("cellTypeProportions returns a data frame", {
    result <- cellTypeProportions(sostaSPE, "structAssign", "cellType")
    expect_s3_class(result, "data.frame")
    expect_true(all(colnames(result) %in% unique(sostaSPE$cellType)))
})

test_that("cellTypeProportions handles incorrect input", {
    expect_error(cellTypeProportions(list(), "structAssign", "cellType"), "must be an object of class 'SpatialExperiment'")
    expect_error(cellTypeProportions(sostaSPE, "wrongColumn", "cellType"), "must exist in colData")
})

# Test minBoundaryDistances
test_that("minBoundaryDistances returns correct numeric vector", {
    result <- minBoundaryDistances(sostaSPE, "imageName", "structAssign", allStructs)
    expect_type(result, "double")
    expect_length(result, ncol(sostaSPE))
})

test_that("minBoundaryDistances handles incorrect input", {
    expect_error(minBoundaryDistances(list(), "imageName", "structAssign", allStructs), "must be an object of class 'SpatialExperiment'")
    expect_error(minBoundaryDistances(sostaSPE, "wrongColumn", "structAssign", allStructs), "must be a character string and exist in colData")
})
