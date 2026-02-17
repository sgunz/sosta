# Calculate the proportion of each cell type within spatial structures

Calculate the proportion of each cell type within spatial structures

## Usage

``` r
cellTypeProportions(spe, structColumn, cellTypeColumn, nCores = 1)
```

## Arguments

- spe:

  SpatialExperiment object

- structColumn:

  character; name of the `colData` column specifying the structure
  assignments

- cellTypeColumn:

  character; name of the `colData` column specifying cell types

- nCores:

  integer; The number of cores to use for parallel processing (default
  is 1).

## Value

A data frame where rows correspond to unique structures and columns
correspond to cell types, containing the proportion of each cell type
within each structure.

## Examples

``` r
library("SpatialExperiment")
data("sostaSPE")
allStructs <- reconstructShapeDensitySPE(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    markSelect = "A", bndw = 3.5, thres = 0.045
)
# The function `assingCellsToStructures` needs colnames so we create them here
colnames(sostaSPE) <- paste0("cell_", c(1:dim(sostaSPE)[2]))
# Assign the structure assignment in the order of the columns in the `SpatialExperiment` object
colData(sostaSPE)$structAssign <- assingCellsToStructures(
    spe = sostaSPE, allStructs = allStructs, imageCol = "imageName"
)[colnames(sostaSPE)]
cellTypeProportions(sostaSPE, "structAssign", "cellType")
#>                  A          B          C
#> image1_1 0.8166667 0.08333333 0.10000000
#> image1_2 0.8416076 0.07801418 0.08037825
#> image1_3 0.8076923 0.13286713 0.05944056
#> image1_4 0.8695652 0.08695652 0.04347826
#> image2_1 1.0000000 0.00000000 0.00000000
#> image2_2 0.8689408 0.08617594 0.04488330
#> image2_3 0.8812709 0.07023411 0.04849498
#> image3_1 0.8333333 0.16666667 0.00000000
#> image3_2 0.6595745 0.23404255 0.10638298
#> image3_3 1.0000000 0.00000000 0.00000000
#> image3_4 0.7894737 0.21052632 0.00000000
#> image3_5 0.8333333 0.11111111 0.05555556
#> image3_6 0.7407407 0.18518519 0.07407407
#> image3_7 0.8156124 0.09555855 0.08882907
```
