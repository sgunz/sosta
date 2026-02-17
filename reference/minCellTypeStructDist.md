# Compute minimum distances from each cell types to structure boundaries per structure

Compute minimum distances from each cell types to structure boundaries
per structure

## Usage

``` r
minCellTypeStructDist(
  spe,
  allStructs,
  structID = "structID",
  cellTypeColumn,
  imageCol,
  nCores = 1
)
```

## Arguments

- spe:

  SpatialExperiment object

- allStructs:

  sf object; contains spatial structures with corresponding image names

- structID:

  character; name of the column in `allStructs` containing structure IDs
  (default: "structID")

- cellTypeColumn:

  character; name of the `colData` column specifying cell types

- imageCol:

  character; name of the `colData` column specifying the image name

- nCores:

  integer; Number of cores for parallel processing (default = 1)

## Value

A data frame where rows are structure IDs and cols are cell types with
minimum distances

## Examples

``` r
data("sostaSPE")
allStructs <- reconstructShapeDensitySPE(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    markSelect = "A", bndw = 3.5, thres = 0.045
)
minCellTypeStructDist(sostaSPE, allStructs, cellTypeColumn = "cellType", imageCol = "imageName")
#>                    A            B           C
#> image1_1 0.339579195 0.0177862558  0.18307652
#> image1_2 0.004753321 0.0057017800  0.09518254
#> image1_3 0.019610020 0.1209811147  1.98995331
#> image1_4 0.024085373 0.1227499022  0.32744306
#> image2_1 0.042182357 3.1526609387 20.17907379
#> image2_2 0.000000000 0.0136124758  0.04010389
#> image2_3 0.000000000 0.0008276229  0.19485275
#> image3_1 0.160624167 0.0018582197  1.23529566
#> image3_2 0.015645028 0.0000000000  0.72834897
#> image3_3 0.067773907 1.8206507161  1.18206828
#> image3_4 0.396058339 0.1933954248  5.97367009
#> image3_5 0.117146590 0.4274221810  0.29128786
#> image3_6 0.139339119 0.0085563439  5.41574547
#> image3_7 0.037328479 0.0125960808  0.07804082
#> image3_8 0.008340147 1.0793095450  1.12111461
```
