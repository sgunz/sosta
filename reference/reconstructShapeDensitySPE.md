# Reconstruct structure from spatial experiment object per image id

Reconstruct structure from spatial experiment object per image id

## Usage

``` r
reconstructShapeDensitySPE(
  spe,
  marks,
  imageCol = NULL,
  markSelect,
  dim = 500,
  bndw = NULL,
  thres = NULL,
  complement = FALSE,
  nCores = 1
)
```

## Arguments

- spe:

  SpatialExperiment; a object of class `SpatialExperiment`

- marks:

  character; name of column in `colData` that will correspond to the
  `ppp` marks

- imageCol:

  character; name of a column in `colData` that corresponds to the image

- markSelect:

  character; name of mark that is to be selected for the reconstruction

- dim:

  numeric; x dimension of the final reconstruction. A lower resolution
  speed up computation but lead to less exact reconstruction. Default =
  500

- bndw:

  numeric; bandwidth of the sigma parameter in the density estimation,
  if no value is given the bandwidth is estimated using cross validation
  with the `bw.diggle` function for each image individually.

- thres:

  numeric; intensity threshold for the reconstruction; if NULL the
  threshold is set as the mean between the mode of the pixel intensity
  distributions estimated for each image individual

- complement:

  logical; reconstruct everything but the mark of interest, default =
  `FALSE`

- nCores:

  numeric; number of cores for parallel processing using `mclapply`.
  Default = 1

## Value

simple feature collection

## Examples

``` r
data("sostaSPE")
allStructs <- reconstructShapeDensitySPE(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    markSelect = "A", bndw = 3.5, thres = 0.005
)
allStructs
#> Simple feature collection with 5 features and 2 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 0.003206789 ymin: 0.01876509 xmax: 127.9648 ymax: 127.9749
#> CRS:           NA
#>          structID imageName                   sostaPolygon
#> image1   image1_1    image1 POLYGON ((52.22127 127.9601...
#> image2.1 image2_1    image2 POLYGON ((0.5238102 13.0702...
#> image2.2 image2_2    image2 POLYGON ((12.29547 127.719,...
#> image2.3 image2_3    image2 POLYGON ((52.47265 127.4631...
#> image3   image3_1    image3 POLYGON ((25.3302 127.8568,...
```
