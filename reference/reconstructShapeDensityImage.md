# Reconstruct structure from spe object with given image id

Reconstruct structure from spe object with given image id

## Usage

``` r
reconstructShapeDensityImage(
  spe,
  marks,
  imageCol = NULL,
  imageId = NULL,
  markSelect,
  dim = 500,
  bndw = NULL,
  thres = NULL,
  complement = FALSE
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

- imageId:

  character; image id, must be present in imageCol

- markSelect:

  character; name of mark that is to be selected for the reconstruction

- dim:

  numeric; x dimension of the final reconstruction. A lower resolution
  speed up computation but lead to less exact reconstruction. Default =
  500

- bndw:

  numeric; smoothing bandwidth in the density estimation, corresponds to
  the `sigma` parameter in the `density.ppp` function, if no value is
  given the bandwidth is estimated using cross validation with the
  `bw.diggle` function.

- thres:

  numeric; intensity threshold for the reconstruction; if NULL the
  threshold is set as the mean between the mode of the pixel intensity
  distributions

- complement:

  logical; reconstruct everything but the mark of interest, default =
  `FALSE`

## Value

sf object of class `POLYGON`

## Examples

``` r
data("sostaSPE")
struct <- reconstructShapeDensityImage(sostaSPE,
    marks = "cellType", imageCol = "imageName", imageId = "image1",
    markSelect = "A", dim = 500
)
plot(struct)
```
