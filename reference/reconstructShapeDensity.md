# Reconstruct polygon from point pattern density

This function estimates the density of a spatial point pattern (`ppp`),
thresholds the density to create a binary image, and then converts it to
a valid `sf` object (polygons).

## Usage

``` r
reconstructShapeDensity(
  ppp,
  markSelect = NULL,
  bndw = NULL,
  thres = NULL,
  complement = FALSE,
  dim
)
```

## Arguments

- ppp:

  point pattern object of class `ppp`

- markSelect:

  character; name of mark that is to be selected for the reconstruction

- bndw:

  bandwidth of kernel density estimator

- thres:

  intensity threshold for the reconstruction

- complement:

  logical; reconstruct everything but the mark of interest, default =
  `FALSE`

- dim:

  numeric; x dimension of the final reconstruction.

## Value

sf object of class `POLYGON`

## Examples

``` r
data("sostaSPE")
ppp <- SPE2ppp(sostaSPE, marks = "cellType", imageCol = "imageName", imageId = "image1")
thres <- findIntensityThreshold(ppp, markSelect = "A", dim = 500)
struct <- reconstructShapeDensity(ppp, markSelect = "A", thres = thres, dim = 500)
plot(struct)
```
