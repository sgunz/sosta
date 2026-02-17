# Estimate the intensity threshold for the reconstruction of spatial structures

Estimate the intensity threshold for the reconstruction of spatial
structures

## Usage

``` r
findIntensityThreshold(ppp, markSelect = NULL, bndw = NULL, dim, steps = 250)
```

## Arguments

- ppp:

  point pattern object of class `ppp`

- markSelect:

  character; name of mark that is to be selected for the reconstruction

- bndw:

  numeric; bandwith of the sigma parameter in the density estimation, if
  no value is given the bandwith is estimated using cross validation
  with the `bw.diggle` function.

- dim:

  numeric; x dimension of the final reconstruction.

- steps:

  numeric; value used to filter the density estimates, where only
  densities greater than the maximum value divided by `threshold` are
  considered. Default is 250.

## Value

numeric; estimated intensity threshold

## Examples

``` r
data(sostaSPE)
ppp <- SPE2ppp(sostaSPE, marks = "cellType", imageCol = "imageName", imageId = "image1")
findIntensityThreshold(ppp, markSelect = "A", dim = 250)
#> [1] 0.04595513
```
