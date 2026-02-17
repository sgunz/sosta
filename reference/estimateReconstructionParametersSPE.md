# Estimate reconstruction parameters from a set of images

Estimate reconstruction parameters from a set of images

## Usage

``` r
estimateReconstructionParametersSPE(
  spe,
  marks,
  imageCol = NULL,
  markSelect = NULL,
  nImages = NULL,
  fun = "bw.diggle",
  dim = 500,
  nCores = 1,
  plotHist = TRUE
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

- nImages:

  integer; number of images for the estimation. Will be randomly sampled

- fun:

  character; function to estimate the kernel density. Default bw.diggle.

- dim:

  numeric; x dimension of the final reconstruction. A lower resolution
  speed up computation but lead to less exact reconstruction. Default =
  500

- nCores:

  numeric; number of cores for parallel processing using `mclapply`.
  Default = 1

- plotHist:

  logical; if histogram of estimated densities and thresholds should be
  plotted (only if `imageCol` is not NULL). Default = TRUE

## Value

tibble; tibble with estimated intensities and thresholds

## Examples

``` r
data("sostaSPE")
estimateReconstructionParametersSPE(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    markSelect = "A", plotHist = TRUE
)

#>                 img     bndw      thres
#> image1 c(113.96.... 3.655770 0.04424890
#> image2 c(84.159.... 3.568156 0.05912054
#> image3 c(29.747.... 4.189788 0.03065195
```
