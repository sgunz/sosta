# Function to estimate the intensity threshold for the reconstruction of spatial structures

Function to estimate the intensity threshold for the reconstruction of
spatial structures

## Usage

``` r
.intensityThreshold(densityImage, steps = 250)
```

## Arguments

- densityImage:

  real-valued pixel image; output from the function `.intensityImage`

- steps:

  numeric; value used to filter the density estimates, where only
  densities greater than the maximum value divided by `threshold` are
  considered. Default is 250.

## Value

numeric; estimated threshold
