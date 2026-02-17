# Function to estimate the intensity image of a point pattern

Function to estimate the intensity image of a point pattern

## Usage

``` r
.intensityImage(ppp, markSelect = NULL, bndw = NULL, dim)
```

## Arguments

- ppp:

  point pattern object of class `ppp`

- markSelect:

  character; name of mark that is to be selected for the reconstruction

- bndw:

  bandwidth of kernel density estimator

- dim:

  numeric; x dimension of the final reconstruction.

## Value

list; list with the intensity image and the bandwidth and dimension
parameters
