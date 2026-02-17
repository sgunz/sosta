# Function to convert spatial coordinates of a `SpatialExperiment` object to a `ppp` object

Function to convert spatial coordinates of a `SpatialExperiment` object
to a `ppp` object

## Usage

``` r
SPE2ppp(spe, marks = NULL, imageCol = NULL, imageId = NULL)
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

## Value

ppp; object of type `ppp`

## Examples

``` r
data(sostaSPE)
SPE2ppp(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    imageId = "image1"
)
#> Marked planar point pattern: 1845 points
#> Multitype, with levels = A, B, C 
#> window: rectangle = [0.10611, 127.83934] x [0.11046, 127.96008] units
```
