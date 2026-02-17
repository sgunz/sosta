# Function to get the dimension based on dim of y axis

Function to get the dimension based on dim of y axis

## Usage

``` r
getDimXY(ppp, ydim)
```

## Arguments

- ppp:

  point pattern object of class `ppp`

- ydim:

  dimension of y axis

## Value

vector; vector with x and y dimension

## Examples

``` r
data(sostaSPE)
pp <- SPE2ppp(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    imageId = "image1"
)
getDimXY(pp, 500)
#> [1] 500 500
```
