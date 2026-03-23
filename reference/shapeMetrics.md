# Calculate a set of shape metrics of a single polygon

Calculate a set of shape metrics of a single polygon

## Usage

``` r
shapeMetrics(sfPoly)
```

## Arguments

- sfPoly:

  POLYGON of class sfc

## Value

list; list of shape metrics

## Details

For multiple polyogns or a MULTIPOLYGON object use the function
[`totalShapeMetrics`](https://sgunz.github.io/sosta/reference/totalShapeMetrics.md)

## Examples

``` r
matrix_R <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 1, 1, 1, 0, 0, 0,
    0, 1, 1, 0, 0, 1, 1, 0, 0,
    0, 1, 1, 0, 0, 1, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 0, 0, 0,
    0, 1, 1, 0, 1, 1, 0, 0, 0,
    0, 1, 1, 0, 0, 1, 1, 0, 0,
    0, 1, 1, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0
), nrow = 9, byrow = TRUE)
polyR <- binaryImageToSF(matrix_R, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
shapeMetrics(polyR)
#> $Area
#> [1] 0.3703704
#> 
#> $Compactness
#> [1] 0.2137138
#> 
#> $Eccentricity
#> [1] 0.8571429
#> 
#> $Circularity
#> [1] 0.583684
#> 
#> $Solidity
#> [1] 0.7228916
#> 
#> $Curl
#> [1] 0.6402552
#> 
#> $fibreLength
#> [1] 2.162026
#> 
#> $fibreWidth
#> [1] 0.1713071
#> 
```
