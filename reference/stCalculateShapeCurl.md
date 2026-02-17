# Calculate curl of a polygon

Calculate curl of a polygon

## Usage

``` r
stCalculateShapeCurl(sfPoly)
```

## Arguments

- sfPoly:

  `POLYGON ` of class `sf`

## Value

numeric; the curl of the polygon

## Examples

``` r
matrixR <- matrix(c(
    1, 1, 1, 1, 1, 0,
    1, 1, 0, 0, 1, 1,
    1, 1, 0, 0, 1, 1,
    1, 1, 1, 1, 1, 0,
    1, 1, 0, 1, 1, 0,
    1, 1, 0, 0, 1, 1,
    1, 1, 0, 0, 1, 1
), nrow = 7, byrow = TRUE)
polyR <- binaryImageToSF(matrixR, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
stCalculateShapeCurl(polyR)
#> $Curl
#> [1] 0.6637659
#> 
#> $fibreLength
#> [1] 2.974119
#> 
#> $fibreWidth
#> [1] 0.2401672
#> 
```
