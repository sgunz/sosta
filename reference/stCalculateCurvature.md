# Calculate curvature of sf object

Calculate curvature of sf object

## Usage

``` r
stCalculateCurvature(sfPoly, smoothness = 5)
```

## Arguments

- sfPoly:

  `POLYGON ` of class `sf`

- smoothness:

  list; curvature measures

## Value

list; list of curvatures values

## References

https://stackoverflow.com/questions/62250151/calculate-curvature-of-a-closed-object-in-r

## Examples

``` r
matrixR <- matrix(c(
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
polyR <- binaryImageToSF(matrixR, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
stCalculateCurvature(polyR)
#> $meanAbsCurv
#> [1] 5.3
#> 
#> $totalAbsCurv
#> [1] 116.6
#> 
#> $minCurv
#> [1] -9
#> 
#> $maxCurv
#> [1] 9
#> 
```
