# Calculate the length of feature axes of a single sf polygon

Calculate the length of feature axes of a single sf polygon

## Usage

``` r
stFeatureAxes(sfPoly)
```

## Arguments

- sfPoly:

  `POLYGON ` of class `sf`

## Value

list; list containing the major and minor axis lengths

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
stFeatureAxes(polyR)
#> $majorAxisLength
#> [1] 0.7777778
#> 
#> $minorAxisLength
#> [1] 0.6666667
#> 
```
