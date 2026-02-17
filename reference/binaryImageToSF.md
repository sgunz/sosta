# Converts a binary matrix to an sf polygon

Converts a binary matrix to an sf polygon

## Usage

``` r
binaryImageToSF(binaryMatrix, xmin, xmax, ymin, ymax)
```

## Arguments

- binaryMatrix:

  matrix; binary matrix

- xmin:

  integer; minimum x coordinate of the coordinate system

- xmax:

  integer; maximum x coordinate of the coordinate system

- ymin:

  integer; minimum y coordinate of the coordinate system

- ymax:

  integer; maximum y coordinate of the coordinate system

## Value

sf object

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
plot(polyR)
```
