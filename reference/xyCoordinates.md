# Function to extract x y coordinates from binary image

Function to extract x y coordinates from binary image

## Usage

``` r
xyCoordinates(inputMatrix)
```

## Arguments

- inputMatrix:

  a binary matrix

## Value

matrix; matrix with x,y coordinates of the cell of the input matrix

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
xyCoordinates(matrixR)
#>       x y
#>  [1,] 2 2
#>  [2,] 3 2
#>  [3,] 4 2
#>  [4,] 5 2
#>  [5,] 6 2
#>  [6,] 7 2
#>  [7,] 8 2
#>  [8,] 2 3
#>  [9,] 3 3
#> [10,] 4 3
#> [11,] 5 3
#> [12,] 6 3
#> [13,] 7 3
#> [14,] 8 3
#> [15,] 2 4
#> [16,] 5 4
#> [17,] 2 5
#> [18,] 5 5
#> [19,] 6 5
#> [20,] 2 6
#> [21,] 3 6
#> [22,] 4 6
#> [23,] 5 6
#> [24,] 6 6
#> [25,] 7 6
#> [26,] 8 6
#> [27,] 3 7
#> [28,] 4 7
#> [29,] 7 7
#> [30,] 8 7
```
