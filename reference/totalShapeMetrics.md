# Calculate a set of shape metrics of a set of polygons

Calculate a set of shape metrics of a set of polygons

## Usage

``` r
totalShapeMetrics(sfInput)
```

## Arguments

- sfInput:

  `MULTIPOLYGON` of class sf

## Value

matrix; matrix of shape metrics

## Details

Calculate a set of shape metrics of a set of polygons. The function
calculates all metrics that are implemented in the function
[`shapeMetrics()`](https://sgunz.github.io/sosta/reference/shapeMetrics.md)

## Examples

``` r
data(sostaSPE)
struct <- reconstructShapeDensityImage(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    imageId = "image1", markSelect = "A", dim = 500
)
totalShapeMetrics(struct)
#>                   struct1      struct2     struct3
#> Area         4230.8784155 3560.7340306 228.1064618
#> Compactness     0.1310661    0.2426745   0.6019973
#> Eccentricity    0.7853205    0.4645111   0.8250754
#> Circularity     0.4707314    0.5949015   0.8843705
#> Solidity        0.5940693    0.8009580   0.9767832
#> Curl            0.6272325    0.4060797   0.2611775
#> fibreLength   304.5610227  196.5875983  25.5873638
#> fibreWidth     13.8917265   18.1127094   8.9148090
```
