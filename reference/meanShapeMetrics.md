# Calculate mean shape metrics of a set of polygons

Calculate mean shape metrics of a set of polygons

## Usage

``` r
meanShapeMetrics(totalShapeMetricMatrix)
```

## Arguments

- totalShapeMetricMatrix:

  matrix of shape metrics

## Value

matrix; matrix of mean shape metrics

## Examples

``` r
data(sostaSPE)
struct <- reconstructShapeDensityImage(sostaSPE,
    marks = "cellType", imageCol = "imageName",
    imageId = "image1", markSelect = "A", dim = 500
)
shapeMetrics <- totalShapeMetrics(struct)
meanShapeMetrics(shapeMetrics)
#>                          [,1]
#> Area             2673.2396360
#> Compactness         0.3252460
#> Eccentricity        0.6916357
#> Circularity         0.6500011
#> Solidity            0.7906035
#> Curl                0.4314966
#> fibreLength       175.5786616
#> fibreWidth         13.6397483
#> numberStructures    3.0000000
#> totalArea        8019.7189079
```
