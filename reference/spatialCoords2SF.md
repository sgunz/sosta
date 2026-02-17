# Function to convert spatialCoords to an sf object

Function to convert spatialCoords to an sf object

## Usage

``` r
spatialCoords2SF(spe)
```

## Arguments

- spe:

  SpatialExperiment; a object of class `SpatialExperiment`

## Value

sf; Simple feature collection of geometry type POINT

## Examples

``` r
data(sostaSPE)
speSel <- sostaSPE[, sostaSPE[["imageName"]] == "image1"]
spatialCoords2SF(speSel)
#> Simple feature collection with 1845 features and 0 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 0.1061074 ymin: 0.1104626 xmax: 127.8393 ymax: 127.9601
#> CRS:           NA
#> First 10 features:
#>                     geometry
#> 1  POINT (113.9648 79.72171)
#> 2  POINT (18.82812 51.33354)
#> 3  POINT (119.7184 57.28533)
#> 4   POINT (38.5573 119.7797)
#> 5  POINT (7.772233 87.46652)
#> 6  POINT (18.21367 63.46682)
#> 7  POINT (74.94187 113.7856)
#> 8  POINT (82.93037 90.34805)
#> 9  POINT (40.93704 116.3168)
#> 10 POINT (39.38816 86.56549)
```
