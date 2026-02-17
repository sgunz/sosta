# Function to convert `SpatialExperiment` object to a data frame

Function to convert `SpatialExperiment` object to a data frame

## Usage

``` r
.SPE2df(spe, imageCol = NULL, marks = NULL, colNames = FALSE)
```

## Arguments

- spe:

  SpatialExperiment; a object of class `SpatialExperiment`

- imageCol:

  character; name of a column in `colData` that corresponds to the image

- marks:

  character; name of column in `colData` with categorical marks

- colNames:

  logical; extract `colnames` from `SpatialExperiment`

## Value

data.frame with x, y coordinates, image, and categorical mark
information

## Examples

``` r
data(sostaSPE)
.SPE2df(sostaSPE, marks = "cellType", imageCol = "imageName") |> head()
#>            x         y cellType imageName
#> 1 113.964828  79.72171        A    image1
#> 2  18.828120  51.33354        A    image1
#> 3 119.718375  57.28533        A    image1
#> 4  38.557299 119.77972        A    image1
#> 5   7.772233  87.46652        A    image1
#> 6  18.213670  63.46682        A    image1
```
