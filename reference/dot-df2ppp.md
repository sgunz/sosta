# Function to convert `data.frame` to `ppp` object

Assumes that the `data.frame` is the output of
[`.SPE2df()`](https://sgunz.github.io/sosta/reference/dot-SPE2df.md).
Column order is important!

## Usage

``` r
.df2ppp(df, xName, yName, marks = NULL)
```

## Arguments

- df:

  data.frame; with x, y coordinates, image, and categorical mark
  information.

- xName:

  character; column name of x coordinate

- yName:

  character; column name of y coordinate

- marks:

  character; column name of the mark variable

## Value

ppp; object of type `ppp`

## See also

[`.SPE2df`](https://sgunz.github.io/sosta/reference/dot-SPE2df.md),
[`as.ppp`](https://rdrr.io/pkg/spatstat.geom/man/as.ppp.html)

## Examples

``` r
data(sostaSPE)
df <- .SPE2df(sostaSPE, marks = "cellType", imageCol = "imageName")
ppp <- .df2ppp(df, xName = "x", yName = "y", marks = "cellType")
```
