# Create a Point Pattern on a Simulated Tissue Image

This function generates a spatial point pattern with different types of
points (`A`, `B`, `C`) distributed over the simulated tissue structure.

## Usage

``` r
createPointPatternTissue(tissueImage, intA, intB, intCInA, intCInB)
```

## Arguments

- tissueImage:

  Matrix; A binary matrix representing the simulated tissue.

- intA:

  Numeric; Intensity of type "A" points (points per unit area) on tissue
  regions.

- intB:

  Numeric; Intensity of type "B" points (points per unit area) on
  non-tissue regions.

- intCInA:

  Numeric; Intensity of type "C" points placed in extended regions
  around tissue.

- intCInB:

  Numeric; Intensity of type "C" points placed within tissue.

## Value

A `ppp` object representing the spatial point pattern.

## Examples

``` r
tissueImage <- simulateTissueBlobs(128, 100, 7)
createPointPatternTissue(tissueImage, 0.1, 0.1, 0.005, 0.005)
#> Marked planar point pattern: 1822 points
#> marks are of storage type  ‘character’
#> window: rectangle = [0, 128] x [0, 128] units
```
