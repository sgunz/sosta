# Example SpatialExperiment Object with Simulated Tissue Images and Point Patterns

This dataset contains a simulated `SpatialExperiment` object
(`sostaSPE`) representing three tissue images, each with a corresponding
spatial point pattern. The point patterns contain different cell types
(`A`, `B`, and `C`), distributed according to simulated tissue
structures.

## Usage

``` r
sostaSPE
```

## Format

A `SpatialExperiment` object with the following structure:

- x:

  Numeric; x-coordinate of each point (cell location).

- y:

  Numeric; y-coordinate of each point (cell location).

- cell_type:

  Factor; Cell type assigned to each point (A, B, or C).

- image_name:

  Factor; Identifier for the tissue image (`image1`, `image2`, or
  `image3`).

## Details

The dataset was generated as follows:

- Three tissue images were simulated using
  [`simulateTissueBlobs()`](https://sgunz.github.io/sosta/reference/simulateTissueBlobs.md).

- Spatial point patterns were created for each tissue using
  [`createPointPatternTissue()`](https://sgunz.github.io/sosta/reference/createPointPatternTissue.md).

- The point pattern data was converted into a `SpatialExperiment` object
  with spatial coordinates.
