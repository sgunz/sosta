# Simulate Tissue Blobs

This function generates a simulated tissue-like structure using a
Gaussian blur technique.

## Usage

``` r
simulateTissueBlobs(size, seedNumber, clumpSize)
```

## Arguments

- size:

  Integer; The size (width and height) of the simulated tissue image.

- seedNumber:

  Integer; The number of random seed points used to generate tissue
  blobs.

- clumpSize:

  Numeric; The standard deviation (sigma) of the Gaussian blur applied
  to generate tissue clumps.

## Value

A binary matrix representing the simulated tissue structure.

## Examples

``` r
tissueImage <- simulateTissueBlobs(128, 100, 7)
image(tissueImage)

```
