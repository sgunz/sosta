# Spatial Omic Structure Analysis

![](inst/extdata/sosta.png)

`sosta` (spatial omics structure analysis) is a framework to
reconstruct, characterize and compare spatial structures from spatial
omics data.

`sosta` builds on existing R packages for spatial analysis such as
`spatstat`, `sf` and the [Bioconductor](https://bioconductor.org/)
environment by integrating with packages `SpatialExperiment` and
`SpatialFeatureExperiment`, which facilitates the implementation of
additional custom metrics that fit the needs of users.

## Installation

Installation of the latest [Bioconductor
version](https://bioconductor.org/packages/sosta) can be done as
follows:

``` r

if (!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
}
BiocManager::install("sosta")
```

You can install the development version of `sosta` from GitHub with:

``` r

# install.packages("devtools")
devtools::install_github("sgunz/sosta")
```

## Feedback

We are happy to get your feedback. Please send it via email to [Samuel
Gunz](https://www.mls.uzh.ch/en/research/robinson/groupmembers/samuel-gunz.html)
or open a issue on [GitHub](https://github.com/sgunz/sosta/issues).

## Citation

If you use `sosta` please consider citing:

> Gunz, Samuel, Helena Lucia Crowell, and Mark D. Robinson. “Analysis of
> Anatomical Multi-Cellular Structures from Spatial Omics Data Using
> Sosta”, p. 2025.10.13.682065. bioRxiv, 14 Oct. 2025,
> <https://doi.org/10.1101/2025.10.13.682065>.
