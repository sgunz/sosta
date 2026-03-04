# sosta 1.3.1

* added argument to simulate noise in `createPointPatternTissue`
* improved intensity value threshold estimation for reconstruction in `.intensityThreshold`

# sosta 1.2.0

* Bioc 3.22 release version

# sosta 1.1.4

* the argument `imageCol` is now not longer needed in reconstruction functions
* `assingCellsToStructures` for large datasets and 1 sample

# sosta 1.1.3

* added sticker and preprint reference

# sosta 1.1.2

* switch from `GenomeInfoDb` to new `Seqinfo` package 

# sosta 1.1.1

* added new function `minCellTypeStructDist`
* added possibility to reconstruct everything around an object (`complement = TRUE`)
* relaxed assertion of having at least two cells of each cell type in reconstruction to having at least one of one cell type

# sosta 1.0.1

* fix problems in when assigning cell to structures

# sosta 0.99.9

* Clarification of the transformation used in the diabetic islets vignette.

# sosta 0.99.8

* Small changes in the vignettes and description.

# sosta 0.99.7

* Small updates in vignettes.
* Updated news file for 0.99.6

# sosta 0.99.6

* Improved parallel computation by reducing memory overhead by using data frames within `mclapply` instead of `SpatialExperiment` objects.
* Introduced new functions `.SPE2df` and `.SPE2ppp` to convert to different objects for reducing memory overhead. 

# sosta 0.99.5

* New Analysis Functions: Added `assignCellsToStructures`, `calculateBorderMetrics`, `calculateCellTypeFractions`, `calculateDistanceToBorder`, `calculateDistanceToTissueRegion`, and `calculateTissueRegionFractions` to improve spatial and structural analysis.  
* Simulated data set: Introduced a new data set for testing and demonstrating package capabilities.  
* Vignette update: Expanded documentation to include multi-sample and condition-based analysis.  
* Bug fixes & optimisations: Addressed various issues and improved performance.  
* Code cleanup: Enhanced readability, maintainability, and structure. Consistent use of camelCase for function arguments.

# sosta 0.99.4

* Set number of cores to 1 to pass build on Windows.

# sosta 0.99.3

Second round of [review](https://github.com/Bioconductor/Contributions/issues/3584#issuecomment-2568132189)

* Changed required R version to >= 4.4.0
* Fixed typos
* Updated news file to match versions

# sosta 0.99.2

Addressed first [review](https://github.com/Bioconductor/Contributions/issues/3584#issuecomment-2404878963)

* Added R version and selective imports
* Functional programming to avoid code repetition
* Added unit tests and input checks
* Updated the vignette with installation info

# sosta 0.99.1

Initial Bioconductor submission
