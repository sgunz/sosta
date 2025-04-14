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
