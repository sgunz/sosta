# Simulate 3 images
set.seed(123)
tissueImage1 <- simulateTissueBlobs(128, 100, 7)
pointPattern1 <- createPointPatternTissue(tissueImage1, 0.1, 0.1, 0.005, 0.005)

tissueImage2 <- simulateTissueBlobs(128, 90, 9)
pointPattern2 <- createPointPatternTissue(tissueImage2, 0.12, 0.08, 0.003, 0.002)

tissueImage3 <- simulateTissueBlobs(128, 110, 6)
pointPattern3 <- createPointPatternTissue(tissueImage3, 0.095, 0.12, 0.006, 0.006)

# Create data frames
df1 <- as.data.frame(pointPattern1)
df1$imageName <- "image1"

df2 <- as.data.frame(pointPattern2)
df2$imageName <- "image2"

df3 <- as.data.frame(pointPattern3)
df3$imageName <- "image3"

# Combine
dfAll <- rbind(df1, df2, df3)
colnames(dfAll)[3] <- "cellType"


sostaSPE <- SpatialExperiment::SpatialExperiment(
    colData = dfAll,
    spatialCoordsNames = c("x", "y")
)


usethis::use_data(sostaSPE, overwrite = TRUE)
