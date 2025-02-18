# Simulate 3 images
set.seed(123)
tissue_image1 <- simulateTissueBlobs(128, 100, 7)
point_pattern1 <- createPointPatternTissue(tissue_image1, 0.1, 0.1, 0.005, 0.005)

tissue_image2 <- simulateTissueBlobs(128, 90, 9)
point_pattern2 <- createPointPatternTissue(tissue_image2, 0.12, 0.08, 0.003, 0.002)

tissue_image3 <- simulateTissueBlobs(128, 110, 6)
point_pattern3 <- createPointPatternTissue(tissue_image3, 0.095, 0.12, 0.006, 0.006)

# Create data frames
df1 <- as.data.frame(point_pattern1)
df1$image_name <- "image1"

df2 <- as.data.frame(point_pattern2)
df2$image_name <- "image2"

df3 <- as.data.frame(point_pattern3)
df3$image_name <- "image3"

# Combine
df_all <- rbind(df1,df2,df3)
colnames(df_all)[3] <- "cell_type"


sostaSPE <- SpatialExperiment::SpatialExperiment(
    colData = df_all,
    spatialCoordsNames = c("x", "y"))


usethis::use_data(sostaSPE, overwrite = TRUE)



