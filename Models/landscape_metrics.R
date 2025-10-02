# Load required libraries
library(landscapemetrics)
library(terra)
library(sf)

# --- 1. File paths ---
# CHANGE THIS YEAR VALUE FOR DIFFERENT YEARS
year_to_process <- 2020
forest_path <- paste0("/Users/rexdevereux/Desktop/Python/HFLD/Data/FNF/", year_to_process, "_fnf.tif")
boundary_path <- "/Users/rexdevereux/Desktop/Python/HFLD/Data/Covar_boundary/covar_boundaries.gpkg"
output_dir <- "/Users/rexdevereux/Desktop/Python/HFLD/Data/FNF/results"
final_output_csv <- paste0("/Users/rexdevereux/Desktop/Python/HFLD/Data/FNF/", year_to_process, "_forest_metrics_with_metadata.csv")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory\n")
}

# --- 2. Load ADM1 boundary  ---
cat("Loading ADM1 boundary ...\n")
boundary <- st_read(boundary_path, quiet = TRUE)

# --- 3. CRS check and alignment ---
cat("Getting raster CRS info...\n")
raster_info <- rast(forest_path, lyrs = 1)
raster_crs <- crs(raster_info)
vector_crs <- st_crs(boundary)

cat("Checking CRS compatibility...\n")
if (!st_crs(boundary) == st_crs(raster_info)) {
  cat("CRS mismatch — reprojecting boundary to match raster...\n")
  boundary <- st_transform(boundary, crs(raster_info))
} else {
  cat("CRS matches. No reprojection needed.\n")
}

# --- 4. Add polygon index for merging ---
boundary$poly_id <- 1:nrow(boundary)

# --- 5. Process polygons and collect results ---
cat("Starting to process polygons one by one...\n")
start_time <- Sys.time()
print(start_time)

# List to store all results
all_results <- list()

# Process all polygons
start_poly <- 1
end_poly <- nrow(boundary)

for (i in start_poly:end_poly) {
  cat(sprintf("Processing polygon %d of %d\n", i, nrow(boundary)))
  
  # Extract single polygon
  poly <- boundary[i, ]
  
  # Get polygon extent
  poly_extent <- ext(st_bbox(poly))
  
  # Crop raster to the polygon extent
  forest_crop <- tryCatch({
    crop(rast(forest_path), poly_extent)
  }, error = function(e) {
    cat(sprintf("Error cropping raster for polygon %d: %s\n", i, e$message))
    all_results[[i]] <- NULL
    next
  })
  
  # Skip if the cropped raster is empty
  if (is.null(forest_crop) || terra::nlyr(forest_crop) == 0 || terra::ncell(forest_crop) == 0) {
    cat(sprintf("Polygon %d produced empty raster. Skipping.\n", i))
    all_results[[i]] <- NULL
    next
  }
  
  # Mask the cropped raster with the polygon
  forest_clip <- tryCatch({
    mask(forest_crop, vect(poly))
  }, error = function(e) {
    cat(sprintf("Error masking raster for polygon %d: %s\n", i, e$message))
    all_results[[i]] <- NULL
    next
  })
  
  # Skip if the clipped raster is empty
  if (is.null(forest_clip) || terra::nlyr(forest_clip) == 0 || terra::ncell(forest_clip) == 0) {
    cat(sprintf("Polygon %d produced empty masked raster. Skipping.\n", i))
    all_results[[i]] <- NULL
    next
  }
  
  # Check if forest class 1 exists in this polygon
  val_table <- terra::freq(forest_clip)
  has_forest <- any(val_table$value == 1)
  
  if (!has_forest) {
    cat(sprintf("Polygon %d has no forest (class 1). Adding zeros.\n", i))
    # Create empty result for this polygon
    empty_result <- data.frame(
      poly_id = poly$poly_id,
      layer = 1,
      level = "class",
      class = 1,
      id = NA,
      metric = c("ed", "pd", "area_mn"),
      value = 0
    )
    all_results[[i]] <- empty_result
  } else {
    # Calculate metrics
    metrics <- tryCatch({
      # Calculate all metrics
      all_metrics <- calculate_lsm(
        landscape = forest_clip,
        what = c("lsm_c_ed", "lsm_c_pd", "lsm_c_area_mn")
      )
      
      # Filter to keep only forest class (class 1)
      all_metrics[all_metrics$class == 1,]
    }, error = function(e) {
      cat(sprintf("Error calculating metrics for polygon %d: %s\n", i, e$message))
      return(NULL)
    })
    
    # Skip if metrics calculation failed
    if (is.null(metrics) || nrow(metrics) == 0) {
      cat(sprintf("Failed to calculate metrics for polygon %d.\n", i))
      all_results[[i]] <- NULL
      next
    }
    
    # Add polygon ID to results
    metrics$poly_id <- poly$poly_id
    
    # Convert area_mn from square meters to hectares for easier interpretation
    if ("area_mn" %in% metrics$metric) {
      area_indices <- which(metrics$metric == "area_mn")
      metrics$value[area_indices] <- metrics$value[area_indices] / 10000  # m² to ha
    }
    
    # Store results for this polygon
    all_results[[i]] <- metrics
  }
  
  # Clean up to free memory
  rm(forest_crop, forest_clip)
  gc()
}

end_time <- Sys.time()
cat("Processing completed in ", difftime(end_time, start_time, units = "mins"), " minutes\n")
cat("Processed all", end_poly - start_poly + 1, "polygons\n")

# --- 6. Combine all results ---
cat("Combining results...\n")

# Remove NULL elements and combine all results
valid_results <- all_results[!sapply(all_results, is.null)]
combined_results <- do.call(rbind, valid_results)

# Skip if no results were found
if (is.null(combined_results) || nrow(combined_results) == 0) {
  cat("No valid results found to combine.\n")
} else {
  # Merge with original metadata
  metrics_with_meta <- merge(
    x = combined_results,
    y = st_drop_geometry(boundary),
    by.x = "poly_id",
    by.y = "poly_id",
    all.x = TRUE
  )
  
  # Save final combined results
  write.csv(metrics_with_meta, final_output_csv, row.names = FALSE)
  cat(sprintf("Saved combined results to %s\n", final_output_csv))
  
  # Create a more readable summary with one row per polygon
  cat("Creating a summary table...\n")
  summary_file <- gsub("\\.csv$", "_summary.csv", final_output_csv)
  
  # Reshape data to wide format - using base R to avoid additional package dependencies
  # First, create a dataframe with one row per polygon
  unique_polys <- unique(metrics_with_meta$poly_id)
  summary_df <- data.frame(poly_id = unique_polys)
  
  # Add each metric as a column
  for (metric_name in unique(metrics_with_meta$metric)) {
    metric_data <- metrics_with_meta[metrics_with_meta$metric == metric_name, c("poly_id", "value")]
    colnames(metric_data)[2] <- metric_name
    summary_df <- merge(summary_df, metric_data, by="poly_id", all.x=TRUE)
  }
  
  # Remove duplicate rows
  summary_df <- unique(summary_df)
  
  # Add any additional columns from the boundary data
  if(ncol(boundary) > 1) {
    boundary_data <- st_drop_geometry(boundary)
    summary_df <- merge(
      x = summary_df,
      y = boundary_data,
      by = "poly_id",
      all.x = TRUE
    )
  }
  
  # Save summary
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat(sprintf("Saved summary to %s\n", summary_file))
}

cat("Script completed.\n")