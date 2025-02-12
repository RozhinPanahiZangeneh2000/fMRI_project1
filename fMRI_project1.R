# Load required libraries
library(neurobase)  # For reading NIfTI fMRI files
library(tidyverse)  # For data manipulation
library(oro.nifti)  # Alternative for handling NIfTI
library(RNifti)     # Another efficient NIfTI handling package

# Define file paths (modify based on your actual directory)
fmri_files <- list.files(path = "path_to_fmri_data", pattern = "*.nii.gz", full.names = TRUE)
atlas_file <- "path_to_craddock_atlas.nii.gz"

# Load the Craddock atlas (assumes it's a 3D NIfTI file)
atlas <- readNIfTI(atlas_file, reorient = FALSE)

# Define Regions of Interest (ROIs)
roi_vPPC <- 52   # vPPC (ventral Posterior Parietal Cortex)
roi_V1 <- 80     # Primary Visual Cortex

# Initialize empty list to store extracted data
roi_data <- list()

# Loop through each participant's fMRI scan
for (i in 1:length(fmri_files)) {
  # Load fMRI data (4D NIfTI)
  fmri_data <- readNIfTI(fmri_files[i], reorient = FALSE)

  # Extract dimensions (should be 4D: x, y, z, time)
  dims <- dim(fmri_data)
  
  # Ensure the 4th dimension is considered (time)
  time_points <- dims[4]

  # Flatten the atlas into a vector (matching spatial coordinates)
  atlas_vec <- as.vector(atlas)
  
  # Flatten fMRI data into a matrix where rows = voxels, columns = time points
  fmri_mat <- matrix(fmri_data, nrow = prod(dims[1:3]), ncol = time_points)

  # Extract voxel indices for each ROI
  idx_vPPC <- which(atlas_vec == roi_vPPC)
  idx_V1 <- which(atlas_vec == roi_V1)

  # Extract mean intensity across time for each ROI
  vPPC_signal <- colMeans(fmri_mat[idx_vPPC, ], na.rm = TRUE)
  V1_signal <- colMeans(fmri_mat[idx_V1, ], na.rm = TRUE)

  # Store extracted data as a data frame
  roi_data[[i]] <- data.frame(
    Participant = paste0("Subject_", i),
    Time = 1:time_points,
    vPPC = vPPC_signal,
    V1 = V1_signal
  )
}

# Combine data from all participants
final_data <- bind_rows(roi_data)

# View first few rows
head(final_data)

# Save data to CSV for later analysis
write.csv(final_data, "fmri_roi_data.csv", row.names = FALSE)

# Plot comparison of average fMRI activity over time
ggplot(final_data, aes(x = Time)) +
  geom_line(aes(y = vPPC, color = "vPPC")) +
  geom_line(aes(y = V1, color = "V1")) +
  facet_wrap(~ Participant) +
  theme_minimal() +
  labs(title = "fMRI Signal in vPPC vs. Primary Visual Cortex",
       x = "Time (Volumes)", y = "Mean Signal Intensity") +
  scale_color_manual(values = c("vPPC" = "red", "V1" = "blue"))
