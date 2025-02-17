# Load required libraries
library(neurobase)  # For reading NIfTI fMRI files
library(tidyverse)  # For data manipulation and visualization
library(oro.nifti)  # Alternative for handling NIfTI
library(RNifti)     # Another efficient NIfTI handling package
library(psych)      # For correlation analysis

# Define file paths (modify based on actual directory)
fmri_files <- list.files(path = "path_to_fmri_data", pattern = "*.nii.gz", full.names = TRUE)
atlas_file <- "path_to_craddock_atlas.nii.gz"

# Load the Craddock atlas (assumes it's a 3D NIfTI file)
atlas <- readNIfTI(atlas_file, reorient = FALSE)

# Define Regions of Interest (ROIs)
roi_vPCC <- 200  # vPCC (ventral Posterior Cingulate Cortex)
roi_V1 <- 80     # Primary Visual Cortex

# Initialize empty list to store extracted data
roi_data <- list()

# Loop through each participant's fMRI scan
for (i in 1:length(fmri_files)) {
  # Load fMRI data (4D NIfTI)
  fmri_data <- readNIfTI(fmri_files[i], reorient = FALSE)

  # Extract dimensions (should be 4D: x, y, z, time)
  dims <- dim(fmri_data)
  time_points <- dims[4]

  # Flatten atlas and fMRI data
  atlas_vec <- as.vector(atlas)
  fmri_mat <- matrix(fmri_data, nrow = prod(dims[1:3]), ncol = time_points)

  # Extract voxel indices for each ROI
  idx_vPCC <- which(atlas_vec == roi_vPCC)
  idx_V1 <- which(atlas_vec == roi_V1)

  # Extract mean intensity across time for each ROI
  vPCC_signal <- colMeans(fmri_mat[idx_vPCC, ], na.rm = TRUE)
  V1_signal <- colMeans(fmri_mat[idx_V1, ], na.rm = TRUE)

  # Compute functional connectivity (Pearson correlation)
  connectivity <- cor(V1_signal, vPCC_signal, use = "complete.obs")
  
  # Store extracted data as a data frame
  roi_data[[i]] <- data.frame(
    Participant = paste0("Subject_", i),
    Time = 1:time_points,
    vPCC = vPCC_signal,
    V1 = V1_signal,
    Connectivity = connectivity
  )
}

# Combine data from all participants
final_data <- bind_rows(roi_data)

# Save data to CSV for later analysis
write.csv(final_data, "fmri_lsd_connectivity.csv", row.names = FALSE)

# Plot V1 and vPCC activity over time
ggplot(final_data, aes(x = Time)) +
  geom_line(aes(y = vPCC, color = "vPCC")) +
  geom_line(aes(y = V1, color = "V1")) +
  facet_wrap(~ Participant) +
  theme_minimal() +
  labs(title = "fMRI Signal in vPCC vs. Primary Visual Cortex under LSD",
       x = "Time (Volumes)", y = "Mean Signal Intensity") +
  scale_color_manual(values = c("vPCC" = "purple", "V1" = "blue"))

# Display correlation results
print("Functional Connectivity between V1 and vPCC (Pearson correlation):")
print(final_data %>% select(Participant, Connectivity) %>% distinct())

