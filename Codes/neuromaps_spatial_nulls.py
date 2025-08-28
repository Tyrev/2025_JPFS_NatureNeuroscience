# -*- coding: utf-8 -*-
"""
Spin Test
=========================

Use spatial null models to test the correlation between two brain
annotations.
"""

# source ~/Documents/bin/activate

import nibabel as nib
import numpy as np
from neuromaps import images, nulls, stats, plotting, transforms

###############################################################################
# Load data

custom_map_path = "Use_3_TSPO_706_mirr_mRNA_ADNI_CorticalMasked.nii"
allen_map_path = "Use_interaction_estimate_zscore_CorticalMasked_Allen.nii"

custom_map = images.load_nifti(custom_map_path)
allen_map = images.load_nifti(allen_map_path)

print("Custom map:", custom_map)
print("Allen map:", allen_map)

###############################################################################
# Verify if resampling is necessary

custom_img = nib.load(custom_map_path)
allen_img = nib.load(allen_map_path)

print("Custom map shape:", custom_img.shape)
print("Allen map shape:", allen_img.shape)
print("Custom map affine:\n", custom_img.affine)
print("Allen map affine:\n", allen_img.affine)

# If the shapes and affine matrices match, proceed without resampling.
if custom_img.shape == allen_img.shape and (custom_img.affine == allen_img.affine).all():
    print("Both maps are already in the same space and resolution. No resampling needed.")
else:
    raise ValueError("Maps have different resolutions or orientations. Resampling may be required.")

###############################################################################
# Correlate the brain maps

corr = stats.compare_images(custom_map, allen_map)
print(f'Correlation: r = {corr:.02f}')

###############################################################################
# Convert volume to mni152_to_mni152 surface
mni152_custom_map = transforms.mni152_to_mni152(custom_map, target = '2mm')
mni152_allen_map = transforms.mni152_to_mni152(allen_map, target = '2mm')

# Check new shapes
print("Custom map (fsaverage) shape:", mni152_custom_map)
print("Allen map (fsaverage) shape:", mni152_allen_map)

###############################################################################
# Correlate the brain maps

corr = stats.compare_images(mni152_custom_map, mni152_allen_map)
print(f'Correlation: r = {corr:.02f}')

###############################################################################
# Extract data from images for statistical testing

custom_map_data = images.load_data(mni152_custom_map)
allen_map_data = images.load_data(mni152_allen_map)

print("Custom map data shape:", custom_map_data.shape)
print("Allen map data shape:", allen_map_data.shape)

###############################################################################
# Generate null models using spatial autocorrelation
# (Burt et al., 2020 method)

rotated = nulls.burt2020(custom_map_data, atlas='mni152', density = '2mm', n_perm=1000, seed=1234)

print("Generated null models shape:", rotated.shape)

###############################################################################
# Compute correlation and non-parametric p-value

corr, pval = stats.compare_images(custom_map_data, allen_map_data, nulls=rotated)
print(f'Correlation: r = {corr:.02f}, p = {pval:.04f}')


###############################################################################
# Visualization
import matplotlib.pyplot as plt
import seaborn as sns

# Compute correlations for null models
null_corrs = [stats.compare_images(rotated[..., i], allen_map_data) for i in range(rotated.shape[-1])]

# Plot the distribution of null correlations
sns.histplot(null_corrs, bins=30, kde=True, label="Null Distributions", color="blue")

# Plot the empirical correlation as a vertical line
plt.axvline(corr, color="red", linestyle="dashed", label=f"Empirical r = {corr:.02f}")

plt.xlabel("Correlation")
plt.ylabel("Density")
plt.legend()
plt.title("Empirical vs. Null Correlations")

# Save the figure
output_path = "null_vs_empirical_correlation.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")

# Close the figure to free memory
plt.close()

print(f"Figure saved as {output_path}")

###############################################################################
# Save table
import pandas as pd

# Create a DataFrame with null correlations
df_null = pd.DataFrame({"Source": "Null", "Correlation": null_corrs})

# Create a DataFrame for the empirical correlation
df_empirical = pd.DataFrame({"Source": ["Empirical"], "Correlation": [corr]})

# Combine both DataFrames
df = pd.concat([df_null, df_empirical], ignore_index=True)

# Save to CSV
output_csv = "correlation_data_long_format.csv"
df.to_csv(output_csv, index=False)

print(f"Table saved as {output_csv}")
