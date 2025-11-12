#!/usr/bin/env python3
"""
Create a multi-library vis.heart.h5ad from the single library version
by duplicating it with different library IDs
"""

import scanpy as sc
import numpy as np
import pandas as pd
import anndata

# Load the single library
print("Loading single library...")
adata_single = sc.read_h5ad("inst/testdata/vis.heart.single.h5ad")
print(f"Single library: {adata_single.n_obs} spots, {adata_single.n_vars} genes")

# Create two copies with different library IDs
lib1_name = "HCAHeartST11702008"
lib2_name = "HCAHeartST11702009"

print(f"\nCreating multi-library version with 2 libraries...")

# Copy 1
adata1 = adata_single.copy()
adata1.obs['library_id'] = lib1_name
adata1.obs_names = [f"{lib1_name}_{i}" for i in range(adata1.n_obs)]

print(f"Library 1 ({lib1_name}): {adata1.n_obs} spots")

# Copy 2 - with slight variation
adata2 = adata_single.copy()
adata2.obs['library_id'] = lib2_name
adata2.obs_names = [f"{lib2_name}_{i}" for i in range(adata2.n_obs)]

# Add some noise to make it slightly different (optional, for realism)
if adata2.X is not None:
    # Add small random variation (5% noise)
    np.random.seed(43)
    noise = np.random.randn(*adata2.X.shape) * 0.05
    adata2.X = adata2.X + noise

print(f"Library 2 ({lib2_name}): {adata2.n_obs} spots")

# Concatenate
print("\nConcatenating libraries...")
adata_multi = anndata.concat([adata1, adata2], join='inner', label='library_id_orig')

# Make sure library_id is set correctly
if 'library_id' not in adata_multi.obs.columns:
    # Use the concatenated labels
    if 'library_id_orig' in adata_multi.obs.columns:
        adata_multi.obs['library_id'] = adata_multi.obs['library_id_orig']
    else:
        # Manually assign based on obs_names
        adata_multi.obs['library_id'] = [name.split('_')[0] for name in adata_multi.obs_names]

print(f"Multi-library: {adata_multi.n_obs} spots, {adata_multi.n_vars} genes")
print(f"Libraries: {adata_multi.obs['library_id'].unique()}")

# Handle spatial data if it exists
if 'spatial' in adata_multi.obsm:
    print("Spatial coordinates found in obsm")

if 'spatial' in adata_multi.uns:
    print("Spatial metadata found in uns")
    print(f"  Keys: {list(adata_multi.uns['spatial'].keys())}")

    # Update spatial metadata for both libraries if needed
    # Make sure both library names are represented
    if lib1_name not in adata_multi.uns['spatial']:
        # Copy the first available spatial metadata
        first_key = list(adata_multi.uns['spatial'].keys())[0]
        adata_multi.uns['spatial'][lib1_name] = adata_multi.uns['spatial'][first_key].copy()

    if lib2_name not in adata_multi.uns['spatial']:
        first_key = list(adata_multi.uns['spatial'].keys())[0]
        adata_multi.uns['spatial'][lib2_name] = adata_multi.uns['spatial'][first_key].copy()

# Save
print("\nSaving multi-library version...")
adata_multi.write_h5ad("inst/testdata/vis.heart.h5ad", compression='gzip')

# Check file size
import os
size_mb = os.path.getsize("inst/testdata/vis.heart.h5ad") / (1024**2)
print(f"Saved: vis.heart.h5ad ({size_mb:.1f} MB)")

print("\n=== Summary ===")
print(f"Single library: inst/testdata/vis.heart.single.h5ad")
print(f"  - {adata_single.n_obs} spots, {adata_single.n_vars} genes")
print(f"Multi library: inst/testdata/vis.heart.h5ad")
print(f"  - {adata_multi.n_obs} spots, {adata_multi.n_vars} genes")
print(f"  - {len(adata_multi.obs['library_id'].unique())} libraries")

print("\nDone!")
