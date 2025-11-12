#!/usr/bin/env python3
"""
Downsample vis.heart.h5ad test data
Goal: Reduce spots to half and genes to 10k for faster vignette building
"""

import scanpy as sc
import numpy as np
import pandas as pd
import os

# Set random seed for reproducibility
np.random.seed(42)

# Path to original data
orig_path = "inst/testdata/vis.heart.h5ad"

if not os.path.exists(orig_path):
    print(f"Error: {orig_path} not found")
    exit(1)

print("Loading original h5ad file...")
adata = sc.read_h5ad(orig_path)

print(f"Original data: {adata.n_obs} spots, {adata.n_vars} genes")
print(f"Layers: {list(adata.layers.keys())}")
print(f"Obsm keys: {list(adata.obsm.keys())}")

# Check if library_id exists
if 'library_id' in adata.obs.columns:
    print(f"\nFound library_id column")
    libraries = adata.obs['library_id'].unique()
    print(f"Libraries: {libraries}")
    n_libs = len(libraries)
else:
    print("\nNo library_id found, treating as single library")
    libraries = None
    n_libs = 1

# --- Create single library version ---
print("\n=== Creating single library version ===")

if libraries is not None and n_libs > 0:
    # Take first library
    lib_name = libraries[0]
    print(f"Selecting library: {lib_name}")
    adata_single = adata[adata.obs['library_id'] == lib_name].copy()
else:
    adata_single = adata.copy()

print(f"Single library: {adata_single.n_obs} spots, {adata_single.n_vars} genes")

# Downsample spots (half)
n_spots = adata_single.n_obs
n_spots_keep = int(np.ceil(n_spots / 2))
print(f"Downsampling spots: {n_spots} -> {n_spots_keep}")

spots_idx = np.random.choice(adata_single.n_obs, size=n_spots_keep, replace=False)
adata_single = adata_single[spots_idx, :].copy()

# Downsample genes (top 10k by variance)
n_genes = adata_single.n_vars
if n_genes > 10000:
    print(f"Selecting top 10k variable genes from {n_genes}...")

    # Calculate gene variance from X (use the data matrix)
    if adata_single.X is not None:
        gene_vars = np.array(adata_single.X.var(axis=0)).flatten()
    else:
        print("Warning: X is None, trying layers")
        if 'data' in adata_single.layers:
            gene_vars = np.array(adata_single.layers['data'].var(axis=0)).flatten()
        else:
            # Use first available layer
            first_layer = list(adata_single.layers.keys())[0]
            print(f"Using layer: {first_layer}")
            gene_vars = np.array(adata_single.layers[first_layer].var(axis=0)).flatten()

    # Get top 10k genes by variance
    top_genes_idx = np.argsort(-gene_vars)[:10000]
    adata_single = adata_single[:, top_genes_idx].copy()
else:
    print(f"Keeping all {n_genes} genes")

print(f"Single library downsampled: {adata_single.n_obs} spots, {adata_single.n_vars} genes")

# Save single library
print("Saving single library...")
adata_single.write_h5ad("inst/testdata/vis.heart.single.h5ad", compression='gzip')
print(f"Saved: vis.heart.single.h5ad")

# --- Create multi-library version ---
print("\n=== Creating multi-library version ===")

if libraries is not None and n_libs > 1:
    # Keep only first 2 libraries
    n_libs_keep = min(2, n_libs)
    libs_to_keep = libraries[:n_libs_keep]
    print(f"Keeping {n_libs_keep} libraries: {libs_to_keep}")

    adata_multi = adata[adata.obs['library_id'].isin(libs_to_keep)].copy()
else:
    # Use the single library version
    adata_multi = adata.copy()

print(f"Multi-library: {adata_multi.n_obs} spots, {adata_multi.n_vars} genes")

# Downsample spots per library (half)
if libraries is not None and n_libs > 1:
    # Sample half from each library
    selected_idx = []
    for lib in libs_to_keep:
        lib_mask = adata_multi.obs['library_id'] == lib
        lib_idx = np.where(lib_mask)[0]
        n_lib_spots = len(lib_idx)
        n_keep = int(np.ceil(n_lib_spots / 2))

        print(f"  Library {lib}: {n_lib_spots} -> {n_keep} spots")

        selected = np.random.choice(lib_idx, size=n_keep, replace=False)
        selected_idx.extend(selected)

    adata_multi = adata_multi[selected_idx, :].copy()
else:
    # Single library or no library info - sample half overall
    n_spots = adata_multi.n_obs
    n_spots_keep = int(np.ceil(n_spots / 2))
    print(f"Downsampling spots: {n_spots} -> {n_spots_keep}")

    spots_idx = np.random.choice(adata_multi.n_obs, size=n_spots_keep, replace=False)
    adata_multi = adata_multi[spots_idx, :].copy()

# Downsample genes (top 10k by variance)
n_genes = adata_multi.n_vars
if n_genes > 10000:
    print(f"Selecting top 10k variable genes from {n_genes}...")

    # Calculate gene variance
    if adata_multi.X is not None:
        gene_vars = np.array(adata_multi.X.var(axis=0)).flatten()
    else:
        if 'data' in adata_multi.layers:
            gene_vars = np.array(adata_multi.layers['data'].var(axis=0)).flatten()
        else:
            first_layer = list(adata_multi.layers.keys())[0]
            gene_vars = np.array(adata_multi.layers[first_layer].var(axis=0)).flatten()

    top_genes_idx = np.argsort(-gene_vars)[:10000]
    adata_multi = adata_multi[:, top_genes_idx].copy()

print(f"Multi-library downsampled: {adata_multi.n_obs} spots, {adata_multi.n_vars} genes")

# Save multi-library
print("Saving multi-library...")
adata_multi.write_h5ad("inst/testdata/vis.heart.h5ad", compression='gzip')
print(f"Saved: vis.heart.h5ad")

# Check final file sizes
print("\n=== Final file sizes ===")
files_to_check = [
    "inst/testdata/vis.heart.h5ad",
    "inst/testdata/vis.heart.single.h5ad"
]

for f in files_to_check:
    if os.path.exists(f):
        size_mb = os.path.getsize(f) / (1024**2)
        print(f"  {os.path.basename(f)}: {size_mb:.1f} MB")

print("\nAll done!")
