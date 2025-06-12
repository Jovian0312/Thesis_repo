#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 23:32:18 2025

@author: u5501584
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Single folder path
folder_path = '/Users/u5501584/documents/results/PEG_avg_3.2'

# Colors and labels for each dataset
colors = ['#2E8B57', '#FF4500']  # Green for 5μm, OrangeRed for 10μm
labels = ['5μm Mean Force', '10μm Mean Force']
sd_labels = ['5μm ±1 SD', '10μm ±1 SD']

# Initialize plot
plt.figure(figsize=(8, 6))

# Process data for 5μm and 10μm separately
for idx, size in enumerate(['5um', '10um']):
    # Filter files containing the size indicator
    files = sorted([f for f in os.listdir(folder_path) if f.endswith('.xlsx') and size in f.lower()])

    # Process each file and aggregate data
    all_data = []
    for i in range(len(files)):
        file_path = os.path.join(folder_path, files[i])
        
        # Load and clean data
        df = pd.read_excel(file_path)
        force_mean = df['Force'].replace(0, np.nan).mean()
        df['Force'] = df['Force'].replace(0, force_mean)
        df.to_excel(file_path, index=False)  # Save cleaned version

        # Reload and apply rolling average
        df = pd.read_excel(file_path)
        df['Filtered_Force'] = df['Force'].rolling(window=2).mean()
        all_data.append(df)

    # Combine all data for the size
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Bin distances
    bin_size = 15  # μm
    combined_df['binned_distance'] = (combined_df['distance'] // bin_size) * bin_size

    # Group by binned distance and calculate mean and std
    stats_df = combined_df.groupby('binned_distance')['Filtered_Force'].agg(['mean', 'std']).reset_index()

    # Plot mean with shaded standard deviation
    plt.plot(stats_df['binned_distance'], stats_df['mean'],
             color=colors[idx], linewidth=2, marker='o', markersize=5, label=labels[idx])
    
    plt.fill_between(stats_df['binned_distance'],
                     stats_df['mean'] - stats_df['std'],
                     stats_df['mean'] + stats_df['std'],
                     color=colors[idx], alpha=0.3, label=sd_labels[idx])

# Customize plot
plt.xlabel('Distance (μm)', fontsize=18)
plt.ylabel('Filtered Force (pN)', fontsize=18)
plt.tick_params(axis='both', labelsize=14)  # Increase tick label size
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()