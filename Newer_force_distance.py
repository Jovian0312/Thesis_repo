#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 19:25:23 2025

@author: u5501584
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

def mixed_decay(x, a, b, c, d):
    """Combined exponential and linear decay function"""
    return a * np.exp(-b * x) + c * x + d

def exp_decay(x, a, b, c):
    """Exponential decay function"""
    return a * np.exp(-b * x) + c

def linear_decay(x, a, b):
    """Linear decay function"""
    return a * x + b

folder_path = '/Users/u5501584/documents/results/'
files = sorted([f for f in os.listdir(folder_path) if f.endswith('.xlsx')])
for file in files:
    file_path = os.path.join(folder_path, file)
    try:
        df = pd.read_excel(file_path)
        
        # Check if required columns exist
        required_columns = ['Force', 'distance']
        if not all(col in df.columns for col in required_columns):
            missing = [col for col in required_columns if col not in df.columns]
            print(f"File {file} is missing required columns: {missing}")
            continue
            
        # Data preprocessing
        force_mean = df['Force'].replace(0, np.nan).mean()
        df['Force'] = df['Force'].replace(0, force_mean)
        
        # Create a unique identifier if Track_ID doesn't exist
        if 'Track_ID' not in df.columns:
            df['Track_ID'] = 1  # Assign all points to the same track if no ID exists
        
        df['Filtered_Force'] = df['Force'].rolling(window=5, min_periods=1).mean()

        unique_ids = df['Track_ID'].unique()
        dfs = {id_: df[df['Track_ID'] == id_] for id_ in unique_ids} 
        
        plt.figure(figsize=(8,6))
        #colours = ['#ff5a5e', '#c31e23', '#0d7d87','#348EBB']
        colours = ['#00b0be','#8fd7d7','#c99b38','#348EBB']
        
        # Collect all data points for fitting
        all_x = []
        all_y = []
        
        for idx, (bead_id, group) in enumerate(dfs.items()):
            colour = colours[idx % len(colours)] 
            x_data = group['distance'].values
            y_data = group['Filtered_Force'].values
            
            # Plot scatter points
            plt.scatter(x_data, y_data, label=f'Track {bead_id}', s=20, alpha=0.7, c=colour)
            
            # Collect all data points
            all_x.extend(x_data)
            all_y.extend(y_data)
        
        # Convert to numpy arrays and sort for fitting
        all_x = np.array(all_x)
        all_y = np.array(all_y)
        sort_idx = np.argsort(all_x)
        all_x = all_x[sort_idx]
        all_y = all_y[sort_idx]
        
        # Try different fitting functions and select the best one
        best_fit = None
        best_r2 = -np.inf
        fit_info = ""
        
        # Try mixed decay fit
        try:
            popt, pcov = curve_fit(mixed_decay, all_x, all_y, 
                                  p0=[max(all_y)-min(all_y), 1, -1, min(all_y)], 
                                  maxfev=5000)
            y_pred = mixed_decay(all_x, *popt)
            ss_res = np.sum((all_y - y_pred)**2)
            ss_tot = np.sum((all_y - np.mean(all_y))**2)
            r_squared = 1 - (ss_res / ss_tot)
            
            if r_squared > best_r2:
                best_r2 = r_squared
                best_fit = (mixed_decay, popt)
                fit_info = (f'Mixed Fit:\n'
                           f'y = {popt[0]:.2f}e^(-{popt[1]:.3f}x) + {popt[2]:.3f}x + {popt[3]:.2f}\n'
                           f'R² = {r_squared:.3f}')
        except Exception as e:
            print(f"Mixed decay fit failed for {file}: {str(e)}")
        
        # Try exponential decay fit
        try:
            popt, pcov = curve_fit(exp_decay, all_x, all_y, 
                                  p0=[max(all_y)-min(all_y), 0.01, min(all_y)], 
                                  maxfev=5000)
            y_pred = exp_decay(all_x, *popt)
            ss_res = np.sum((all_y - y_pred)**2)
            ss_tot = np.sum((all_y - np.mean(all_y))**2)
            r_squared = 1 - (ss_res / ss_tot)
            
            if r_squared > best_r2:
                best_r2 = r_squared
                best_fit = (exp_decay, popt)
                fit_info = (f'Exp Fit:\n'
                           f'y = {popt[0]:.2f}e^(-{popt[1]:.3f}x) + {popt[2]:.2f}\n'
                           f'R² = {r_squared:.3f}')
        except Exception as e:
            print(f"Exponential decay fit failed for {file}: {str(e)}")
        
        # Try linear fit
        try:
            popt, pcov = curve_fit(linear_decay, all_x, all_y, 
                                  p0=[-1, max(all_y)], 
                                  maxfev=5000)
            y_pred = linear_decay(all_x, *popt)
            ss_res = np.sum((all_y - y_pred)**2)
            ss_tot = np.sum((all_y - np.mean(all_y))**2)
            r_squared = 1 - (ss_res / ss_tot)
            
            if r_squared > best_r2:
                best_r2 = r_squared
                best_fit = (linear_decay, popt)
                fit_info = (f'Linear Fit:\n'
                           f'y = {popt[0]:.3f}x + {popt[1]:.2f}\n'
                           f'R² = {r_squared:.3f}')
        except Exception as e:
            print(f"Linear fit failed for {file}: {str(e)}")
        
        # Plot the best fit if found
        if best_fit is not None:
            fit_func, popt = best_fit
            x_fit = np.linspace(min(all_x), max(all_x), 100)
            y_fit = fit_func(x_fit, *popt)
            
            plt.plot(x_fit, y_fit, 'k--', linewidth=2, label='Best Fit')
            
            # Position the equation in the upper middle with a white background
            plt.text(0.5, 0.95, fit_info, transform=plt.gca().transAxes, 
                    fontsize=15, color='black', verticalalignment='top',
                    horizontalalignment='center',
                    bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray', boxstyle='round'))
        else:
            plt.text(0.5, 0.95, "Could not fit any curve to data", 
                    transform=plt.gca().transAxes, fontsize=15, color='red',
                    horizontalalignment='center')
        
        plt.xlabel('Distance (μm)', fontsize=20) 
        plt.ylabel('Force (pN)', fontsize=20)
        plt.legend(title='Key: Track IDs', loc='upper right', fontsize='small', title_fontsize='medium')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"Error processing file {file}: {str(e)}")