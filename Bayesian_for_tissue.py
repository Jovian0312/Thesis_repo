#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  8 23:15:29 2025

@author: u5501584
"""
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold
from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
from tqdm import tqdm

# Load data
df = pd.read_excel('/Users/u5501584/documents/Part_Results/Dangelo_embryo.xlsx')
U_data = ((df['Centroid Position'].values - df['Centroid Position'].values[0])*1e-6)  # Convert to micrometers
t_data = np.linspace(0, 59.7, len(U_data))
U_truth = ((df['Filtered'].values - df['Filtered'].values[0])*1e-6)

# Force parameters
F0 = 204e-12  # Force magnitude [N]
T = 10        # Total cycle period [s]
Ton = 5       # Force ON duration per cycle [s]

# Modified spring-dashpot model with numerical stability
def spring_dashpot_model(t, K, mu1, mu2, F0=450e-12, T=10, Ton=5):
    U = np.zeros_like(t)
    for i, ti in enumerate(t):
        t_cycle = ti % T
        exp_term1 = np.clip(np.exp(-K * t_cycle / mu1), 1e-10, 1e10)
        exp_term2 = np.clip(np.exp(-K * (t_cycle - Ton) / mu1), 1e-10, 1e10)
        exp_term3 = np.clip(np.exp(-K * Ton / mu1), 1e-10, 1e10)
        if t_cycle < Ton:
            U[i] = (F0 / mu2) * (t_cycle + (mu2 / K) * (1 - exp_term1))
        else:
            U[i] = (F0 / mu2) * (
                Ton - (mu2 / K) * exp_term2 * (1 - exp_term3)
            )
    return U

# Define the objective function
def objective(params, lambda_reg=1):
    K, mu1, mu2 = params['K'], params['mu1'], params['mu2']
    # Penalty for very small parameters to avoid numerical instability
    if K < 1e-6 or mu1 < 1e-6 or mu2 < 1e-6:
        return {'loss': 1e10 - 1e6 * (np.log(K + 1e-10) + np.log(mu1 + 1e-10) + np.log(mu2 + 1e-10)), 'status': STATUS_OK}
    try:
        U_pred = spring_dashpot_model(t_data, K, mu1, mu2, F0, T, Ton)
        mse = np.mean((U_pred - U_data) ** 2)
        # Quadratic regularization to penalize large parameters
        regularisation = lambda_reg * (K**2 + mu1**2 + mu2**2)
        # Logarithmic barrier for small parameters
        regularisation += lambda_reg * (-np.log(K) - np.log(mu1) - np.log(mu2))
        # Variance penalty to encourage reasonable model output
        variance_penalty = 1e-3 / (np.var(U_pred) + 1e-6)
        pbar.update(1)
        return {'loss': mse + regularisation + variance_penalty, 'status': STATUS_OK}
    except (OverflowError, ValueError):
        return {'loss': 1e10, 'status': STATUS_OK}

# Define the search space with uniform distribution
space = {
    'K': hp.uniform('K', 1e-4, 1e-2),    # µN/µm
    'mu1': hp.uniform('mu1', 1e-4, 1e-2), # µN·s/µm
    'mu2': hp.uniform('mu2', 1e-4, 1e-2)  # µN·s/µm
}

# Run Hyperopt optimization with progress bar
trials = Trials()
with tqdm(total=20000, desc="Optimizing") as pbar:
    best = fmin(
        fn=lambda params: objective(params, lambda_reg=1),
        space=space,
        algo=tpe.suggest,
        max_evals=20000,
        trials=trials
    )

# Extract best parameters
K_opt, mu1_opt, mu2_opt = best['K'], best['mu1'], best['mu2']
print("\nOptimized parameters (Hyperopt):")
print(f"K = {K_opt:.6f} µN/µm")
print(f"μ1 = {mu1_opt:.6f} µN·s/µm")
print(f"μ2 = {mu2_opt:.6f} µN·s/µm")

# Generate best-fit curve
t_fit = np.linspace(t_data.min(), t_data.max(), 500)
U_fit = spring_dashpot_model(t_fit, K_opt, mu1_opt, mu2_opt, F0, T, Ton)

# Calculate R-squared
residuals = U_data - spring_dashpot_model(t_data, K_opt, mu1_opt, mu2_opt, F0, T, Ton)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((U_data - np.mean(U_data))**2)
r_squared = 1 - (ss_res / ss_tot)
print(f"R-squared = {r_squared:.4f}")

# Cross-validation
kf = KFold(n_splits=5, shuffle=True, random_state=42)
mse_scores = []
for train_idx, test_idx in kf.split(t_data):
    t_train, t_test = t_data[train_idx], t_data[test_idx]
    U_train, U_test = U_data[train_idx], U_data[test_idx]
    U_pred = spring_dashpot_model(t_test, K_opt, mu1_opt, mu2_opt, F0, T, Ton)
    mse_scores.append(mean_squared_error(U_test, U_pred))
print(f"Cross-validated MSE: {np.mean(mse_scores):.6f} ± {np.std(mse_scores):.6f}")

# Plot results
plt.figure(figsize=(10, 6))
plt.scatter(t_data, U_data, label='Experimental Data', color='black', alpha=0.7)
plt.plot(t_fit, U_fit, label=f'Hyperopt Fit (R²={r_squared:.3f})', color='red', linewidth=2)
plt.xlabel('Time (s)', fontsize=15)
plt.ylabel('Displacement (µm)', fontsize=15)
plt.title('Spring-Dashpot Model Optimization with Hyperopt', fontsize=16)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot residuals
plt.figure(figsize=(10, 6))
plt.scatter(t_data, residuals, color='blue', alpha=0.7)
plt.axhline(0, color='black', linestyle='--')
plt.xlabel('Time (s)', fontsize=15)
plt.ylabel('Residuals (µm)', fontsize=15)
plt.title('Residuals of Optimized Fit', fontsize=16)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot force profile
force = np.array([F0 if (ti % T) < Ton else 0 for ti in t_fit])
plt.figure(figsize=(10, 6))
plt.plot(t_fit, force, label='Applied Force', color='green')
plt.xlabel('Time (s)', fontsize=15)
plt.ylabel('Force (N)', fontsize=15)
plt.title('Force Profile', fontsize=16)
plt.legend(fontsize=12)
plt.tick_params(axis='both', labelsize=12)
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot optimization loss trajectory
plt.figure(figsize=(10, 6))
losses = [trial['result']['loss'] for trial in trials.trials]
plt.plot(losses)
plt.xlabel('Iteration', fontsize=15)
plt.ylabel('Loss', fontsize=15)
plt.title('Optimization Loss Trajectory', fontsize=16)
plt.tick_params(axis='both', labelsize=12)
plt.grid(True)
plt.tight_layout()
plt.show()