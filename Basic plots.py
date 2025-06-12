#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 01:14:45 2025

@author: u5501584
"""
import pandas as pd
import matplotlib.pyplot as plt

# Read the Excel file
df = pd.read_excel('~/Documents/Results/Vout_current.xlsx')

# Create the scatter plot without a legend label
plt.scatter(df['v_out'], df['current'], color='#33777c', marker='o', s=10)



# Add labels and title
plt.xlabel('Current (Amps)', fontsize=15)
plt.ylabel('Temperature(Degrees celcius)', fontsize=15)
plt.legend()

# Show the plot
plt.show()