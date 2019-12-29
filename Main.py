'''
Alex Krueger
2019
Numerical Solutions for Modeling H2O Adsorption on a Silica Column

Key Concepts:
    Numerical Solutions
    Finite Difference Methods
    Dependent Partial Differential Equations
'''

### Background:
'''
This program uses a finite difference method to solve two linked partial differential equations describing the gas 
phase (diffusion) and solid phase (adsorption) behavior.

Equation 1: dG/dtau = -dG/dZ - kappa*(G - S)
Equation 2: dS/dtau = psi*kappa*(G - S)

'''
import numpy as np
import pandas as pd

### Reading the data in.
# The experiment consists of four configurations, each conducted with trials in triplicate.
# HF = high flow, LF = low flow, C = with cooling jacket, N = no cooling jacket, 1/2/3 = trial number
configs = ['HFC', 'HFN', 'LFC', 'LFN']
trials = ['1', '2', '3']

# Reading these in as pandas DataFrames and storing each as an object element in a list.
# This method is better suited than MultiIndex, xarray, or 3D numpy array because the matrices have varying dimensions.
s_matrices = [[pd.read_csv('Resources\\'+configs[k]+trials[i]+'.csv') for i in range(0,3)] for k in range(0,4)]







