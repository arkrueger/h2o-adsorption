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
This project uses a finite difference method to solve two linked partial differential equations describing the gas 
phase (diffusion) and solid phase (adsorption) behavior.

Equation 1: dG/dtau = -dG/dZ - kappa*(G - S)
Equation 2: dS/dtau = psi*kappa*(G - S)

Equations 1 and 2 are in a dimensionless form, scaled such that tau and Z take on values between 0 and 1. This is 
critical to the stability of the finite difference. 

The dimensionless term tau is defined as:
tau === t*N_a/(porosity*C*V)
Where t = time, N_a = gas molar flowrate, porosity = silica porosity, C = molar gas concentration from gas law (PV=nRT)
    V = silica bed volume


'''
import numpy as np
import pandas as pd
import numerical_tools_adsorption as nta


### Setting up variables for physical quantities.
# Bed Height and Volume
# Follows this order: rows are configuration, columns are trials.
# All heights, diameters, and volumes are on a meter basis.
bed_height = np.array([[5.30, 5.30, 5.70], [5.90, 5.50, 5.00], [4.40, 5.30, 4.80],
                       [5.80, 5.70, 5.10]]).flatten() * 1E-2  # meters
bed_volume = np.array([[0.000001896593738, 0.000001896593738, 0.000002039732888],
                       [0.000002000201035, 0.000001864594185, 0.000001695085623],
                       [0.000001500770909, 0.000001807746776, 0.000001637204628],
                       [0.000001848407641, 0.000001816538544, 0.00000162532396]]).flatten()  # m^3
bed_diameter = np.array([[6.75, 6.75, 6.75], [6.75, 6.57, 6.57], [6.59, 6.59, 6.59],
                         [6.37, 6.37, 6.37]]).flatten() / 1000  # meters
flowrate = np.array([[6.2, 6.6, 7.2], [6.0, 6.3, 6.7], [2.0, 2.0, 2.0], [2.0, 2.0, 2.0]]).flatten() / 6E4  # m^3 / s
R = 8.206E-5  # m3*atm*K&^-1?mol^-1
temp = np.array([297, 299, 299, 307, 308, 308, 293, 292, 292, 302, 305, 306])  # Kelvin
porosity = 0.74
delta_t = 15  # seconds (every 15 seconds, an image was taken of the silica column)
# tau_factor is used to convert time to the dimensionless tau.
tau_factor = np.divide(flowrate, np.multiply(porosity, bed_volume))


### Reading the data in.
# The experiment consists of four configurations, each conducted with trials in triplicate.
# HF = high flow, LF = low flow, C = with cooling jacket, N = no cooling jacket, 1/2/3 = trial number
# For each trial, there is a .csv (2-D) that contains
configs = ['HFC', 'HFN', 'LFC', 'LFN']
trials = ['1', '2', '3']
trial_name = [i + k for i in configs for k in trials]
trial_num = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

# Dict of dicts to organize data. Outer dict maps the experimental trials, inner dict contains the associated data.
attributes = dict.fromkeys(['S_matrix', 'trial_num', 'trial_name', 'bed_height', 'bed_volume', 'bed_diameter',
                            'flowrate', 'tau_factor', 'tau_range', 'psi', 'Z_50', 'Z_lower', 'Z_upper', 'kappa'])
# If the variable or object already exists, its name is placed here to set up the data dict by comprehension notation.
attribute_keys = ['trial_num', 'trial_name', 'bed_height', 'bed_volume', 'bed_diameter', 'flowrate', 'tau_factor']
# Here we iterate through the trial_name array (containing the outer dict keys) and the arrays containing the physical
# quantities at the same rate and order. The k index is just to access the arrays those quantities are stored in. This
# ensures that the right quantities are associated with the right trial.
# That's confusing and there's probably a better way to map it. Change the following line at your own peril.
data = {trial_name[i]: {k: globals()[k][i] for k in attribute_keys} for i in range(0, 12)}

# Populate the dict with the S matrices from the resource files.
for trial in data:
    data[trial]['S_matrix'] = np.genfromtxt('Resources\\' + trial + '.csv', delimiter=",")
    data[trial]['porosity'] = porosity

### Analysis
# Extracting psi values and collecting Z_35, Z_50, and Z_65.
for i in data:
    # Temporary reference to the current trial to make code more readable.
    trial = data[i]
    # Initialize the time axis (tau) - used for graphing and curve fitting.
    num_of_frames = len(trial['S_matrix'][:, 0])
    trial['tau_range'] = np.linspace(0.0, num_of_frames * delta_t * trial['tau_factor'], num_of_frames)

    # find_mean_Z to get the average Z position at each time step where saturation was 50% (call this Z_50).
    # (Capital Z denotes the dimensionless scaled value, as opposed to lowercase z, which is pixels.)
    trial['Z_50'] = nta.find_mean_Z(trial, 0.5, show_graphs=1)

    # psi is the slope of Z_50 vs. tau, so find a linear fit using only the finite values (there are typically NaNs).
    idx = np.isfinite(trial['Z_50'])
    trial['psi'] = np.polyfit(trial['tau_range'][idx], trial['Z_50'][idx], 1)[0]
    # While we're at it, let's collect the Z_35 and Z_65.
    trial['Z_35'] = nta.find_mean_Z(trial, 0.35, show_graphs=0)
    trial['Z_65'] = nta.find_mean_Z(trial, 0.65, show_graphs=0)

# This could all be lumped into a single loop. Here I'm breaking it up into two loops in the interest of displaying
# similar graphs together.
# Fit kappa.
for i in data:
    # Temporary reference to the current trial to make code more readable.
    trial = data[i]
    trial['kappa'] = nta.fit_kappa(trial, 4, show_graphs=1)

# TODO Plot Petrovic-Thodos correlation.


