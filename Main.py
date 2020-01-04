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

### Reading the data in.
# The experiment consists of four configurations, each conducted with trials in triplicate.
# HF = high flow, LF = low flow, C = with cooling jacket, N = no cooling jacket, 1/2/3 = trial number
# For each trial, there is a .csv (2-D) that contains
configs = ['HFC', 'HFN', 'LFC', 'LFN']
trials = ['1', '2', '3']

# Reading these in as ndarrays and storing each as an object element in a list. This method is better suited than
# pandas MultiIndex, xarray, or 3D numpy array because the matrices have varying dimensions.
# Rows: configurations (x4)      Cols: trials (x3)
S_matrices = [[np.genfromtxt('Resources\\' + configs[k] + trials[i] + '.csv', delimiter=",") for i in range(0, 3)] for k
              in range(0, 4)]

### Setting up variables for physical quantities.
porosity = 0.74

# Bed Height and Volume
# Follows the same order as the trials array below.
# All heights, diameters, and volumes are on a meter basis.
bed_height = np.array([[5.30, 5.30, 5.70], [5.90, 5.50, 5.00], [4.40, 5.30, 4.80], [5.80, 5.70, 5.10]]) * 1E-2  # meters
bed_volume = np.array([[0.000001896593738, 0.000001896593738, 0.000002039732888],
                       [0.000002000201035, 0.000001864594185, 0.000001695085623],
                       [0.000001500770909, 0.000001807746776, 0.000001637204628],
                       [0.000001848407641, 0.000001816538544, 0.00000162532396]])  # m^3
bed_diameter = np.array([[6.75, 6.75, 6.75], [6.75, 6.57, 6.57], [6.59, 6.59, 6.59],
                         [6.37, 6.37, 6.37]]) / 1000  # meters
flowrate = np.array([[6.2, 6.6, 7.2], [6.0, 6.3, 6.7], [2.0, 2.0, 2.0], [2.0, 2.0, 2.0]]) / 6E4  # m^3 / s
R = 8.206E-5  # m3*atm*K&^-1?mol^-1
T = 273.15 + 24  # Kelvin, assumed constant
# tau_factor is used to convert time to the dimensionless tau
tau_factor = np.divide(flowrate, np.multiply(porosity, bed_volume))

# Getting tau ranges for each trial in case we want them again for plotting.
tau_range = []
mean_Z = []  # using capital Z to denote dimensionless scaled values

# Extracting Psi values.
psi = np.full((4, 3), -1)
i = 1
for n in range(4):
    mean_Z.append([])
    tau_range.append([])
    for k in range(3):
        mean_Z[n].append(find_mean_Z(S_matrices[n][k], 0.5, 1))
        tau_range[n].append(np.arange(1, len(S_matrices[n][k][:, 0]) * 16 * tau_factor[n][k], tau_factor[n][k] * 15))
        # TODO clean the mean_Z and corresponding elements of tau_range, then use polyfit with deg = 1 for psi (slope)
        # populate the corresponding loc in meanZ
        # with the Z array
        # the Z array is a 1D array (axis=time) where each element contains the average Z-value along the silica column
        # that corresponds to 50% saturation at that time step
        # here the mean Z-value is determined by taking the mean value of the Z-values observed where
        # S = 50%, within a small threshold

#         meanZ{i} = MOD_extract_z(sCellMatrix{k,n},tauFactor(i),50,1);
#         tauRange(i) = {1:15*tauFactor(i):(length(sCellMatrix{k,n}(:,1)).*15.*tauFactor(i))};
#         % Clean out the NaN values from the meanZ array and augment the
#         % tauRange to reflect this, store in new arrays.
#         cleanMeanZMachine = [];
#         cleanMeanTauRange = [];
#         for x = 1:length(meanZ{i})
#             if ~isnan(meanZ{i}(x))
#                 cleanMeanZMachine = [cleanMeanZMachine, meanZ{i}(x)];
#                 cleanMeanTauRange = [cleanMeanTauRange, tauRange{i}(x)];
#             end
#         end


#         fit = polyfit(cleanMeanTauRange, cleanMeanZMachine, 1);
#         psi(i) = fit(1);
#         i = i+1;
#     end
# end
# i = 1;
