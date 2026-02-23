import numpy as np
import matplotlib.pyplot as plt

def SupernovaFitting():

    # Load the data from the fitting

    filename = "results_supernovafitting.txt"
    data = np.loadtxt(filename, skiprows=200)   # Skip the burnin of the first chains

    # Extract parameters
    chi2 = data[0, :]
    h = data[1,:]
    Omega_M = data[2,:]
    Omega_K = data[3,:]

    #=====================================#
    # Now the plotting 
    #=====================================#

    # Set the style preferences:
    









