import numpy as np
import matplotlib.pyplot as plt

def plot_style():

    # Set the style preferences:
    plt.style.use("seaborn-v0_8-darkgrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "legend.fontsize": 12,
        "figure.figsize": (8,6),
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "text.usetex": True
    })


def load_background_data():
    # Load the data from cosmology.txt
    filename = "./cosmology.txt"
    data = np.loadtxt(filename)

    return data


def load_mcmc_results():
    # Load the results from the MCMC fitting

    filename = "results_supernovafitting.txt"
    data = np.loadtxt(filename, skiprows=200)   # Skip the burnin of the first chains

    # Extract parameters
    chi2    = data[:,0]
    h       = data[:,1]
    Omega_M = data[:,2]
    Omega_K = data[:,3]

    return chi2, h, Omega_M, Omega_K
        
#===========================================#
# Now the functions for the different plots 
#===========================================#


def plot_luminosity_distance(filename):
    data        = np.loadtxt(filename, skiprows=1)     # Skip the header
    z           = data[:,0]                            # Redshift
    d_L         = data[:,1]                            # in Gyr
    errorbars   = data[:,2]                            # in Gyr

    cosmo_data = load_background_data()
    x           = cosmo_data[:,0]
    z_model     = np.exp(-x) - 1
    dL_model    = cosmo_data[:,11]   

    # We wish to plot d_L / z (remember the errorbars!):
    dL_over_z = d_L / z
    err_over_z = errorbars / z

    plt.figure()

    plt.errorbar(
        z,
        dL_over_z,
        yerr=err_over_z,
        fmt='o',
        markersize=4,
        capsize=3,
        label="Supernova data",
        color='blue',
    )

    plt.plot(
        z_model,
        dL_model/z_model,
        label="Fiducial model",
        color="red")

    plt.xlabel(r"$z$")
    plt.ylabel(r"$d_L(z)/z$ [Gyr]")
    plt.title(r"Luminosity Distance")
    plt.legend()

    plt.tight_layout()
    plt.show()




# Setting the style and plotting
if __name__ == "__main__":
    plot_style()
    plot_luminosity_distance("data/supernovadata.txt")

  








