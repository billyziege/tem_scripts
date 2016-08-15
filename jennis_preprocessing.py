import numpy as np
from scipy import constants
import argparse

def calcEmittance(mass, x, p):
    """ Function to calculate normalized emittance
    Args:
        mass (real) : mass of particles
        x (real array) : array of positions of particles
        p (real array) : array of momenta of particles
    Returns:
        emittance (real): value of emittance
    """
    mean_x = np.mean(x)
    mean_p = np.mean(p)
    mean_xp = np.mean(x * p)
    var_x = np.var(x, ddof = 1)
    var_p = np.var(p, ddof = 1)
    emittance = np.sqrt(var_x * var_p - (mean_xp - mean_x * mean_p)**2)/mass
    return emittance

def calcEnergySpread(mass, x, p):
    """ Function to calculate energy spread
    Args:
        mass (real) : mass of particles
        x (real array) : array of positions of particles
        p (real array) : array of momenta of particles
    Returns:
        gamma (real): value of Lorentz gamma
        E_spread (real): value of energy spread
    """
    mean_x = np.mean(x)
    mean_p = np.mean(p)
    mean_xp = np.mean(x * p)
    var_x = np.var(x, ddof = 1)
    var_p = np.var(p, ddof = 1)
    gamma = np.sqrt(1 + (mean_p/mass)**2)
    delta_p = np.sqrt(var_p - (mean_xp - mean_x * mean_p)**2/var_x)
    delta_gamma = np.sqrt(gamma**2 -1)/gamma * delta_p/mass
    E_spread = delta_gamma/gamma
    return gamma, E_spread

def calcBunchLength(mass, x, p):
    """ Function to calculate electron bunch length 
    Args:
        mass (real) : mass of particles
        x (real array) : array of positions of particles
        p (real array) : array of momenta of particles

    Returns:
        np.std(t_dist) (real): standard deviation of temporal distribution
    """
    speed_light = constants.physical_constants["speed of light in vacuum"][0]
    vel = (p/mass)/np.sqrt(1 + (p/mass)**2) * speed_light
    t_dist = x / vel
    return np.std(t_dist)

def analyze_cosy_run(filename, n_macro_part = 100):
    """ Function to analyze data from cosy simulation
    Args:
        filename (string): name of simulation file to open
        n_macro_part (int): how many electrons in one macro particle
    """
    me = constants.physical_constants["electron mass energy equivalent in MeV"][0]
    macro_mass = me * n_macro_part

    data = np.loadtxt(filename)
    emittance_T = calcEmittance(macro_mass, data[:,0], data[:,3])
    emittance_z = calcEmittance(macro_mass, data[:,2], data[:,5])
    gamma, energy_spread = calcEnergySpread(macro_mass, data[:,2], data[:,5])
    K = me * (gamma - 1)
    bunch_length = calcBunchLength(macro_mass, data[:,2], data[:,5])
    FWHM_bunch_length = bunch_length * 2 * np.sqrt(2.0 * np.log(2.0))

    print("Kinetic energy  = %.6f MeV" % K)
    print("Long emittance = %.3f nm " % (emittance_z * 1e9))
    print("Transv emittance = %.3f nm " % (emittance_T * 1e9))
    print("Energy spread \delta \gamma/ \gamma = %.3e " % energy_spread)
    print("Bunch length = %.3f ps" % (bunch_length * 1e12))
    print("FWHM Bunch length = %.3f ps" % (FWHM_bunch_length * 1e12))

    h = constants.physical_constants["Planck constant in eV s"][0]
    speed_light = constants.physical_constants["speed of light in vacuum"][0]
    constant = h*speed_light/(2.0*me*1e6)
    Ne = data.shape[0]*n_macro_part
    print("Ne %e" %Ne)
    coh_flux = constant**2 * Ne/emittance_T**2
    print("Coh. flux = %e " % coh_flux)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the analyze_cosy_run function Jenni wrote.')
    parser.add_argument('phase_volume_file', type=str, help='The path to the simultion output files containing the phase space.')
    parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)

    args = parser.parse_args()

    analyze_cosy_run(args.phase_volume_file,args.number_of_electrons_per_macroparticle)
