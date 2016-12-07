import sys
import os
import argparse
import numpy as np
from scipy import constants
from cosy_output import CosyOutput
from coordinate_vector_3d import Cartesian3DVector
from phase_volume_6d import Phase6DVolume
from my_covariance_matrix import MyEditableCovarianceMatrix
from binned_phase_volume import BinnedPhase6DVolume, CylindricalBinnedPhase6DVolume, SlicedPhase6DVolume

parser = argparse.ArgumentParser(description='Write the mean of the phase volume file(s).')
parser.add_argument('phase_volume_files', nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('-b','--number_bins', dest="number_bins", type=int, help='Number of bins per dimension. Default is none.', default=None)
parser.add_argument('-w','--warp', dest="warp", action="store_true", help='Tells the program that the output is from Warp.  The default is Cosy output.', default=False)
parser.add_argument('--number_bins_x', dest="number_bins_x", type=int, help='Number of bins in x-direction. Default is none.', default=None)
parser.add_argument('--number_bins_y', dest="number_bins_y", type=int, help='Number of bins in y-direction. Default is none.', default=None)
parser.add_argument('--number_bins_z', dest="number_bins_z", type=int, help='Number of bins in z-direction. Default is none.', default=None)
parser.add_argument('--number_bins_rho', dest="number_bins_rho", type=int, help='Number of bins in rho-direction. Default is none.', default=None)
parser.add_argument('--cylindrical', dest="cylindrical", action="store_true", help='Calculate and report the temperature in cylindrical coordinates.', default=False)
parser.add_argument('--header', dest="header", action="store_true", help='Prints out the names of the columns', default=False)
parser.add_argument('--sliced', dest="sliced", action="store_true", help='Calculate and report the temperature in sliced coordinates.', default=False)


args = parser.parse_args()

if args.warp:
  mass_of_electron = constants.physical_constants["electron mass"][0]
  #kb = 1.38064852e-23#in J/K
  kb = constants.physical_constants["elementary charge"][0]#in J/eV
else:
  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  #kb = 8.6173324e-11#in MeV/K
  kb=1e-6#in MeV/eV
mass_in_kg = constants.physical_constants["electron mass"][0]
epsilon = 8.854187817e-12#F/m
charge_of_electron = constants.physical_constants["elementary charge"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

phase_volume = Phase6DVolume(mass_of_electron)
for filepath in args.phase_volume_files:
  filename = os.path.basename(filepath)
  if args.warp:
    phase_volume.injectFile(filepath,mass=mass_of_electron,fieldnames=["x","y","z","px","py","pz","vx","vy","vz"])
  else:
    phase_volume.injectFile(filepath,mass=mass_of_electron,momentum_weight=args.number_of_electrons_per_macroparticle)
if args.cylindrical:
  binned_phase_volume = CylindricalBinnedPhase6DVolume(phase_volume,number_bins_1d=args.number_bins,number_bins_rho=args.number_bins_rho,number_bins_z=args.number_bins_z)
elif args.sliced:
  if args.number_bins_z is None:
    args.number_bins_z = args.number_bins
  binned_phase_volume = SlicedPhase6DVolume(phase_volume,number_bins_z=args.number_bins_z)
else:
  binned_phase_volume = BinnedPhase6DVolume(phase_volume,number_bins_1d=args.number_bins,number_bins_x=args.number_bins_x,number_bins_y=args.number_bins_y,number_bins_z=args.number_bins_z,x_mean=0,y_mean=0)
binned_phase_volume.k_b = kb
binned_phase_volume.electrons_per_macroparticle = args.number_of_electrons_per_macroparticle
binned_phase_volume.epsilon = epsilon
binned_phase_volume.charge = charge_of_electron
binned_phase_volume.mass_in_kg = mass_in_kg
binned_phase_volume.printAverageStatistics(args.header)
