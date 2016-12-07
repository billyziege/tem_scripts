import sys
import os
import argparse
import numpy as np
from scipy import constants
from phase_volume_6d import Phase6DVolume

parser = argparse.ArgumentParser(description='Write the mean of the phase volume file(s).')
parser.add_argument('phase_volume_files', nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('--header', dest="header", action='store_true', help='If flag is present, will print out the name of the columns.  Default is off.', default=False)
parser.add_argument('--filter_r_max', dest="r_max", type=float, help='Tells the script to filter the phase volume returning particles with r smaller than or eaual to r_max.', default=None)

args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

c = constants.physical_constants["speed of light in vacuum"][0] #m/sec by default
h = constants.physical_constants["Planck constant in eV s"][0]/1e6 #Planck's constant in MeV/s

phase_volume = Phase6DVolume()
for filepath in args.phase_volume_files:
  filename = os.path.basename(filepath)
  phase_volume.injectFile(filepath,mass=mass_of_electron,momentum_weight=args.number_of_electrons_per_macroparticle)

if args.r_max is not None:
  phase_volume = phase_volume.filterByRadius(max_r=args.r_max)

output = [len(phase_volume)*args.number_of_electrons_per_macroparticle]
header = ["N_e"]

cov_matrix = phase_volume.getCovarianceMatrix()
e_x = np.sqrt(cov_matrix.getSubDeterminant(["x","px"]))/(mass_of_electron)
output.append(e_x*1e6)
header.append("emittance_x")
e_y = np.sqrt(cov_matrix.getSubDeterminant(["y","py"]))/(mass_of_electron)
output.append(e_y*1e6)
header.append("emittance_y")
e_z = np.sqrt(cov_matrix.getSubDeterminant(["z","pz"]))/(mass_of_electron)
output.append(e_z*1e6)
header.append("emittance_z")

r_perp_sq = phase_volume.calcRperpSq()
output.append(np.sqrt(r_perp_sq))
header.append("sigma_r")

L_t = np.sqrt(r_perp_sq)*h*c/(mass_of_electron*e_x*2.0)
output.append(L_t)
header.append("L_t")

if args.header:
  print ",".join(header)
print ",".join([str(o) for o in output])
