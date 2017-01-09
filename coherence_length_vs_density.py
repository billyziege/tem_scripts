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
parser.add_argument('--gamma', dest="gamma", type=float, help='The lorentz factor for the distribution.  Defaults to 1.2.', default=1.2)
parser.add_argument('-w','--warp', dest="warp", action="store_true", help='Tells the program that the output is from Warp.  The default is Cosy output.', default=False)
parser.add_argument('--input_format', dest="input_format", type=str, 
                    help='Describes the format of the input file.  Right now, only "xv" is supported.' +
                     '"xp" is supported by default as is "xpv".',default=None)

args = parser.parse_args()

if args.warp:
  mass_of_electron = constants.physical_constants["electron mass"][0]
else:
  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

c = constants.physical_constants["speed of light in vacuum"][0] #m/sec by default
h = constants.physical_constants["Planck constant in eV s"][0]/1e6 #Planck's constant in MeV/s

phase_volume = Phase6DVolume()
for filepath in args.phase_volume_files:
  filename = os.path.basename(filepath)
  if args.input_format == "xv":
    phase_volume.injectFile(filepath,mass=mass_of_electron,fieldnames=["x","y","z","vx","vy","vz"],velocity_weight=1./c)
  elif args.warp:
    phase_volume.injectFile(filepath,mass=mass_of_electron,momentum_weight=args.number_of_electrons_per_macroparticle,fieldnames=["x","y","z","px","py","pz","vx","vy","vz"])
  else:
    phase_volume.injectFile(filepath,mass=mass_of_electron,momentum_weight=args.number_of_electrons_per_macroparticle)

if args.r_max is not None:
  phase_volume = phase_volume.filterByRadius(max_r=args.r_max)

output = [len(phase_volume)*args.number_of_electrons_per_macroparticle]
header = ["N_e"]

cov_matrix = phase_volume.getCovarianceMatrix()

sd_x =  np.sqrt(cov_matrix.getCovarianceElement("x","x"))
output.append(sd_x*1e3)
header.append("sd_x")
sd_vx =  np.sqrt(cov_matrix.getCovarianceElement("px","px"))/(args.gamma*mass_of_electron)
output.append(sd_vx)
header.append("sd_vx")
e_x = np.sqrt(cov_matrix.getSubDeterminant(["x","px"]))/(mass_of_electron*c)
output.append(e_x*1e6)
header.append("ex")

sd_y =  np.sqrt(cov_matrix.getCovarianceElement("y","y"))
output.append(sd_y*1e3)
header.append("sd_y")
sd_vy =  np.sqrt(cov_matrix.getCovarianceElement("py","py"))/(args.gamma*mass_of_electron)
output.append(sd_vy)
header.append("sd_vy")
e_y = np.sqrt(cov_matrix.getSubDeterminant(["y","py"]))/(mass_of_electron*c)
output.append(e_y*1e6)
header.append("ey")

sd_z =  np.sqrt(cov_matrix.getCovarianceElement("z","z"))
output.append(sd_z*1e3)
header.append("sd_z")
sd_vz =  np.sqrt(cov_matrix.getCovarianceElement("pz","pz"))/(args.gamma*mass_of_electron)
output.append(sd_vz)
header.append("sd_vz")
e_z = np.sqrt(cov_matrix.getSubDeterminant(["z","pz"]))/(mass_of_electron*c)
output.append(e_z*1e6)
header.append("ez")

output.append(e_x*e_y*e_z*1e18)
header.append("exeyez")

r_perp_sq = phase_volume.calcRperpSq()
output.append(np.sqrt(r_perp_sq))
header.append("sigma_r")

L_t = np.sqrt(r_perp_sq)*h*c/(mass_of_electron*e_x*2.0)
output.append(L_t)
header.append("L_t")

B4D = len(phase_volume)*args.number_of_electrons_per_macroparticle/(e_x*e_y)
output.append(B4D)
header.append("B4D")

B6D = B4D/e_z
output.append(B6D)
header.append("B6D")

if args.header:
  print ", ".join(header)
print ",".join([str(o) for o in output])
