import sys
import os
import argparse
import numpy as np
from scipy import constants
from cosy_output import CosyOutput
from coordinate_vector_3d import Cartesian3DVector
from phase_volume_6d import Phase6DVolume
from my_covariance_matrix import MyEditableCovarianceMatrix

parser = argparse.ArgumentParser(description='Write the mean of the phase volume file(s).')
parser.add_argument('phase_volume_files', nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('-o','--output_name', dest="output_name", type=str, help='Specify the output file name.  Do not include directory paths in this option.  It will be ignored.  See --output_dir to rectify this.  This defaults to the string before "-" of the first phase volume file plus "_mean.txt".', default=None)
parser.add_argument('-d','--output_dir', dest="output_dir", type=str, help='Specify the output directory.  This defaults to the dirname of the first phase volume file.', default=None)



args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

output_dir = args.output_dir
if output_dir is None:
  output_dir = os.path.dirname(args.phase_volume_files[0])
if args.output_name is None:
  pieces = os.path.basename(args.phase_volume_files[0]).split('-')
  output_name = pieces[0]+"_mean.txt"
else:
  output_name = os.path.basename(args.output_name)
output_path = os.path.join(output_dir,output_name)

phase_volume = Phase6DVolume()
for filepath in args.phase_volume_files:
  filename = os.path.basename(filepath)
  phase_volume.injectFile(filepath,mass=mass_of_macroparticle)
x_mean = Cartesian3DVector(x=phase_volume.getMean(["x"]),y=phase_volume.getMean(["y"]),z=phase_volume.getMean(["z"]))
p_mean = Cartesian3DVector(x=phase_volume.getMean(["px"]),y=phase_volume.getMean(["py"]),z=phase_volume.getMean(["pz"]))
with open(output_path,'w') as f:
  output = []
  output.append(x_mean.x)
  output.append(x_mean.y)
  output.append(x_mean.z)
  output.append(p_mean.x)
  output.append(p_mean.y)
  output.append(p_mean.z)
  f.write(" ".join([str(o) for o in output])+"\n")
