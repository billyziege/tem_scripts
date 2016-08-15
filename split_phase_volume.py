import sys
import os
import argparse
import numpy as np
from scipy import constants
from cosy_output import CosyOutput
from coordinate_vector_3d import Cartesian3DVector
from phase_volume_6d import Phase6DVolume
from my_covariance_matrix import MyEditableCovarianceMatrix
from binned_phase_volume import BinnedPhase6DVolume, CylindricalBinnedPhase6DVolume

parser = argparse.ArgumentParser(description='Write files for the binned phase volumes aggregated over z.')
parser.add_argument('phase_volume_files', nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space.')
parser.add_argument('-o', '--output',dest="output",type=str,help='The beginning of the path to the output file to which _iz.txt will be added.  Defaut is output.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('-b','--number_bins', dest="number_bins", type=int, help='Number of bins per dimension. Default is none.', default=None)
#parser.add_argument('--number_bins_x', dest="number_bins_x", type=int, help='Number of bins in x-direction. Default is none.', default=None)
#parser.add_argument('--number_bins_y', dest="number_bins_y", type=int, help='Number of bins in y-direction. Default is none.', default=None)
parser.add_argument('--number_bins_z', dest="number_bins_z", type=int, help='Number of bins in z-direction. Default is none.', default=None)
parser.add_argument('--number_bins_rho', dest="number_bins_rho", type=int, help='Number of bins in rho-direction. Default is none.', default=None)
parser.add_argument('--cylindrical', dest="cylindrical", action="store_false", help='Use cylindrical coordinates.  Default is true.', default=True)

args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

phase_volume = Phase6DVolume(mass_of_macroparticle)
for filepath in args.phase_volume_files:
  filename = os.path.basename(filepath)
  phase_volume.injectFile(filepath,mass=mass_of_macroparticle)
if args.cylindrical:
  binned_phase_volume = CylindricalBinnedPhase6DVolume(phase_volume,number_bins_1d=args.number_bins,number_bins_rho=args.number_bins_rho,number_bins_z=args.number_bins_z)
for ir in range(0,len(binned_phase_volume.phase_volumes)):
  output_file = args.output+"_"+str(ir)+".txt"
  with open(output_file,'w') as f:
    for iz in range(0,len(binned_phase_volume.phase_volumes[ir])):
      if len(binned_phase_volume.phase_volumes[ir][iz]) > 0: 
        f.write(str(binned_phase_volume.phase_volumes[ir][iz])+ "\n")
