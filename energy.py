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

parser = argparse.ArgumentParser(description='Write the mean of the phase volume file(s).')
parser.add_argument('phase_volume_files', nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('-w','--warp', dest="warp", action="store_true", help='Tells the program that the output is from Warp.  The default is Cosy output.', default=False)
parser.add_argument('-e','--extraction_field', dest="extraction_field", type=float, help="The applied extraction field strength in MV/m.  Default is 1.", default=1)

args = parser.parse_args()

if args.warp:
  mass_of_electron = constants.physical_constants["electron mass"][0]
else:
  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron
e_charge = constants.physical_constants["elementary charge"][0]


phase_volume = Phase6DVolume(mass_of_macroparticle)
for filepath in args.phase_volume_files:
  filename = os.path.basename(filepath)
  phase_volume.injectFile(filepath,mass=mass_of_macroparticle)
ke = phase_volume.getKineticEnergy()
fpe = phase_volume.getExtractionFieldPotentialEnergy(e_charge,args.extraction_field*1000000)
spe = phase_volume.getSelfPotentialEnergy(e_charge)
output = [ke,fpe,spe,ke+fpe+spe]
print ",".join([str(o) for o in output])
