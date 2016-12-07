import sys
import os
import argparse
from scipy import constants
import numpy
from folder_convention import parse_folder_for_metadata
from cosy_output import CosyOutput

parser = argparse.ArgumentParser(description='Write the two columns for the velocity data: z,v/c')
parser.add_argument('screen_file', type=str, help='The path screen.txt are stored.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
args = parser.parse_args()

print "z,std_z";
mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron
z_conv = 10**3;
std_z_conv = 10**6;
cosy_output = CosyOutput()
cosy_output.injectFile(args.screen_file)
#cosy_output.addGammaField(mass_of_macroparticle)
for row in cosy_output.rows:
  output = []
  lorentz_gamma = numpy.sqrt(1 + (row.getCellWithFieldname("p_z").getValue()/mass_of_macroparticle)**2)
  output.append(str(row.getCellWithFieldname("z").getValue()*z_conv))
  output.append(str(row.getCellWithFieldname("std_z").getValue()*std_z_conv/lorentz_gamma))
  print ",".join(output)
