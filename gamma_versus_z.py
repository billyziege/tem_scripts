import sys
import os
import argparse
from scipy import constants
from folder_convention import parse_folder_for_metadata
from cosy_output import CosyOutput
from simulation_output import SimulationOutput

parser = argparse.ArgumentParser(description='Write the two columns for the velocity data: z,gamma_z')
parser.add_argument('directory', type=str, help='The path where the step files and screen.txt are stored.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
args = parser.parse_args()

print "z,gamma_z,gamma_z_divided_by_std_z";
mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron
cosy_output = CosyOutput()
cosy_output.injectFile(os.path.join(args.directory,"screen.txt"))
cosy_output.addVelocityField(mass_of_macroparticle)
for row in cosy_output.rows:
  step_number = row.getCellWithFieldname("step number").getValue()
  if step_number % 10 == 0:
    output = []
    output.append(str(row.getCellWithFieldname("z").getValue()))
    simulation_output = SimulationOutput()
    for filename in os.listdir(args.directory):
      if filename.startswith(str(step_number)+"-x-"):
        pathname = os.path.join(args.directory,filename)
        if os.path.isfile(pathname):
          simulation_output.injectFile(pathname)
    gamma_z = simulation_output.calcGamma("z")
    output.append(str(gamma_z))
    std_z = row.getCellWithFieldname("std_z").getValue()*10**-6
    output.append(str(gamma_z/std_z))
    print ",".join(output)
