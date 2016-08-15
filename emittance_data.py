import sys
import os
import argparse
from scipy import constants
from folder_convention import parse_folder_for_metadata
from cosy_output import CosyOutput
from simulation_output import SimulationOutput

parser = argparse.ArgumentParser(description='Write the three columns for the numbers data, N0, Ne, Fa')
parser.add_argument('list_of_directories', type=str, help='The path to a file containing all of directories in which the 2-Cosy/screen.txt are stored.', default=None)
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron
print "Ne,ex,ez,Fa";
with open(args.list_of_directories) as f:
  for line in f:
    line = line.rstrip()
    output =[]
    folder_metadata = parse_folder_for_metadata(line)
    cosy_output = CosyOutput()
    cosy_output.injectFile(os.path.join(line,"2-Cosy/screen.txt"))
    row_numbers = cosy_output.getRowNumbersWithKeyValuePair("time",120e-12,comparison = 'gt')
    if row_numbers == []:
      continue
    step = int(cosy_output.returnTableValue("step number",row_numbers[0]))
    number_of_macroparticles = int(cosy_output.returnTableValue("number macroparticles",row_numbers[0]))
    number_of_electrons = number_of_macroparticles*args.number_of_electrons_per_macroparticle
    files = []
    for filename in os.listdir(os.path.join(line,"2-Cosy")):
      if filename.startswith(str(step)+'-x'):
        files.append(os.path.join(line,"2-Cosy/"+filename))
    simulation_output = SimulationOutput()
    for path in files:
      simulation_output.injectFile(path)
    simulation_output.convertCoordinatesToBetterUnits(mass_of_macroparticle)
    simulation_output = simulation_output.boostCoordinates("z")
    ex = simulation_output.calcEmittance("x")/args.number_of_electrons_per_macroparticle
    ez = simulation_output.calcEmittance("z")/args.number_of_electrons_per_macroparticle
    if ex == 0 or ez == 0:
      continue
    output.append(str(number_of_electrons))
    output.append(str(ex))
    output.append(str(ez))
    output.append(folder_metadata["applied_field"])
    print ",".join(output)
    
