import sys
import os
import argparse
from folder_convention import parse_folder_for_metadata
from cosy_output import CosyOutput
from simulation_output import SimulationOutput

parser = argparse.ArgumentParser(description='Write the three columns for the numbers data, N0, Ne, Fa')
parser.add_argument('list_of_directories', type=str, help='The path to a file containing all of directories in which the 2-Cosy/screen.txt are stored.', default=None)
args = parser.parse_args()

print "Ne,ex,ez,Fa";
with open(args.list_of_directories) as f:
  for line in f:
    line = line.rstrip()
    output =[]
    folder_metadata = parse_folder_for_metadata(line)
    cosy_output = CosyOutput()
    cosy_output.injectScreenFile(os.path.join(line,"2-Cosy/screen.txt"))
    row_numbers = cosy_output.returnRowNumbersWithKeyValuePair("time",120e-12,comparison = 'gt')
    if row_numbers == []:
      continue
    step = int(cosy_output.returnTableValue("step number",row_numbers[0]))
    files = []
    for filename in os.listdir(os.path.join(line,"2-Cosy")):
      if filename.startswith(str(step)+'-x'):
        files.append(os.path.join(line,"2-Cosy/"+filename))
    simulation_output = SimulationOutput()
    for path in files:
      simulation_output.injectCoordinatesFile(path)
    ex = simulation_output.returnEmmitance("x")
    ez = simulation_output.returnEmmitance("z")
    if ex == 0 or ez == 0:
      continue
    output.append(str(100 * int(cosy_output.returnTableValue("number macroparticles",row_numbers[0]))))
    output.append(str(ex))
    output.append(str(ez))
    output.append(folder_metadata["applied_field"])
    print ",".join(output)
    
