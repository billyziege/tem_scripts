import sys
import os
import argparse
from folder_convention import parse_folder_for_metadata
from cosy_output import CosyOutput

parser = argparse.ArgumentParser(description='Write the three columns for the numbers data, N0, Ne, Fa')
parser.add_argument('list_of_directories', type=str, help='The path to a file containing all of directories in which the 2-Cosy/screen.txt are stored.', default=None)
args = parser.parse_args()

print "N0,Ne,Fa";
with open(args.list_of_directories) as f:
  for line in f:
    line = line.rstrip()
    output =[]
    folder_metadata = parse_folder_for_metadata(line)
    output.append(folder_metadata["total_electrons"])
    cosy_output = CosyOutput()
    cosy_output.injectScreenFile(os.path.join(line,"2-Cosy/screen.txt"))
    row_numbers = cosy_output.returnRowNumbersWithKeyValuePair("time",120e-12,comparison = 'gt')
    if row_numbers == []:
      continue
    output.append(str(100 * int(cosy_output.returnTableValue("number macroparticles",row_numbers[0]))))
    output.append(folder_metadata["applied_field"])
    print ",".join(output)
    
