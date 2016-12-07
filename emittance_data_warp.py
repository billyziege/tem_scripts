import sys
import os
import argparse
import numpy as np
from scipy import constants
from folder_convention import parse_folder_for_metadata
from cosy_output import CosyOutput
from phase_volume_6d import Phase6DVolume
from com_boost import ConventionalCOMBoost, EVarCOMBoost, PTCOMBoost

parser = argparse.ArgumentParser(description='Write the columns for the emittance in time: time, z, std_z, boosted std_z, det COV, con exeyez, det con boosted COV, evar exeyez, det evar boosted COV, pt exeyez, det pt boosted COV')
parser.add_argument('list_of_directories', type=str, help='The path to a file containing all of directories in which the 2-Cosy/screen.txt are stored.', default=None)
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)

args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

print "Ne,ex,ez,Fa";
with open(args.list_of_directories) as f:
  for line in f:
    line = line.rstrip()
    output =[]
    folder_metadata = parse_folder_for_metadata(line)
    #cosy_output = CosyOutput()
    #cosy_output.injectFile(os.path.join(line,"2-Cosy/screen.txt"))
    #row_numbers = cosy_output.getRowNumbersWithKeyValuePair("time",85e-12,comparison = 'gt')
    #if row_numbers == []:
     # print cosy_output.rows
     # print "skipping"
     # exit()
     # continue
    #step = int(cosy_output.returnTableValue("step number",row_numbers[0]))
    step = 450 + int(np.floor(int(folder_metadata["steps"])/10) * 10)
    phase_volume = Phase6DVolume()
    for filename in os.listdir(line):
      if filename == str(step)+"-warp_uem.txt":
        phase_volume.injectFile(os.path.join(line,filename),mass=mass_of_macroparticle,fieldnames=["x","y","z","px","py","pz","vx","vy","vz"])

    normal_cov_matrix = phase_volume.getCovarianceMatrix()
    if len(normal_cov_matrix) != 49:
      continue
    ex = np.sqrt(normal_cov_matrix.getSubDeterminant(["x","px"]))/mass_of_macroparticle
    ez = np.sqrt(normal_cov_matrix.getSubDeterminant(["z","pz"]))/mass_of_macroparticle
    if ex == 0 or ez == 0:
      continue
    output.append(str(int(folder_metadata["e_per_particle"]) * int(folder_metadata["n_particles"])))
    output.append(str(ex))
    output.append(str(ez))
    output.append(folder_metadata["applied_field"])
    print ",".join(output)
