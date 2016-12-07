import sys
import os
import argparse
import numpy as np
import math
import random
from scipy import constants
from cosy_output import CosyOutput
from coordinate_vector_3d import Cartesian3DVector
from phase_volume_6d import Phase6DVolume
from my_covariance_matrix import MyEditableCovarianceMatrix

parser = argparse.ArgumentParser(description='Write the columns for the emittance versus time.')
parser.add_argument('list_of_directories', type=str, help='The path to a where file containing the list of the directories to be used to calculated the stats..', default=None)
parser.add_argument('-s','--step', type=int, help='The timestep value for which to calculate the emittance.', default=100)
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('-n','--number_macroparticle', dest="number_of_particles", type=int, help='The number of macroparticles in the simulation.  This defaults to 5000 unless specified.', default=5000)



args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron


header = []
header.append("number_trials")
header.append("det_ref_cov_matrix")
header.append("det_cov_matrix")
header.append("log_emittance")
print ",".join(header)
step = args.step #Only calculate once all particles are present.
in_loop = True
with open(args.list_of_directories) as f:
  directories = [line.rstrip() for line in f]

trials_step_size = 10
last_trials_step = int(math.floor(len(directories)/trials_step_size))
for i in range(1,last_trials_step):
  number_trials = i*trials_step_size
  trial_dirs = random.sample(directories,number_trials)
  sum_covariance_matrix = MyEditableCovarianceMatrix()
  mean_phase_volume=Phase6DVolume(mass=mass_of_macroparticle)
  for directory in trial_dirs:
    cov_matrix = MyEditableCovarianceMatrix()
    cov_matrix.readCovarianceMatrix(os.path.join(directory,str(args.step)+"_cov_matrix.txt"))
    sum_covariance_matrix.cov_matrix += cov_matrix.cov_matrix
    mean_phase_volume.injectFile(os.path.join(directory,str(args.step)+"_mean.txt"),mass=mass_of_macroparticle)
  sum_covariance_matrix.cov_matrix = sum_covariance_matrix.cov_matrix/len(mean_phase_volume)
  #sum_covariance_matrix.printCovarianceMatrix()
  det_cov_matrix = sum_covariance_matrix.getSubDeterminant()
  ref_cov_matrix = mean_phase_volume.getCovarianceMatrix()
  det_ref_cov_matrix = ref_cov_matrix.getSubDeterminant()
  #ref_cov_matrix.printCovarianceMatrix()
  log_emittance = np.log10(det_ref_cov_matrix) + (args.number_of_particles-1)*np.log10(det_cov_matrix)
  output = []
  output.append(number_trials)
  output.append(det_ref_cov_matrix)
  output.append(det_cov_matrix)
  output.append(log_emittance)
  print ",".join([str(o) for o in output])
