import sys
import os
import argparse
import numpy as np
from scipy import constants
from cosy_output import CosyOutput
from coordinate_vector_3d import Cartesian3DVector
from phase_volume_6d import Phase6DVolume
from my_covariance_matrix import MyEditableCovarianceMatrix

parser = argparse.ArgumentParser(description='Write the columns for the emittance versus time.')
parser.add_argument('list_of_directories', type=str, help='The path to a where file containing the list of the directories to be used to calculated the stats..', default=None)
parser.add_argument('-s','--stepsize', type=int, help='The number by which *-x-*.txt were reported.', default=10)
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)



args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

scaling_factor_for_x_mean = 1e8
scaling_factor_for_p_mean = 1e4

header = []
header.append("time")
header.append("det_ref_cov_matrix")
header.append("det_cov_matrix")
header.append("log_emittance")
print ",".join(header)
step = 100 #Only calculate once all particles are present.
in_loop = True
while in_loop is True:
  output = []
  sum_covariance_matrix = MyEditableCovarianceMatrix()
  mean_phase_volume=Phase6DVolume(mass=mass_of_macroparticle)
  with open(args.list_of_directories) as f:
    initialized = False
    for directory in f:
      directory = directory.rstrip()
      if not initialized:
        cosy_output = CosyOutput()
        cosy_output.injectFile(os.path.join(directory,"screen.txt"))
        row_numbers = cosy_output.getRowNumbersWithKeyValuePair("step number",step,comparison = 'eq')
        if row_numbers == []:
          in_loop = False
          break
        time = float(cosy_output.returnTableValue("time",row_numbers[0]))
        output.append(time)
        initialized = True
      phase_volume = Phase6DVolume()
      for filename in os.listdir(directory):
        if filename.startswith(str(step)+"-x-") and filename.endswith(".txt"):
          phase_volume.injectFile(os.path.join(directory,filename),mass=mass_of_macroparticle)
      sum_covariance_matrix.cov_matrix += phase_volume.getCovarianceMatrix().cov_matrix
      x_mean = Cartesian3DVector(x=scaling_factor_for_x_mean*phase_volume.getMean(["x"]),y=scaling_factor_for_x_mean*phase_volume.getMean(["y"]),z=scaling_factor_for_x_mean*phase_volume.getMean(["z"]))
      #print x_mean
      p_mean = Cartesian3DVector(x=scaling_factor_for_p_mean*phase_volume.getMean(["px"]),y=scaling_factor_for_p_mean*phase_volume.getMean(["py"]),z=phase_volume.getMean(["pz"]))
      #print p_mean
      mean_phase_volume.addParticle(mean_phase_volume.mass,x=x_mean,p=p_mean)
  sum_covariance_matrix.cov_matrix = sum_covariance_matrix.cov_matrix/len(mean_phase_volume)
  #sum_covariance_matrix.printCovarianceMatrix()
  det_cov_matrix = sum_covariance_matrix.getSubDeterminant()
  ref_cov_matrix = mean_phase_volume.getCovarianceMatrix()
  det_ref_cov_matrix = ref_cov_matrix.getSubDeterminant()/(scaling_factor_for_x_mean**6*scaling_factor_for_p_mean**4)
  #ref_cov_matrix.printCovarianceMatrix()
  log_emittance = np.log10(det_ref_cov_matrix) + (len(phase_volume)-1)*np.log10(det_cov_matrix)
  output.append(det_ref_cov_matrix)
  output.append(det_cov_matrix)
  output.append(log_emittance)
  print ",".join([str(o) for o in output])
  step += args.stepsize
