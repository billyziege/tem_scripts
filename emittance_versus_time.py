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
parser.add_argument('directory', type=str, help='The path to a where screen.txt and stepsize-x-*.txt are stored.', default=None)
parser.add_argument('-s','--stepsize', type=int, help='The number by which *-x-*.txt were reported.', default=10)
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('--subsample', dest="subsample", action="store_true", help='Only use the first file of particles instead of the entire ensemble.', default=False)


args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

cosy_output = CosyOutput()
cosy_output.injectFile(os.path.join(args.directory,"screen.txt"))
header = []
header.append("time")
header.append("z")
header.append("std_z")
header.append("det_COV")
header.append("con_std_z")
header.append("con_gamma_z")
header.append("con_eta_z")
header.append("con_exeyez")
header.append("con_det_COV")
header.append("evar_std_z")
header.append("evar_gamma_z")
header.append("evar_eta_z")
header.append("evar_std_px")
header.append("evar_std_x")
header.append("evar_gamma_x")
header.append("evar_eta_x")
header.append("evar_ex")
header.append("evar_ey")
header.append("evar_ez")
header.append("evar_exeyez")
header.append("evar_det_COV")
header.append("pt_std_z")
header.append("pt_gamma_z")
header.append("pt_eta_z")
header.append("pt_exeyez")
header.append("pt_det_COV")
print ",".join(header)
step = args.stepsize
while True:
  output = []
  row_numbers = cosy_output.getRowNumbersWithKeyValuePair("step number",step,comparison = 'eq')
  if row_numbers == []:
    break
  time = float(cosy_output.returnTableValue("time",row_numbers[0]))
  output.append(time)
  z = float(cosy_output.returnTableValue("z",row_numbers[0]))
  output.append(z)
  phase_volume = Phase6DVolume()
  for filename in os.listdir(args.directory):
    if filename.startswith(str(step)+"-x-") and filename.endswith(".txt"):
      if args.subsample:
        if filename != str(step)+"-x-1.txt":
          continue
      phase_volume.injectFile(os.path.join(args.directory,filename),mass=mass_of_macroparticle)


  normal_cov_matrix = phase_volume.getCovarianceMatrix()
  output.append(np.sqrt(normal_cov_matrix.getCovarianceElement("z","z")))
  output.append(np.sqrt(normal_cov_matrix.getSubDeterminant())/(phase_volume.particles[0].mass**3))
  con_boost = ConventionalCOMBoost(phase_volume)
  con_boosted_cov_matrix = con_boost.getBoostedCovarianceMatrix()
  con_var_z = con_boosted_cov_matrix.getCovarianceElement("z","z")
  output.append(np.sqrt(con_var_z))
  con_gamma_z = con_boosted_cov_matrix.getCovarianceElement("z","pz")
  output.append(con_gamma_z)
  con_var_pz = con_boosted_cov_matrix.getCovarianceElement("pz","pz")
  con_eta_z = con_var_pz - con_gamma_z**2/con_var_z
  output.append(np.sqrt(con_eta_z))
  
  con_ex = np.sqrt(con_boosted_cov_matrix.getSubDeterminant(["x","px"]))
  con_ey = np.sqrt(con_boosted_cov_matrix.getSubDeterminant(["y","py"]))
  con_ez = np.sqrt(con_boosted_cov_matrix.getSubDeterminant(["z","pz"]))
  con_exeyez = con_ex*con_ey*con_ez
  output.append(con_exeyez/(con_boost.mass**3))
  
  output.append(np.sqrt(con_boosted_cov_matrix.getSubDeterminant())/(con_boost.mass**3))
  
  evar_boost = EVarCOMBoost(phase_volume)
  evar_boost.phase_volume.getCovarianceMatrix(recalculate=True)
  evar_boosted_cov_matrix = evar_boost.getBoostedCovarianceMatrix()
  evar_var_z = evar_boosted_cov_matrix.getCovarianceElement("z","z")
  output.append(np.sqrt(evar_var_z))
  evar_gamma_z = evar_boosted_cov_matrix.getCovarianceElement("z","pz")
  output.append(evar_gamma_z)
  evar_var_pz = evar_boosted_cov_matrix.getCovarianceElement("pz","pz")
  evar_eta_z = evar_var_pz - evar_gamma_z**2/evar_var_z
  output.append(np.sqrt(evar_eta_z))
  evar_var_px = evar_boosted_cov_matrix.getCovarianceElement("px","px")
  output.append(np.sqrt(evar_var_px))
  evar_var_x = evar_boosted_cov_matrix.getCovarianceElement("x","x")
  output.append(np.sqrt(evar_var_x))
  evar_gamma_x = evar_boosted_cov_matrix.getCovarianceElement("x","px")
  output.append(evar_gamma_x)
  evar_eta_x = evar_var_px - evar_gamma_x**2/evar_var_x
  output.append(np.sqrt(evar_eta_x))
  
  evar_ex = np.sqrt(evar_boosted_cov_matrix.getSubDeterminant(["x","px"]))/evar_boost.mass
  evar_ey = np.sqrt(evar_boosted_cov_matrix.getSubDeterminant(["y","py"]))/evar_boost.mass
  evar_ez = np.sqrt(evar_boosted_cov_matrix.getSubDeterminant(["z","pz"]))/evar_boost.mass
  evar_exeyez = evar_ex*evar_ey*evar_ez
  output.append(evar_ex)
  output.append(evar_ey)
  output.append(evar_ez)
  output.append(evar_exeyez)

  output.append(np.sqrt(evar_boosted_cov_matrix.getSubDeterminant())/(evar_boost.mass**3))

  pt_boost = PTCOMBoost(phase_volume)
  pt_boosted_cov_matrix = pt_boost.getBoostedCovarianceMatrix()
  pt_var_z = pt_boosted_cov_matrix.getCovarianceElement("z","z")
  output.append(np.sqrt(pt_var_z))
  pt_gamma_z = pt_boosted_cov_matrix.getCovarianceElement("z","pz")
  output.append(pt_gamma_z)
  pt_var_pz = pt_boosted_cov_matrix.getCovarianceElement("pz","pz")
  pt_eta_z = pt_var_pz - pt_gamma_z**2/pt_var_z
  output.append(np.sqrt(pt_eta_z))
  
  pt_ex = np.sqrt(pt_boosted_cov_matrix.getSubDeterminant(["x","px"]))
  pt_ey = np.sqrt(pt_boosted_cov_matrix.getSubDeterminant(["y","py"]))
  pt_ez = np.sqrt(pt_boosted_cov_matrix.getSubDeterminant(["z","pz"]))
  pt_exeyez = pt_ex*pt_ey*pt_ez
  output.append(pt_exeyez/(pt_boost.mass**3))

  output.append(np.sqrt(pt_boosted_cov_matrix.getSubDeterminant())/(pt_boost.mass**3))

  print ",".join([str(o) for o in output])
  step += args.stepsize
