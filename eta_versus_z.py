import sys
import os
import numpy
import argparse
from scipy import constants
from folder_convention import parse_folder_for_metadata
from cosy_output import CosyOutput
from simulation_output import SimulationOutput

parser = argparse.ArgumentParser(description='Write the two columns for the velocity data: z,gamma_z')
parser.add_argument('directory', type=str, help='The path where the step files and screen.txt are stored.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
args = parser.parse_args()

print "z,eta_z,eta_z_2,ez,gamma_z,mean_z,mean_pz,mean_pzxz,mean_z_squared,mean_pz_squared,var_z,var_pz"#,e_spread,var_pz,gamma_z_squared,var_z";
mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron
cosy_output = CosyOutput()
cosy_output.injectFile(os.path.join(args.directory,"screen.txt"))
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
    simulation_output.convertCoordinatesToBetterUnits(mass_of_macroparticle)
    simulation_output = simulation_output.boostCoordinates("z")
    eta_z = simulation_output.calcEta("z")
    ez = simulation_output.calcEmittance("z")
    position = simulation_output.returnAllColumnValues("z")
    momentum = simulation_output.returnAllColumnValues("p_z")
    mean_position = numpy.mean(position)
    mean_momentum = numpy.mean(momentum)
    mean_positionxmomentum = numpy.mean([x*p for x,p in zip(position,momentum)])
    mean_position_squared = numpy.mean([x*x for x in position])
    mean_momentum_squared = numpy.mean([p*p for p in momentum])
    var_position = numpy.var(position,ddof=0)
    var_momentum = numpy.var(momentum,ddof=0)
    eta_z_2 = ez**2/var_position
    #lorentz_gamma, e_spread = simulation_output.calcEnergySpread("z")
    gamma_z = simulation_output.calcGamma("z")
    #momentum = simulation_output.returnAllColumnValues("p_z")
    output.append(str(numpy.sqrt(eta_z)*args.number_of_electrons_per_macroparticle))#to move it into the m of macroparticle instead of me
    output.append(str(numpy.sqrt(eta_z_2)*args.number_of_electrons_per_macroparticle))#to move it into the m of macroparticle instead of me
    output.append(str(ez))#to move it into the m of macroparticle instead of me
    output.append(str(gamma_z))
    output.append(str(mean_position))
    output.append(str(mean_momentum))
    output.append(str(mean_positionxmomentum))
    output.append(str(mean_position_squared))
    output.append(str(mean_momentum_squared))
    output.append(str(var_position));#to move it into the m of macroparticle instead of me
    output.append(str(var_momentum));#to move it into the m of macroparticle instead of me
    #output.append(str(e_spread))
    #output.append(str(var_momentum))
    #output.append(str(var_position))
    print ",".join(output)
