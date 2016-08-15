import os
import argparse
import random
from scipy import constants
from phase_volume_6d import Phase6DVolume

parser = argparse.ArgumentParser(description='Write the interparticle distances for the given phase volume.')
parser.add_argument('phase_volume_files', nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space.')
parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
parser.add_argument('-s','--sample', dest='sample', type=int, help="The number of pairs to randomly sample. Default is 10000.", default=10000)
args = parser.parse_args()

mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

phase_volume = Phase6DVolume()
for filepath in args.phase_volume_files:
  filename = os.path.basename(filepath)
  phase_volume.injectFile(filepath,mass=mass_of_macroparticle)
for i in range(1,args.sample):
  particle_i = random.sample(phase_volume.particles,1).pop()#Randomly select a particle.
  particle_j = random.sample(phase_volume.particles,1).pop()
  dist = particle_i.x - particle_j.x
  if abs(dist) == 0:
    continue
  output = []
  output.append(abs(dist.x))
  output.append(abs(dist.y))
  output.append(abs(dist.z))
  output.append(abs(dist))
  print " ".join([str(o) for o in output])
