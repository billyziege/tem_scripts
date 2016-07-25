import sys
import os
import math
import argparse
import numpy as np
from scipy import constants
from conventional_com_boost import ConventionalCOMBoost

class EVarCOMBoost(ConventionalCOMBoost):
  """
  A covariance matrix wiht the evar boost, v/c = c<p^lab>/<E^lab>, built into it.
  """

  def getCOMLorentzGamma(self,direction="z"):
    """
    Wraps the calculation of the lorentz gamma
    and returns it.
    """
    if not hasattr(self,"lorentz_gamma"):
      self.lorentz_gamma = {}
    if not direction in self.lortenz_gamma:
      v_over_c = self.getCOMVelocityOverC(direction)
      self.lorentz_gamma[direction] = np.sqrt(1/(1-v_over_c**2))
    return self.lorentz_gamma[direction]

  def getCOMVelocityOverC(self,direction):
    """
    Wraps the calculation of the v/c
    and returns it.
    """
    if not hasattr(self,"v_over_c"):
      self.v_over_c = {}
    if not direction in self.v_over_c:
      self.v_over_c[direction] = self.calcCOMVelocityOverC(direction)
    return self.v_over_c[direction]

  def calcCOMVelocityOverC(self,direction):
    """
    Returns v/c = c <p_direction>/<E>.
      direction must be one of "x", "y", or "z".
    """
    momentum_in_direction = []
    energy = []
    for particle in self:
      momentum_in_direction.append(getattr(particle.p,direction))
      energy.append(particle.getEnergy())
    mean_momentum = np.mean(momentum_in_direction)
    mean_energy = np.energy(energy)
    return mean_momentum/mean_energy

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package and use simple commands to get some of the straightforward methods.')
  parser.add_argument('coordinates_files', nargs='+', type=str, help='The path(s) to the phase volume files containing the phase space and therefore specifying that the injection function should be tested.')
  parser.add_argument('-c','--covariance_matrix', dest="covariance_matrix", action="store_true", help='Prints out the upper triangular form of the covariance matrix.', default=False)
  parser.add_argument('-r','--correlation_matrix', dest="correlation_matrix", action="store_true", help='Prints out the upper triangular form of the correlation matrix.', default=False)
  parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
  args = parser.parse_args()

  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

  phase_volume = Phase6DVolume()
  for path in args.coordinates_files:
    phase_volume.injectFile(path,mass=mass_of_macroparticle)

  if args.covariance_matrix:
    phase_volume.printCovarianceMatrix()

  if args.correlation_matrix:
    phase_volume.printCorrelationMatrix()

