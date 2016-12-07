import sys
import os
import math
import argparse
import numpy as np
from scipy import constants
from scipy import special
from coordinate_vector import CoordinateException
from coordinate_vector_3d import Cartesian3DVector
from phase_volume_6d import Phase6DVolume

class DisretizedPhase6DVolume(Phase6DVolume):
  """
  An ensemble of 3D position and 3D momentum vectors
  with their 3D position within the min and max vectors.
  """

  def __init__(self,min_vector=None,max_vector=None):
    """
    Sets the min and max vectors (defaults are 0,0,0)
    and passes control to Phase6DVolume.
    """
    if min_vector is None:
      min_vector = Cartesian3DVector()
    if max_vector is None:
      max_vector = Cartesian3DVector()
    if not isinstance(min_vector,Cartesian3DVector):
      raise CoordinateException("The minimum vector for a binned phase space needs to be a Cartesian3DVector.")
    if not isinstance(max_vector,Cartesian3DVector):
      raise CoordinateException("The maximum vector for a binned phase space needs to be a Cartesian3DVector.")
    self.min = min_vector
    self.max = max_vector
    Phase6DVolume.__init__()

class BinnedPhase6DVolume():
  """
  An ensemble of 3D position and 3D momentum vectors
  organized by their extents in position (xyz) space.
  """

  def __init__(self,phase_volume,number_bins_1d=None,number_bins_x=None,number_bins_y=None,
               number_bins_z=None,x_mean=None,y_mean=None,z_mean=None):
    """
    Splits the given phase_volume into bins and stores it.
    """
    self.original_phase_volume = phase_volume
    self.phase_volumes = []
    self.mass = phase_volume.mass
    if number_bins_1d is not None:
      number_bins_x = number_bins_1d
      number_bins_y = number_bins_1d
      number_bins_z = number_bins_1d
    if number_bins_x is None or number_bins_y is None or number_bins_z is None:
      raise CoordinateException("Binning the phase space requires the bins in each direction be provided.")
    self.number_bins_x = number_bins_x
    self.number_bins_y = number_bins_y
    self.number_bins_z = number_bins_z
    self.x_max = phase_volume.getMax("x")
    self.y_max = phase_volume.getMax("y")
    self.z_max = phase_volume.getMax("z")
    self.x_min = phase_volume.getMin("x")
    self.y_min = phase_volume.getMin("y")
    self.z_min = phase_volume.getMin("z")
    if x_mean is not None:
      x_diff = self.x_max - x_mean
      x_diff_lower = x_mean - self.x_min
      if x_diff < x_diff_lower:
        x_diff = x_diff_lower
      self.x_max = x_mean + x_diff
      self.x_min = x_mean - x_diff
    if y_mean is not None:
      y_diff = self.y_max - y_mean
      y_diff_lower = y_mean - self.y_min
      if y_diff < y_diff_lower:
        y_diff = y_diff_lower
      self.y_max = y_mean + y_diff
      self.y_min = y_mean - y_diff
    if z_mean is not None:
      z_diff = self.z_max - z_mean
      z_diff_lower = z_mean - self.z_min
      if z_diff < z_diff_lower:
        z_diff = z_diff_lower
      self.z_max = z_mean + z_diff
      self.z_min = z_mean - z_diff
    self.x_len = self.x_max - self.x_min
    self.y_len = self.y_max - self.y_min
    self.z_len = self.z_max - self.z_min
    self.x_step = self.x_len/number_bins_x
    self.y_step = self.y_len/number_bins_y
    self.z_step = self.z_len/number_bins_z
    self.buildPhaseVolumesStructure()
    self.binOriginalPhaseVolume()

  def buildPhaseVolumesStructure(self):
    """"
    Builds the array structure of the phase volumes attribute.
    """
    for i in range(0,self.number_bins_x):
      self.phase_volumes.append([])
      for j in range(0,self.number_bins_y):
        self.phase_volumes[i].append([])
        for k in range(0,self.number_bins_z):
          self.phase_volumes[i][j].append(Phase6DVolume())

  def binOriginalPhaseVolume(self):
    """
    Distributes the particles from the original phase volume through the
    binned structure.
    """
    for particle in self.original_phase_volume:
      ix = int(np.floor((particle.x.x - self.x_min)/self.x_step))
      iy = int(np.floor((particle.x.y - self.y_min)/self.y_step))
      iz = int(np.floor((particle.x.z - self.z_min)/self.z_step))
      if ix == self.number_bins_x:
        ix -= 1
      if iy == self.number_bins_y:
        iy -= 1
      if iz == self.number_bins_z:
        iz -= 1
      self.phase_volumes[ix][iy][iz].addParticle(particle)

  def getAverageTemperature(self,direction):
    """
    Checks to see if the average temperature has already been caclulated.  If so, returns
    the previously identified temperature.  Otherwise, calculates the temperature and stores
    if for later retrieval and returning.  Note, direction is "x", "y" or "z".
    """
    if not hasattr(self,"average_temperature"):
      self.average_temperature = {}
    if not direction in self.average_temperature:
      self.average_temperature[direction] = self.calcAverageTemperature(direction)
    return self.average_temperature[direction]

  def calcAverageTemperature(self,direction):
    """
    This finds the average temperature across the grid points' ensembles
    of particles.  
    """
    temperatures = self.getTemperatures()
    average_temperature = 0
    total_particles = 0
    for i in range(0,len(temperatures[direction])):
      for j in range(0,len(temperatures[direction][i])):
        for k in range(0,len(temperatures[direction][i][j])):
          n = len(self.phase_volumes[i][j][k])
          if n > 1:
            total_particles += n
            average_temperature += temperatures[direction][i][j][k]*n
    return average_temperature/total_particles

  def getAverageDebyeLengths(self):
    """
    Checks to see if the average Debye length has already been caclulated.  If so, return
    the previously identified value.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_debye_lengths"):
      self.average_debye_lengths = self.calcAverageDebyeLengths()
    return self.average_debye_lengths

  def calcAverageDebyeLengths(self):
    """
    This finds the debye lengths.
    """
    coefficient = np.sqrt(self.epsilon/self.charge) #Charge in square root because of temp in eV!

    average_density = self.original_phase_volume.getAverageDensity()
    average_temperature_x = self.getAverageTemperature("x")
    average_temperature_y = self.getAverageTemperature("y")
    average_temperature_perp = 0.5*(average_temperature_x+average_temperature_y)
    average_temperature_z = self.getAverageTemperature("z")
    average_temperature = (average_temperature_x+average_temperature_y+average_temperature_z)/3.0

    average_debye_lengths = {}
    average_debye_lengths["x"] = coefficient*np.sqrt(average_temperature_x/average_density)
    average_debye_lengths["y"] = coefficient*np.sqrt(average_temperature_x/average_density)
    average_debye_lengths["perp"] = coefficient*np.sqrt(average_temperature_perp/average_density)
    average_debye_lengths["z"] = coefficient*np.sqrt(average_temperature_z/average_density)
    average_debye_lengths["total"] = coefficient*np.sqrt(average_temperature/average_density)
    
    return average_debye_lengths

  def getAveragePlasmaFrequency(self):
    """
    Checks to see if the average plasma frequency has already been caclulated.  If so, return
    the previously identified value.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_plasma_frequency"):
      self.average_plasma_frequency = self.calcAveragePlasmaFrequency()
    return self.average_plasma_frequency

  def calcAveragePlasmaFrequency(self):
    """
    This finds the plasma frequency.
    """
    average_density = self.original_phase_volume.getAverageDensity()
    return self.charge*np.sqrt(average_density/(self.epsilon * self.mass_in_kg))

  def getAverageCoulombParameter(self):
    """
    Checks to see if the coulomb parameter has already been caclulated.  If so, return
    the parameter.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_colulomb_parameter"):
      self.average_coulomb_parameter = self.calcAverageCoulombParameter()
    return self.average_coulomb_parameter

  def calcAverageCoulombParameter(self):
    """
    This finds the average coulomb parameter.
    """
    pi = 3.14159
    coefficient = 4*pi/3 

    average_density = self.original_phase_volume.getAverageDensity()
    average_debye_lengths = self.getAverageDebyeLengths()

    average_coulomb_parameter = {}
    average_coulomb_parameter["total"] = coefficient*average_density*average_debye_lengths["total"]**3
    average_coulomb_parameter["elipsoid"] = coefficient*average_density*average_debye_lengths["perp"]**2*average_debye_lengths["z"]
    average_coulomb_parameter["perp"] = coefficient*average_density*average_debye_lengths["perp"]**3
    average_coulomb_parameter["z"] = coefficient*average_density*average_debye_lengths["z"]**3
    return average_coulomb_parameter

  def getAverageCouplingParameter(self):
    """
    Checks to see if the coupling parameter has already been caclulated.  If so, return
    the parameter.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_colulomb_parameter"):
      self.average_coupling_parameter = self.calcAverageCouplingParameter()
    return self.average_coupling_parameter

  def calcAverageCouplingParameter(self):
    """
    This finds the average coupling parameter.
    """
    pi = 3.14159
    coefficient = self.charge/(4*pi*self.epsilon) #Charge has power one because of temp in eV!

    average_density = self.original_phase_volume.getAverageDensity()
    tx = self.getAverageTemperature('x')
    ty = self.getAverageTemperature('y')
    tz = self.getAverageTemperature('z')
    tperp = 0.5*(tx+ty)
    ttotal = 1.0/3.0*(tx+ty+tz)

    average_coupling_parameter = {}
    average_coupling_parameter["total"] = coefficient*special.cbrt(average_density)/ttotal
    average_coupling_parameter["perp"] = coefficient*special.cbrt(average_density)/tperp
    average_coupling_parameter["z"] = coefficient*special.cbrt(average_density)/tz
    return average_coupling_parameter

  def getTemperatures(self):
    """
    Checks to see if the temperatures have already been caclulated.  If so, return
    the previously identified temperatre.  Otherwise, calculate the temperatures and store
    if for later retrieval and returning.
    """
    if not hasattr(self,"temperatures"):
      self.temperatures = self.calcTemperatures()
    return self.temperatures

  def calcTemperatures(self):
    """
    This finds the temperatures across the grid points' ensembles
    of particles.  
    """
    temperature = {}
    temperature["x"] = []
    temperature["y"] = []
    temperature["z"] = []
    denominator = self.mass*self.k_b
    for i in range(0,len(self.phase_volumes)):
      temperature["x"].append([])
      temperature["y"].append([])
      temperature["z"].append([])
      for j in range(0,len(self.phase_volumes[i])):
        temperature["x"][i].append([])
        temperature["y"][i].append([])
        temperature["z"][i].append([])
        for k in range(0,len(self.phase_volumes[i][j])):
          if len(self.phase_volumes[i][j][k]) <= 1:
            temperature["x"][i][j].append(-1)
            temperature["y"][i][j].append(-1)
            temperature["z"][i][j].append(-1)
            continue
          cov_matrix = self.phase_volumes[i][j][k].getCovarianceMatrix()
          temperature["x"][i][j].append(cov_matrix.getCovarianceElement("px","px")/denominator)
          temperature["y"][i][j].append(cov_matrix.getCovarianceElement("py","py")/denominator)
          temperature["z"][i][j].append(cov_matrix.getCovarianceElement("pz","pz")/denominator)
    return temperature

  def getDebyeLengths(self):
    """
    Checks to see if the Debye lengths have already been caclulated.  If so, return
    the previously identified Debye lengths.  Otherwise, calculate it.
    """
    if not hasattr(self,"debye_lengths"):
      self.debye_lengths = self.calcDebyeLengths()
    return self.debye_lengths

  def calcDebyeLengths(self):
    """
    This finds the Debye lengths across the grid spaces.
    """
    debye_lengths = {}
    debye_lengths["x"] = []
    debye_lengths["y"] = []
    debye_lengths["z"] = []
    coefficient = np.sqrt(self.epsilon/self.charge) #Charge in square root because of temp in eV!
    temperatures = self.getTemperatures()
    volume = self.z_step * self.y_step * self.x_step
    for i in range(0,len(temperatures["x"])):
      debye_lengths["x"].append([])
      debye_lengths["y"].append([])
      debye_lengths["z"].append([])
      for j in range(0,len(temperatures["x"][i])):
        debye_lengths["x"][i].append([])
        debye_lengths["y"][i].append([])
        debye_lengths["z"][i].append([])
        for k in range(0,len(temperatures["x"][i][j])):
          n = len(self.phase_volumes[i][j][k])
          if n <= 0 or temperatures["x"][i][j][k] == -1:
            debye_lengths["x"][i][j].append(-1)
            debye_lengths["y"][i][j].append(-1)
            debye_lengths["z"][i][j].append(-1)
            continue
          density = n*self.electrons_per_macroparticle/volume
          debye_lengths["x"][i][j].append(coefficient*np.sqrt(temperatures["x"][i][j][k]/density))
          debye_lengths["y"][i][j].append(coefficient*np.sqrt(temperatures["y"][i][j][k]/density))
          debye_lengths["z"][i][j].append(coefficient*np.sqrt(temperatures["z"][i][j][k]/density))
    return debye_lengths

  def getPlasmaFrequencies(self):
    """
    Checks to see if the plasma frequencies have already been caclulated.  If so, return
    the previously identified Debye lengths.  Otherwise, calculate it.
    """
    if not hasattr(self,"plasma_frequencies"):
      self.plasma_frequencies = self.calcPlasmaFrequencies()
    return self.plasma_frequencies

  def calcPlasmaFrequencies(self):
    """
    This finds the plasma frequencies across the grid spaces.
    """
    plasma_frequencies = []
    coefficient = self.charge/np.sqrt(self.epsilon * self.mass_in_kg)
    volume = self.z_step * self.y_step * self.x_step
    for i in range(0,len(self.phase_volumes)):
      plasma_frequencies.append([])
      for j in range(0,len(self.phase_volumes[i])):
        plasma_frequencies[i].append([])
        for k in range(0,len(self.phase_volumes[i][j])):
          n = len(self.phase_volumes[i][j][k])
          if n <= 0:
            plasma_frequencies[i][j].append(-1)
            continue
          density = n*self.electrons_per_macroparticle/volume
          plasma_frequencies[i][j].append(coefficient*np.sqrt(density))
    return plasma_frequencies

  def getCoulombParameters(self):
    """
    Checks to see if the coulomb parameter has already been caclulated.  If so, return
    the parameter.  Otherwise, calculate it.
    """
    if not hasattr(self,"colulomb_parameters"):
      self.coulomb_parameters = self.calcCoulombParameters()
    return self.coulomb_parameters

  def calcCoulombParameters(self):
    """
    This finds the local coulomb parameters.
    """
    coulomb_parameters = {}
    coulomb_parameters["total"] = []
    coulomb_parameters["perp"] = []
    coulomb_parameters["z"] = []
    pi = 3.14159
    coefficient = 4*pi/3 
    debye_lengths = self.getDebyeLengths()

    volume = self.z_step * self.y_step * self.x_step
    for i in range(0,len(self.phase_volumes)):
      coulomb_parameters["total"].append([])
      coulomb_parameters["perp"].append([])
      coulomb_parameters["z"].append([])
      for j in range(0,len(self.phase_volumes[i])):
        coulomb_parameters["total"][i].append([])
        coulomb_parameters["perp"][i].append([])
        coulomb_parameters["z"][i].append([])
        for k in range(0,len(self.phase_volumes[i][j])):
          n = len(self.phase_volumes[i][j][k])
          if n <= 0:
            coulomb_parameters["total"][i][j].append(-1)
            coulomb_parameters["perp"][i][j].append(-1)
            coulomb_parameters["z"][i][j].append(-1)
            continue
          density = n*self.electrons_per_macroparticle/volume
          coulomb_parameters["total"][i][j].append(coefficient*density
                                                     *debye_lengths["x"][i][j][k]
                                                     *debye_lengths["y"][i][j][k]
                                                     *debye_lengths["z"][i][j][k])
          debye_length_perp = 0.5*(debye_lengths["x"][i][j][k] + debye_lengths["y"][i][j][k])
          coulomb_parameters["perp"][i][j].append(coefficient*density*debye_length_perp**3)
          coulomb_parameters["z"][i][j].append(coefficient*density*debye_lengths["z"][i][j][k]**3)
    return coulomb_parameters

  def getCouplingParameters(self):
    """
    Checks to see if the coupling parameter has already been caclulated.  If so, return
    the parameter.  Otherwise, calculate it.
    """
    if not hasattr(self,"coupling_parameters"):
      self.coupling_parameters = self.calcCouplingParameters()
    return self.coupling_parameters

  def calcCouplingParameters(self):
    """
    This finds the local coulomb parameters.
    """
    coupling_parameters = {}
    coupling_parameters["total"] = []
    coupling_parameters["perp"] = []
    coupling_parameters["z"] = []
    pi = 3.14159
    coefficient = self.charge/(4*pi*self.epsilon) #Charge has power one because of temp in eV!
    temperatures = self.getTemperatures()

    volume = self.z_step * self.y_step * self.x_step
    for i in range(0,len(self.phase_volumes)):
      coupling_parameters["total"].append([])
      coupling_parameters["perp"].append([])
      coupling_parameters["z"].append([])
      for j in range(0,len(self.phase_volumes[i])):
        coupling_parameters["total"][i].append([])
        coupling_parameters["perp"][i].append([])
        coupling_parameters["z"][i].append([])
        for k in range(0,len(self.phase_volumes[i][j])):
          n = len(self.phase_volumes[i][j][k])
          if n <= 0:
            coupling_parameters["total"][i][j].append(-1)
            coupling_parameters["perp"][i][j].append(-1)
            coupling_parameters["z"][i][j].append(-1)
            continue
          density = n*self.electrons_per_macroparticle/volume
          temperature_total = 1.0/3.0*( temperatures["x"][i][j][k]
                                       + temperatures["y"][i][j][k]
                                       + temperatures["z"][i][j][k] )
          if temperature_total > 0:
            coupling_parameters["total"][i][j].append(coefficient*special.cbrt(density)/temperature_total)
          else:
            coupling_parameters["total"][i][j].append(-1)
          temperature_perp = 0.5*(temperatures["x"][i][j][k] + temperatures["y"][i][j][k])
          if temperature_perp > 0:
            coupling_parameters["perp"][i][j].append(coefficient*special.cbrt(density)/temperature_perp)
          else:
            coupling_parameters["perp"][i][j].append(-1)
          if temperatures["z"][i][j][k] > 0:
            coupling_parameters["z"][i][j].append(coefficient*special.cbrt(density)/temperatures["z"][i][j][k])
          else:
            coupling_parameters["z"][i][j].append(-1)
    return coupling_parameters

  def getEmittances(self):
    """
    Checks to see if the emittances have already been caclulated.  If so, return
    the previously identified temperatre.  Otherwise, calculate the temperatures and store
    if for later retrieval and returning.
    """
    if not hasattr(self,"emittances"):
      self.emittances = self.calcEmittances()
    return self.emittances

  def calcEmittances(self):
    """
    This finds the rms normalized emittances across the grid points' ensembles
    of particles.  
    """
    emittances = {}
    emittances["x"] = []
    emittances["y"] = []
    emittances["z"] = []
    for i in range(0,len(self.phase_volumes)):
      emittances["x"].append([])
      emittances["y"].append([])
      emittances["z"].append([])
      for j in range(0,len(self.phase_volumes[i])):
        emittances["x"][i].append([])
        emittances["y"][i].append([])
        emittances["z"][i].append([])
        for k in range(0,len(self.phase_volumes[i][j])):
          if len(self.phase_volumes[i][j][k]) <= 2:
            emittances["x"][i][j].append(-1)
            emittances["y"][i][j].append(-1)
            emittances["z"][i][j].append(-1)
            continue
          cov_matrix = self.phase_volumes[i][j][k].getCovarianceMatrix()
          det_cov_x_px = cov_matrix.getSubDeterminant(subelements=["x","px"])
          det_cov_y_py = cov_matrix.getSubDeterminant(subelements=["y","py"])
          det_cov_z_pz = cov_matrix.getSubDeterminant(subelements=["z","pz"])
          if det_cov_x_px > 0:
            emittances["x"][i][j].append(np.sqrt(det_cov_x_px/self.mass))
          else:
            emittances["x"][i][j].append(0)
          if det_cov_y_py > 0:
            emittances["y"][i][j].append(np.sqrt(det_cov_y_py/self.mass))
          else:
            emittances["y"][i][j].append(0)
          if det_cov_z_pz > 0:
            emittances["z"][i][j].append(np.sqrt(det_cov_z_pz/self.mass))
          else:
            emittances["z"][i][j].append(0)
    return emittances

  def printTemperatures(self):
    """
    Prints to stdout the x, y, z, Tx, Ty, Tz data.
    """
    temperatures = self.getTemperatures()
    print "x,y,z,Tx,Ty,Tz,n,density"
    volume = self.z_step * self.y_step * self.x_step
    for i in range(0,len(temperatures["x"])):
      x = self.x_min + (i+0.5)*self.x_step
      for j in range(0,len(temperatures["x"][i])):
        y = self.y_min + (j+0.5)*self.y_step
        for k in range(0,len(temperatures["x"][i][j])):
          z = self.z_min + (k+0.5)*self.z_step
          n = len(self.phase_volumes[i][j][k])
          density = n*self.electrons_per_macroparticle/volume
          output = []
          output.append(x)
          output.append(y)
          output.append(z)
          output.append(temperatures["x"][i][j][k])
          output.append(temperatures["y"][i][j][k])
          output.append(temperatures["z"][i][j][k])
          output.append(n)
          output.append(density)
          print ",".join([str(o) for o in output])

  def printDebyeLengths(self):
    """
    Prints to stdout the x,y,z,Lx,Ly,Lz,n,density
    """
    debye_lengths = self.getDebyeLengths()
    print "x,y,z,Lx,Ly,Lz,n,density"
    volume = self.z_step * self.y_step * self.x_step
    for i in range(0,len(debye_lengths["x"])):
      x = self.x_min + (i+0.5)*self.x_step
      for j in range(0,len(debye_lengths["x"][i])):
        y = self.y_min + (j+0.5)*self.y_step
        for k in range(0,len(debye_lengths["x"][i][j])):
          z = self.z_min + (k+0.5)*self.z_step
          n = len(self.phase_volumes[i][j][k])
          density = n*self.electrons_per_macroparticle/volume
          output = []
          output.append(x)
          output.append(y)
          output.append(z)
          output.append(debye_lengths["x"][i][j][k])
          output.append(debye_lengths["y"][i][j][k])
          output.append(debye_lengths["z"][i][j][k])
          output.append(n)
          output.append(density)
          print ",".join([str(o) for o in output])

  def printBinnedStats(self):
    """
    Prints to stdout the x, y, z, Tx, Ty, Tz data.
    """
    pi = 3.14159
    temperatures = self.getTemperatures()
    debye_lengths = self.getDebyeLengths()
    emittances = self.getEmittances()
    plasma_frequencies = self.getPlasmaFrequencies()
    coulomb_parameters = self.getCoulombParameters()
    coupling_parameters = self.getCouplingParameters()
    print "x,y,z,Tx,Ty,Tz,Lx,Ly,Lz,ex,ey,ez,Ctotal,Cperp,Cz,Gtotal,Gperp,Gz,w,plasma_period,n,density"
    volume = self.z_step * self.y_step * self.x_step
    for i in range(0,len(temperatures["x"])):
      x = self.x_min + (i+0.5)*self.x_step
      for j in range(0,len(temperatures["x"][i])):
        y = self.y_min + (j+0.5)*self.y_step
        for k in range(0,len(temperatures["x"][i][j])):
          z = self.z_min + (k+0.5)*self.z_step
          n = len(self.phase_volumes[i][j][k])
          density = n*self.electrons_per_macroparticle/volume
          output = []
          output.append(x)
          output.append(y)
          output.append(z)
          if temperatures["x"][i][j][k] == -1:
            continue
          output.append(temperatures["x"][i][j][k])
          output.append(temperatures["y"][i][j][k])
          output.append(temperatures["z"][i][j][k])
          output.append(debye_lengths["x"][i][j][k])
          output.append(debye_lengths["y"][i][j][k])
          output.append(debye_lengths["z"][i][j][k])
          output.append(emittances["x"][i][j][k])
          output.append(emittances["y"][i][j][k])
          output.append(emittances["z"][i][j][k])
          output.append(coulomb_parameters["total"][i][j][k])
          output.append(coulomb_parameters["perp"][i][j][k])
          output.append(coulomb_parameters["z"][i][j][k])
          output.append(coupling_parameters["total"][i][j][k])
          output.append(coupling_parameters["perp"][i][j][k])
          output.append(coupling_parameters["z"][i][j][k])
          output.append(plasma_frequencies[i][j][k])
          output.append(2.0*pi/plasma_frequencies[i][j][k])
          output.append(n)
          output.append(density)
          print ",".join([str(o) for o in output])

  def printAverageStatistics(self,header=False):
    """
    Prints out the average statistics for the phase volume.
    """
    tx = binned_phase_volume.getAverageTemperature('x')
    ty = binned_phase_volume.getAverageTemperature('y')
    tz = binned_phase_volume.getAverageTemperature('z')
    output = [tx, ty, 0.5*(tx+ty),tz,1.0/3.0*(tx+ty+tz)]
    output.append(binned_phase_volume.original_phase_volume.getAverageDensity())
    w = binned_phase_volume.getAveragePlasmaFrequency()
    output.append(w)
    output.append(2.0*3.14159/w)
    average_debye_lengths = binned_phase_volume.getAverageDebyeLengths()
    output.append(average_debye_lengths["x"])
    output.append(average_debye_lengths["y"])
    output.append(average_debye_lengths["perp"])
    output.append(average_debye_lengths["z"])
    output.append(average_debye_lengths["total"])
    average_coulomb_parameter = binned_phase_volume.getAverageCoulombParameter()
    output.append(average_coulomb_parameter["total"])
    output.append(average_coulomb_parameter["elipsoid"])
    output.append(average_coulomb_parameter["perp"])
    output.append(average_coulomb_parameter["z"])
    average_coupling_parameter = binned_phase_volume.getAverageCouplingParameter()
    output.append(average_coupling_parameter["total"])
    output.append(average_coupling_parameter["perp"])
    output.append(average_coupling_parameter["z"])
    if header:
      print "tx,ty,tperp,tz,ttotal,density,plasma_freq,plasma_period,lx,ly,lperp,lz,ltotal,Coulombtotal,Coulombelipsoid,Coulombperp,Coulombz,Gammatotal,Gammaperp,Gammaz"
    print ",".join([str(o) for o in output])

class CylindricalBinnedPhase6DVolume():
  """
  An ensemble of 3D position and 3D momentum vectors
  organized by their extents in position (xyz) space.
  """

  def __init__(self,phase_volume,number_bins_1d=None,number_bins_rho=None,number_bins_z=None,rho_min=0,z_mean=None):
    """
    Splits the given phase_volume into bins and stores it.
    """
    self.orginal_phase_volume = phase_volume
    self.phase_volumes = []
    self.mass = phase_volume.mass
    if number_bins_1d is not None:
      number_bins_rho = number_bins_1d
      number_bins_z = number_bins_1d
    if number_bins_rho is None or number_bins_z is None:
      raise CoordinateException("Binning the phase space requires the bins in each direction be provided.")
    self.number_bins_rho = number_bins_rho
    self.number_bins_z = number_bins_z
    self.rho_max = phase_volume.getMaxR()
    self.z_max = phase_volume.getMax("z")
    self.rho_min = rho_min
    self.z_min = phase_volume.getMin("z")
    self.rho_len = self.rho_max - self.rho_min
    self.z_len = self.z_max - self.z_min
    self.rho_step = self.rho_len/number_bins_rho
    self.z_step = self.z_len/number_bins_z
    for i in range(0,number_bins_rho):
      self.phase_volumes.append([])
      for k in range(0,number_bins_z):
        self.phase_volumes[i].append(Phase6DVolume())
    for particle in phase_volume:
      coordinate = particle.x.convertToCylindrical()
      irho = int(np.floor((coordinate.rho - self.rho_min)/self.rho_step))
      iz = int(np.floor((particle.x.z - self.z_min)/self.z_step))
      if irho == number_bins_rho:
        irho -= 1
      if iz == number_bins_z:
        iz -= 1
      self.phase_volumes[irho][iz].addParticle(particle)

  def getTemperatures(self):
    """
    Checks to see if the temperatures have already been caclulated.  If so, return
    the previously identified temperatre.  Otherwise, calculate the temperatures and store
    if for later retrieval and returning.
    """
    if not hasattr(self,"temperatures"):
      self.temperatures = self.calcTemperatures()
    return self.temperatures

  def calcTemperatures(self):
    """
    This finds the temperatures across the grid points' ensembles
    of particles.  
    """
    temperature = {}
    temperature["rho"] = []
    temperature["phi"] = []
    temperature["z"] = []
    denominator = self.mass*self.k_b
    for i in range(0,len(self.phase_volumes)):
      temperature["rho"].append([])
      temperature["phi"].append([])
      temperature["z"].append([])
      for k in range(0,len(self.phase_volumes[i])):
        temperature["rho"][i].append([])
        temperature["phi"][i].append([])
        temperature["z"][i].append([])
        if len(self.phase_volumes[i][k]) <= 1:
          temperature["rho"][i][k] = -1
          temperature["phi"][i][k] = -1
          temperature["z"][i][k] = -1
          continue
        cov_matrix = self.phase_volumes[i][k].getCylindricalCovarianceMatrix()
        temperature["rho"][i][k] = cov_matrix.getCovarianceElement("prho","prho")/denominator
        temperature["phi"][i][k] = cov_matrix.getCovarianceElement("pphi","pphi")/denominator
        temperature["z"][i][k] = cov_matrix.getCovarianceElement("pz","pz")/denominator
    return temperature

  def getDebyeLengths(self):
    """
    Checks to see if the Debye lengths have already been caclulated.  If so, return
    the previously identified Debye lengths.  Otherwise, calculate it.
    """
    if not hasattr(self,"debye_lengths"):
      self.debye_lengths = self.calcDebyeLengths()
    return self.debye_lengths

  def calcDebyeLengths(self):
    """
    This finds the Debye lengths across the grid spaces.
    """
    pi=3.14159
    debye_lengths = {}
    debye_lengths["rho"] = []
    debye_lengths["phi"] = []
    debye_lengths["z"] = []
    coefficient = np.sqrt(self.epsilon)/self.charge
    temperatures = self.getTemperatures()
    for i in range(0,len(temperatures["rho"])):
      rho = self.rho_min + i*self.rho_step
      next_rho = self.rho_min + (i+1)*self.rho_step
      area = pi*(next_rho**2 - rho**2)
      volume = area*self.z_step
      debye_lengths["rho"].append([])
      debye_lengths["phi"].append([])
      debye_lengths["z"].append([])
      for k in range(0,len(temperatures["rho"][i])):
        n = len(self.phase_volumes[i][k])
        if n == 0:
          debye_lengths["rho"][i].append(-1)
          debye_lengths["phi"][i].append(-1)
          debye_lengths["z"][i].append(-1)
          continue
        if temperatures["rho"][i][k] == -1:
          debye_lengths["rho"][i].append(-1)
          debye_lengths["phi"][i].append(-1)
          debye_lengths["z"][i].append(-1)
          continue
        density = n*self.electrons_per_macroparticle/volume
        debye_lengths["rho"][i].append(coefficient*np.sqrt(temperatures["rho"][i][k]/density))
        debye_lengths["phi"][i].append(coefficient*np.sqrt(temperatures["phi"][i][k]/density))
        debye_lengths["z"][i].append(coefficient*np.sqrt(temperatures["z"][i][k]/density))
    return debye_lengths

  def printTemperatures(self):
    """
    Prints to stdout the rho, z, Trho, Tz data.
    """
    pi=3.14159
    temperatures = self.getTemperatures()
    print "rho,z,Trho,Tphi,Tz,n,density"
    for i in range(0,len(temperatures["rho"])):
      rho = self.rho_min + i*self.rho_step
      next_rho = self.rho_min + (i+1)*self.rho_step
      area = pi*(next_rho**2 - rho**2)
      volume = area*self.z_step
      for k in range(0,len(temperatures["rho"][i])):
        z = self.z_min + (k+0.5)*self.z_step
        n = len(self.phase_volumes[i][k])
        density = n*self.electrons_per_macroparticle/volume
        output = []
        output.append(rho)
        output.append(z)
        output.append(temperatures["rho"][i][k])
        output.append(temperatures["phi"][i][k])
        output.append(temperatures["z"][i][k])
        output.append(n)
        output.append(density)
        print ",".join([str(o) for o in output])

  def printDebyeLengths(self):
    """
    Prints to stdout the rho, z, Lrho, Lz
    """
    debye_lengths = self.getDebyeLengths()
    print "rho,z,Lrho,Lphi,Lz,n"
    for i in range(0,len(debye_lengths["rho"])):
      rho = self.rho_min + i*self.rho_step
      for k in range(0,len(debye_lengths["rho"][i])):
        z = self.z_min + (k+0.5)*self.z_step
        n = len(self.phase_volumes[i][k])
        output = []
        output.append(rho)
        output.append(z)
        output.append(debye_lengths["rho"][i][k])
        output.append(debye_lengths["phi"][i][k])
        output.append(debye_lengths["z"][i][k])
        output.append(n)
        print ",".join([str(o) for o in output])

class SlicedPhase6DVolume():
  """
  An ensemble of 3D position and 3D momentum vectors
  organized by their extents in the z direction.  This only 
  provides statistics in the z direction.
  """

  def __init__(self,phase_volume,number_bins_z=None,z_mean=None):
    """
    Splits the given phase_volume into bins and stores it.
    """
    self.original_phase_volume = phase_volume
    self.phase_volumes = []
    self.mass = phase_volume.mass
    if number_bins_z is None:
      raise CoordinateException("Binning the phase space requires the bins in each direction be provided.")
    self.number_bins_z = number_bins_z
    self.z_max = phase_volume.getMax("z")
    self.z_min = phase_volume.getMin("z")
    if z_mean is not None:
      z_diff = self.z_max - z_mean
      z_diff_lower = z_mean - self.z_min
      if z_diff < z_diff_lower:
        z_diff = z_diff_lower
      self.z_max = z_mean + z_diff
      self.z_min = z_mean - z_diff
    self.z_len = self.z_max - self.z_min
    self.z_step = self.z_len/number_bins_z
    for k in range(0,number_bins_z):
      self.phase_volumes.append(Phase6DVolume())
    for particle in phase_volume:
      iz = int(np.floor((particle.x.z - self.z_min)/self.z_step))
      if iz == number_bins_z:
        iz -= 1
      self.phase_volumes[iz].addParticle(particle)

  def getAverageTemperature(self,direction):
    """
    Checks to see if the average temperature has already been caclulated.  If so, returns
    the previously identified temperature.  Otherwise, calculates the temperature and stores
    if for later retrieval and returning.  
    """
    if not hasattr(self,"average_temperature"):
      self.average_temperature = {}
    if not direction in self.average_temperature:
      self.average_temperature[direction] = self.calcAverageTemperature(direction)
    return self.average_temperature[direction]

  def calcAverageTemperature(self,direction):
    """
    This finds the average temperature across the grid points' ensembles
    of particles.  
    """
    temperatures = self.getTemperatures()
    average_temperature = 0
    total_particles = 0
    for k in range(0,len(temperatures[direction])):
      n = len(self.phase_volumes[k])
      if n > 1:
        total_particles += n
        average_temperature += temperatures[direction][k]*n
    return average_temperature/total_particles

  def getAverageDebyeLengths(self,direction):
    """
    Checks to see if the average Debye length has already been caclulated.  If so, return
    the previously identified value.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_debye_lengths"):
      self.average_debye_lengths = {}
    if not direction in self.average_debye_lengths:
      self.average_debye_lengths[direction] = self.calcAverageDebyeLengths(direction)
    return self.average_debye_lengths[direction]

  def calcAverageDebyeLengths(self,direction):
    """
    This finds the debye lengths.
    """
    coefficient = np.sqrt(self.epsilon/self.charge) #Charge in square root because of temp in eV!

    average_density = self.original_phase_volume.getAverageDensity()

    average_debye_length = coefficient*np.sqrt(self.getAverageTemperature(direction)/average_density)
    
    return average_debye_length

  def getAveragePlasmaFrequency(self):
    """
    Checks to see if the average plasma frequency has already been caclulated.  If so, return
    the previously identified value.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_plasma_frequency"):
      self.average_plasma_frequency = self.calcAveragePlasmaFrequency()
    return self.average_plasma_frequency

  def calcAveragePlasmaFrequency(self):
    """
    This finds the plasma frequency.
    """
    average_density = self.original_phase_volume.getAverageDensity()
    return self.charge*np.sqrt(average_density/(self.epsilon * self.mass_in_kg))

  def getAverageCoulombParameter(self,direction):
    """
    Checks to see if the coulomb parameter has already been caclulated.  If so, return
    the parameter.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_colulomb_parameter"):
      self.average_coulomb_parameter = {}
    if not direction in self.average_coulomb_parameter:
      self.average_coulomb_parameter[direction] = self.calcAverageCoulombParameter(direction)
    return self.average_coulomb_parameter[direction]

  def calcAverageCoulombParameter(self,direction):
    """
    This finds the average coulomb parameter.
    """
    pi = 3.14159
    coefficient = 4*pi/3 

    average_density = self.original_phase_volume.getAverageDensity()
    average_debye_length = self.getAverageDebyeLengths(direction)

    average_coulomb_parameter = coefficient*average_density*average_debye_length**3
    return average_coulomb_parameter

  def getAverageCouplingParameter(self,direction):
    """
    Checks to see if the coupling parameter has already been caclulated.  If so, return
    the parameter.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_colulomb_parameter"):
      self.average_coupling_parameter = {}
    if not direction in self.average_coupling_parameter:
      self.average_coupling_parameter[direction] = self.calcAverageCouplingParameter(direction)
    return self.average_coupling_parameter[direction]

  def calcAverageCouplingParameter(self,direction):
    """
    This finds the average coupling parameter.
    """
    pi = 3.14159
    coefficient = self.charge/(4*pi*self.epsilon) #Charge has power one because of temp in eV!

    average_density = self.original_phase_volume.getAverageDensity()

    average_coupling_parameter = coefficient*special.cbrt(average_density)/self.getAverageTemperature(direction)
    return average_coupling_parameter

  def getTemperatures(self):
    """
    Checks to see if the temperatures have already been caclulated.  If so, return
    the previously identified temperatre.  Otherwise, calculate the temperatures and store
    if for later retrieval and returning.
    """
    if not hasattr(self,"temperatures"):
      self.temperatures = self.calcTemperatures()
    return self.temperatures

  def calcTemperatures(self):
    """
    This finds the temperatures across the grid points' ensembles
    of particles.  
    """
    temperature = {}
    temperature["x"] = []
    temperature["y"] = []
    temperature["perp"] = []
    temperature["z"] = []
    temperature["total"] = []
    denominator = self.mass*self.k_b
    for k in range(0,len(self.phase_volumes)):
      if len(self.phase_volumes[k]) <= 1:
        temperature["x"].append(-1)
        temperature["y"].append(-1)
        temperature["perp"].append(-1)
        temperature["z"].append(-1)
        temperature["total"].append(-1)
        continue
      cov_matrix = self.phase_volumes[k].getCovarianceMatrix()
      temperature["x"].append(cov_matrix.getCovarianceElement("px","px")/denominator)
      temperature["y"].append(cov_matrix.getCovarianceElement("py","py")/denominator)
      temperature["perp"].append(0.5*(temperature["x"][k]+temperature["y"][k]))
      temperature["z"].append(cov_matrix.getCovarianceElement("pz","pz")/denominator)
      temperature["total"].append(1.0/3.0*(temperature["x"][k]+temperature["y"][k]+temperature["z"][k]))
    return temperature

  def getDebyeLengths(self):
    """
    Checks to see if the Debye lengths have already been caclulated.  If so, return
    the previously identified Debye lengths.  Otherwise, calculate it.
    """
    if not hasattr(self,"debye_lengths"):
      self.debye_lengths = self.calcDebyeLengths()
    return self.debye_lengths

  def calcDebyeLengths(self):
    """
    This finds the Debye lengths across the grid spaces.
    """
    debye_lengths = {}
    debye_lengths["x"] = []
    debye_lengths["y"] = []
    debye_lengths["perp"] = []
    debye_lengths["z"] = []
    debye_lengths["total"] = []
    coefficient = np.sqrt(self.epsilon/self.charge) #Charge in square root because of temp in eV!
    temperatures = self.getTemperatures()
    for k in range(0,len(self.phase_volumes)):
      n = len(self.phase_volumes[k])
      if n <= 1 or temperatures == -1:
        debye_lengths["x"].append(-1)
        debye_lengths["y"].append(-1)
        debye_lengths["perp"].append(-1)
        debye_lengths["z"].append(-1)
        debye_lengths["total"].append(-1)
        continue
      volume = self.z_step * self.phase_volumes[k].calcRperpSq()
      density = n*self.electrons_per_macroparticle/volume
      debye_lengths["x"].append(coefficient*np.sqrt(temperatures["x"][k]/density))
      debye_lengths["y"].append(coefficient*np.sqrt(temperatures["y"][k]/density))
      debye_lengths["perp"].append(coefficient*np.sqrt(temperatures["perp"][k]/density))
      debye_lengths["z"].append(coefficient*np.sqrt(temperatures["z"][k]/density))
      debye_lengths["total"].append(coefficient*np.sqrt(temperatures["total"][k]/density))
    return debye_lengths

  def getPlasmaFrequencies(self):
    """
    Checks to see if the plasma frequencies have already been caclulated.  If so, return
    the previously identified Debye lengths.  Otherwise, calculate it.
    """
    if not hasattr(self,"plasma_frequencies"):
      self.plasma_frequencies = self.calcPlasmaFrequencies()
    return self.plasma_frequencies

  def calcPlasmaFrequencies(self):
    """
    This finds the plasma frequencies across the grid spaces.
    """
    plasma_frequencies = []
    coefficient = self.charge/np.sqrt(self.epsilon * self.mass_in_kg)
    for k in range(0,len(self.phase_volumes)):
      n = len(self.phase_volumes[k])
      if n <= 1:
        plasma_frequencies.append(-1)
        continue
      volume = self.z_step * self.phase_volumes[k].calcRperpSq()
      density = n*self.electrons_per_macroparticle/volume
      plasma_frequencies.append(coefficient*np.sqrt(density))
    return plasma_frequencies

  def getCoulombParameters(self):
    """
    Checks to see if the coulomb parameter has already been caclulated.  If so, return
    the parameter.  Otherwise, calculate it.
    """
    if not hasattr(self,"colulomb_parameters"):
      self.coulomb_parameters = self.calcCoulombParameters()
    return self.coulomb_parameters

  def calcCoulombParameters(self):
    """
    This finds the local coulomb parameters.
    """
    coulomb_parameters = {}
    coulomb_parameters["x"] = []
    coulomb_parameters["y"] = []
    coulomb_parameters["perp"] = []
    coulomb_parameters["z"] = []
    coulomb_parameters["total"] = []
    pi = 3.14159
    coefficient = 4*pi/3 
    debye_lengths = self.getDebyeLengths()

    for k in range(0,len(self.phase_volumes)):
      n = len(self.phase_volumes[k])
      if n <= 1:
        coulomb_parameters["x"].append(-1)
        coulomb_parameters["y"].append(-1)
        coulomb_parameters["perp"].append(-1)
        coulomb_parameters["z"].append(-1)
        coulomb_parameters["total"].append(-1)
        continue
      volume = self.z_step * self.phase_volumes[k].calcRperpSq()
      density = n*self.electrons_per_macroparticle/volume
      coulomb_parameters["x"].append(coefficient*density*debye_lengths["x"][k]**3)
      coulomb_parameters["y"].append(coefficient*density*debye_lengths["y"][k]**3)
      coulomb_parameters["perp"].append(coefficient*density*debye_lengths["perp"][k]**3)
      coulomb_parameters["z"].append(coefficient*density*debye_lengths["z"][k]**3)
      coulomb_parameters["total"].append(coefficient*density*debye_lengths["total"][k]**3)
    return coulomb_parameters

  def getCouplingParameters(self):
    """
    Checks to see if the coupling parameter has already been caclulated.  If so, return
    the parameter.  Otherwise, calculate it.
    """
    if not hasattr(self,"coupling_parameters"):
      self.coupling_parameters = self.calcCouplingParameters()
    return self.coupling_parameters

  def calcCouplingParameters(self):
    """
    This finds the local coulomb parameters.
    """
    coupling_parameters = {}
    coupling_parameters["x"] = []
    coupling_parameters["y"] = []
    coupling_parameters["perp"] = []
    coupling_parameters["z"] = []
    coupling_parameters["total"] = []
    pi = 3.14159
    coefficient = self.charge/(4*pi*self.epsilon) #Charge has power one because of temp in eV!
    temperatures = self.getTemperatures()

    for k in range(0,len(self.phase_volumes)):
      n = len(self.phase_volumes[k])
      if n <= 1:
        coupling_parameters["x"].append(-1)
        coupling_parameters["y"].append(-1)
        coupling_parameters["perp"].append(-1)
        coupling_parameters["z"].append(-1)
        coupling_parameters["total"].append(-1)
        continue
      volume = self.z_step * self.phase_volumes[k].calcRperpSq()
      density = n*self.electrons_per_macroparticle/volume
      coupling_parameters["z"].append(coefficient*special.cbrt(density)/temperatures["x"][k])
      coupling_parameters["y"].append(coefficient*special.cbrt(density)/temperatures["y"][k])
      coupling_parameters["perp"].append(coefficient*special.cbrt(density)/temperatures["perp"][k])
      coupling_parameters["z"].append(coefficient*special.cbrt(density)/temperatures["z"][k])
      coupling_parameters["total"].append(coefficient*special.cbrt(density)/temperatures["total"][k])
    return coupling_parameters

  def getEmittances(self):
    """
    Checks to see if the emittances have already been caclulated.  If so, return
    the previously identified temperatre.  Otherwise, calculate the temperatures and store
    if for later retrieval and returning.
    """
    if not hasattr(self,"emittances"):
      self.emittances = self.calcEmittances()
    return self.emittances

  def calcEmittances(self):
    """
    This finds the rms normalized emittances across the grid points' ensembles
    of particles.  
    """
    emittances = {}
    emittances["x"] = []
    emittances["y"] = []
    emittances["z"] = []
    for k in range(0,len(self.phase_volumes)):
      if len(self.phase_volumes[k]) <= 2:
        emittances["x"].append(-1)
        emittances["y"].append(-1)
        emittances["z"].append(-1)
        continue
      volume = self.z_step * self.phase_volumes[k].calcRperpSq()
      cov_matrix = self.phase_volumes[k].getCovarianceMatrix()
      emittances["x"].append(np.sqrt(cov_matrix.getSubDeterminant(subelements=["x","px"]))/self.mass)
      emittances["y"].append(np.sqrt(cov_matrix.getSubDeterminant(subelements=["y","py"]))/self.mass)
      emittances["z"].append(np.sqrt(cov_matrix.getSubDeterminant(subelements=["z","pz"]))/self.mass)
    return emittances

  def printTemperatures(self):
    """
    Prints to stdout the z, Tx, Ty, Tz data.
    """
    temperatures = self.getTemperatures()
    print "z,Tx,Ty,Tz,n,density"
    volume = self.z_step * self.y_step * self.x_step
    for k in range(0,len(self.phase_volumes)):
      volume = self.z_step * self.phase_volumes[k].calcRperpSq()
      z = self.z_min + (k+0.5)*self.z_step
      n = len(self.phase_volumes[k])
      density = n*self.electrons_per_macroparticle/volume
      output = []
      output.append(z)
      output.append(temperatures["x"][k])
      output.append(temperatures["y"][k])
      output.append(temperatures["z"][k])
      output.append(n)
      output.append(density)
      print ",".join([str(o) for o in output])

  def printDebyeLengths(self):
    """
    Prints to stdout the z,Lx,Ly,Lz,n,density
    """
    debye_lengths = self.getDebyeLengths()
    print "z,Lx,Ly,Lz,n,density"
    for k in range(0,len(self.phase_volumes)):
      z = self.z_min + (k+0.5)*self.z_step
      n = len(self.phase_volumes[k])
      if n<=1:
        density = -1
      else:
        volume = self.z_step * self.phase_volumes[k].calcRperpSq()
        density = n*self.electrons_per_macroparticle/volume
      output = []
      output.append(z)
      output.append(debye_lengths["x"][k])
      output.append(debye_lengths["y"][k])
      output.append(debye_lengths["z"][k])
      output.append(n)
      output.append(density)
      print ",".join([str(o) for o in output])

  def printBinnedStats(self):
    """
    Prints to stdout the z, Tx, Ty, Tz data.
    """
    pi = 3.14159
    temperatures = self.getTemperatures()
    debye_lengths = self.getDebyeLengths()
    emittances = self.getEmittances()
    plasma_frequencies = self.getPlasmaFrequencies()
    coulomb_parameters = self.getCoulombParameters()
    coupling_parameters = self.getCouplingParameters()
    print "z,Tx,Ty,Tz,Lx,Ly,Lz,ex,ey,ez,Ctotal,Cperp,Cz,Gtotal,Gperp,Gz,w,plasma_period,n,density"
    for k in range(0,len(self.phase_volumes)):
      z = self.z_min + (k+0.5)*self.z_step
      n = len(self.phase_volumes[k])
      if n<=1:
        density = -1
      else:
        volume = self.z_step * self.phase_volumes[k].calcRperpSq()
        density = n*self.electrons_per_macroparticle/volume
      output = []
      output.append(z)
      output.append(temperatures["x"][k])
      output.append(temperatures["y"][k])
      output.append(temperatures["z"][k])
      output.append(debye_lengths["x"][k])
      output.append(debye_lengths["y"][k])
      output.append(debye_lengths["z"][k])
      output.append(emittances["x"][k])
      output.append(emittances["y"][k])
      output.append(emittances["z"][k])
      output.append(coulomb_parameters["total"][k])
      output.append(coulomb_parameters["perp"][k])
      output.append(coulomb_parameters["z"][k])
      output.append(coupling_parameters["total"][k])
      output.append(coupling_parameters["perp"][k])
      output.append(coupling_parameters["z"][k])
      output.append(plasma_frequencies[k])
      output.append(2.0*pi/plasma_frequencies[k])
      output.append(n)
      output.append(density)
      print ",".join([str(o) for o in output])

  def printAverageStatistics(self,header=False):
    """
    Prints out the average statistics for the phase volume.
    """
    tx = self.getAverageTemperature('x')
    ty = self.getAverageTemperature('y')
    tz = self.getAverageTemperature('z')
    output = [tx, ty, 0.5*(tx+ty),tz,1.0/3.0*(tx+ty+tz)]
    output.append(self.original_phase_volume.getAverageDensity())
    w = self.getAveragePlasmaFrequency()
    output.append(w)
    output.append(2.0*3.14159/w)
    output.append(self.getAverageDebyeLengths("x"))
    output.append(self.getAverageDebyeLengths("y"))
    output.append(self.getAverageDebyeLengths("perp"))
    output.append(self.getAverageDebyeLengths("z"))
    output.append(self.getAverageDebyeLengths("total"))
    output.append(self.getAverageCoulombParameter("total"))
    output.append(self.getAverageCoulombParameter("perp"))
    output.append(self.getAverageCoulombParameter("z"))
    output.append(self.getAverageCouplingParameter("total"))
    output.append(self.getAverageCouplingParameter("perp"))
    output.append(self.getAverageCouplingParameter("z"))
    if header:
      print "tx,ty,tperp,tz,ttotal,density,plasma_freq,plasma_period,lx,ly,lperp,lz,ltotal,Coulombtotal,Coulombperp,Coulombz,Gammatotal,Gammaperp,Gammaz"
    print ",".join([str(o) for o in output])

class KeyedBinnedPhase6DVolume():
  """
  Replaces the lists with dicts so that empty bins can be ignored.
  """

  def __init__(self,phase_volume,number_bins_1d=None,number_bins_x=None,number_bins_y=None,
               number_bins_z=None,x_mean=None,y_mean=None,z_mean=None):
    """
    Initializes the keyed phase volume with an empty phase volume and .
    """
    self.phase_volumes = {}
    self.mass = phase_volume.mass
    if number_bins_1d is not None:
      number_bins_x = number_bins_1d
      number_bins_y = number_bins_1d
      number_bins_z = number_bins_1d
    if number_bins_x is None or number_bins_y is None or number_bins_z is None:
      raise CoordinateException("Binning the phase space requires the bins in each direction be provided.")
    self.number_bins_x = number_bins_x
    self.number_bins_y = number_bins_y
    self.number_bins_z = number_bins_z
    self.x_max = phase_volume.getMax("x")
    self.y_max = phase_volume.getMax("y")
    self.z_max = phase_volume.getMax("z")
    self.x_min = phase_volume.getMin("x")
    self.y_min = phase_volume.getMin("y")
    self.z_min = phase_volume.getMin("z")
    if x_mean is not None:
      x_diff = self.x_max - x_mean
      x_diff_lower = x_mean - self.x_min
      if x_diff < x_diff_lower:
        x_diff = x_diff_lower
      self.x_max = x_mean + x_diff
      self.x_min = x_mean - x_diff
    if y_mean is not None:
      y_diff = self.y_max - y_mean
      y_diff_lower = y_mean - self.y_min
      if y_diff < y_diff_lower:
        y_diff = y_diff_lower
      self.y_max = y_mean + y_diff
      self.y_min = y_mean - y_diff
    if z_mean is not None:
      z_diff = self.z_max - z_mean
      z_diff_lower = z_mean - self.z_min
      if z_diff < z_diff_lower:
        z_diff = z_diff_lower
      self.z_max = z_mean + z_diff
      self.z_min = z_mean - z_diff
    self.x_len = self.x_max - self.x_min
    self.y_len = self.y_max - self.y_min
    self.z_len = self.z_max - self.z_min
    self.x_step = self.x_len/number_bins_x
    self.y_step = self.y_len/number_bins_y
    self.z_step = self.z_len/number_bins_z
    self.binOriginalPhaseVolume()

  def binOriginalPhaseVolume(self):
    """
    Distributes the particles from the original phase volume through the
    binned structure.
    """
    for particle in self.original_phase_volume:
      ix = int(np.floor((particle.x.x - self.x_min)/self.x_step))
      iy = int(np.floor((particle.x.y - self.y_min)/self.y_step))
      iz = int(np.floor((particle.x.z - self.z_min)/self.z_step))
      if ix == number_bins_x:
        ix -= 1
      if iy == number_bins_y:
        iy -= 1
      if iz == number_bins_z:
        iz -= 1
      try:
        self.phase_volumes[ix][iy][iz].addParticle(particle)
      except KeyError:
        if ix not in phase_volumes.keys:
          phase_volumes[ix] = {}
        if iy not in phase_volumes[ix].keys:
          phase_volumes[ix][iy] = {}
        if iz not in phase_volumes[ix][iy].keys:
          phase_volumes[ix][iy][iz] = Phase6DVolume()
        self.phase_volumes[ix][iy][iz].addParticle(particle)
      
  def getSelfPotentialEnergy(self,particle_charge):
    """
    Checks to see if the self potential energy has already been caclulated.  If so, return
    the previously identified value.  Otherwise, calculate the value and store
    if for later retrieval and returning.
    """
    if not hasattr(self,"self_potential_energy"):
      self.self_potential_energy = self.approximateSelfPotentialEnergy(particle_charge)
    return self.self_potential_energy

  def approximateSelfPotentialEnergy(self,particle_charge):
    """
    This finds the find the exact potential energy from particles in near-by boxes (0,+/-1)
    and approximates the rest by trating the number of particles in the boxes at the center of the
    box.
    """
    energy = 0
    k = 8.99e9 #N m^2 / C^2
    numerator = k*particle_charge**2
    for ix in self.phase_volumes.keys():
      for ix2 in self.phase_volumes.keys():
        diff_ix = abs(ix-ix2)
        dist_x = diff_ix*self.x_step
        for iy in self.phase_volumes[ix].keys():
          for iy2 in self.phase_volumes[ix2].keys():
            diff_iy = abs(iy-iy2)
            dist_y = diff_iy*self.y_step
            for iz in self.phase_volumes[ix][iy].keys():
              for iz2 in self.phase_volumes[ix2][iy2].keys():
                diff_iz = abs(iz-iz2)
                dist_z = diff_iz*self.z_step
                if diff_ix == 0 and diff_iy == 0 and diff_iz == 0: #Correct energy of all of the particles in the same box.
                  energy += self.phase_volumes[ix][iy][iz].getSelfPotentialEnergy(particle_charge)
                elif diff_ix <= 1 and diff_iy <= 1 and diff_iz <= 1: #Correct energy between particle in box 2 and neighboring box particles. 
                  for particle_2 in self.phase_volumes[ix2][iy2][iz2]:
                    energy += self.phase_volumes[ix][iy][iz].getSelfPotentialEnergy(particle_charge,particle_2=particle_2)
                else: #the approximation using the potential between boxes with charge located at center (the speed up)
                  energy += numerator*len(self.phase_volumes[ix][iy][iz])*len(self.phase_volumes[ix2][iy2][iz2])/np.sqrt(dist_x**2+dist_y**2+dist_z**2)
    return energy
