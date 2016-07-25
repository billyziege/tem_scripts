import sys
import os
import math
import argparse
import numpy as np
from scipy import constants
from coordinate_vector import CoordinateException
from coordinate_vector_3d import Cartesian3DVector
from my_covariance_matrix import MyCovarianceMatrix

class ParticlePhaseCoordinates():
  """
  The phase coordinates for a particle in 6D.
  """

  def __init__(self,mass,x=None,p=None,v=None):
    """
    Uses the two input Cartesian3DVectors or 0,0,0 for each to define the
    attributes x and p. 
    """
    self.setPosition(x)
    self.setMass(mass)
    if p is not None and v is not None:
      raise CoordinateException("Initializing a particle can only have momentum or velocity, not both.")
    elif p is None:
      self.setVelocity(v)
      self.calcMomentumFromVelocity()
    elif v is None:
      self.setMomentum(p)
      self.calcVelocityFromMomentum()

  def asDict(self,keys=["x","y","z","px","py","pz"]):
    """
    Returns a dict with the given keys taken from x, y, z, px, py, pz, vx, vy, vz.
    """
    out_dict = {}
    for key in keys:
      vector = self.getPosition()
      direction = key
      if key.startswith("v"):
        vector = self.getVelocity()
        direction = key.replace("v","")
      elif key.startswith("p"):
        vector = self.getMomentum()
        direction = key.replace("p","")
      out_dict[key] = getattr(vector,direction)
    return out_dict
       

  def setPosition(self,x):
    """
    The standard set function for position.
    """
    if x is None:
      self.x = Cartesian3DVector()
    else:
      if isinstance(x,Cartesian3DVector):
        self.x = Cartesian3DVector(x.x,x.y,x.z)
      else:
        raise CoordinateException("Initializing a particle with the incorrect position vector type.")

  def getPosition(self):
    """
    The standard get function for position.
    """
    return self.x

  def setMomentum(self,p):
    """
    The standard set function for momentum.
    """
    if p is None:
      self.p = Cartesian3DVector()
    else:
      if isinstance(p,Cartesian3DVector):
        self.p = Cartesian3DVector(p.x,p.y,p.z)
      else:
        raise CoordinateVector("Initializing a particle with the incorrect momentum vector type.")
  
  def getMomentum(self):
    """
    The standard get function for momentum.
    """
    return self.p

  def setMass(self,mass):
    """
    The standard set function for mass.
    """
    self.mass = mass

  def getMass(self):
    """
    The standard get function for mass.
    """
    return self.mass

  def setVelocity(self,v):
    """
    The standard set function.
    """
    if v is None:
      self.v = Cartesian3DVector()
    else:
      if isinstance(v,Cartesian3DVector):
        self.v = Cartesian3DVector(v.x,v.y,v.z)
      else:
        raise CoordinateVector("Initializing a particle with the incorrect velocity vector type.")

  def getVelocity(self):
    """
    The standard get fucntion for velocity.
    """
    return self.v

  def calcVelocityFromMomentum(self):
    """
    Calculates the cartesian 3d vector of velocity from the momentum and mass of the particle.
    """
    if self.mass is None:
      raise CoordinateVector("The particle mass needs to be specified to calculate the particle velocity from momentum.")
    values = {}
    for direction in self.p.order:
      gamma = self.calcLorentzGammaFromMomentum(direction)
      values[direction] = getattr(self.p,direction)/(gamma*self.mass) 
    self.setVelocity(Cartesian3DVector(**values))
    return self.getVelocity()

  def calcMomentumFromVelocity(self):
    """
    Calculates the cartesian 3d vector of momentum from the velocity and mass of the particle.
    """
    if self.mass is None:
      raise CoordinateVector("The particle mass needs to be specified to calculate the particle momentum from velocity.")
    values = {}
    for direction in self.v.order:
      gamma = self.calcLorentzGammaFromVelocity(direction)
      values[direction] = getattr(self.v,direction)*gamma*self.mass
    self.setMomentum(Cartesian3DVector(**values))
    return self.getMomentum()

  def calcLorentzGammaFromMomentum(self,direction):
    """
    Calculates the lorenzt gamma in the provided direction from the momentum and mass of the particle.
    """
    if self.mass is None:
      raise CoordinateVector("The particle mass needs to be specified to calculate the lorentz gamma.")
    if direction not in self.x.order:   
      raise CoordinateVector("The direction, "+str(direction)+ " needs to be one of " +",".join(self.x.order) + " to calculated the lorentz gamma.")
    speed_light = constants.physical_constants["speed of light in vacuum"][0]#m/sec by default
    return np.sqrt(1 + (getattr(self.p,direction)/(self.mass*speed_light))**2)

  def calcLorentzGammaFromVelocity(self,direction):
    """
    Calculates the lorenzt gamma in the provided direction from the velocity of the particle expressed as a fraction of c.
    """
    if direction not in self.v.order:   
      raise CoordinateVector("The direction, "+str(direction)+ " needs to be one of " +",".join(self.x.order) + " to calculated the lorentz gamma.")
    return np.sqrt(1 /(1 - (getattr(self.v,direction))**2))

  def advancePosition(self,time):
    """
    Assuming no forces, advances the particle's position over the provided time.
    """
    speed_light = constants.physical_constants["speed of light in vacuum"][0]#m/sec by default
    velocity = speed_light*self.getVelocity()#conversion from c units to m/sec
    return self.x + time*velocity

  def translate(self,translation_vector):
    """
    Returns a new particle with the translated position.
    """
    if isinstance(translation_vector,Cartesian3DVector):
      new_particle = ParticlePhaseCoordinates(self.mass,self.x-translation_vector,self.p)
      return new_particle
    raise CoordinateException("Translating a particle with the incorrect translation vector type.")

  def getEnergy(self):
    """
    Wraps the calc energy method skipping it if the energy is already stored.
    """
    if not hasattr(self,"energy"):
      self.energy = self.calcEnergy()
    return self.energy

  def calcEnergy(self):
    """
    Calculates the energy of the particle in MeV.
    """
    speed_light = constants.physical_constants["speed of light in vacuum"][0]#m/sec by default
    if self.mass is None:
      raise CoordinateVector("The particle mass needs to be specified to calculate the energy.")
    return np.sqrt(self.p*self.p + self.mass**2)

  def getValueFromFieldname(self,fieldname):
    """
    Provides an interface between the "x", "y", "z", "px", "py", "pz", and "E"
    syntax and the way the data is stored in this object.
    """
    if fieldname == "E":
      return self.getEnergy()
    momentum_direction = fieldname.replace("p","")
    if fieldname.startswith("p") and momentum_direction in ["x","y","z"]:
      return getattr(self.p,momentum_direction)
    elif fieldname in ["x","y","z"]:
      return getattr(self.x,fieldname)
    raise Exception("The given field, "+fieldname+", is not defined for the particle.")

  def __str__(self):
    """
    Passes the string argument to the poistion and then momentum.
    """
    output = []
    output.append(str(self.x))
    output.append(str(self.p))
    return " ".join(output)

class Phase6DVolume():
  """
  An ensemble of 3D position and 3D momentum vectors.
  """

  def __init__(self, mass=None):
    """
    Initializes the phase volume to an empty list.
    """
    self.mass = mass
    self.particles = []

  def __len__(self):
    return len(self.particles)

  def addParticle(self,mass,**kwargs):
    """
    If x is a particle, then it is added, otherwise, a new particle 
    is initialized with the cartesian coordinates x and p.
    """
    if self.mass is None:
      self.mass = mass
    if self.mass != mass:
      raise Exception("Right now phase volumes are only supported forparticles of the same mass.")
    if isinstance(mass,ParticlePhaseCoordinates):
      self.particles.append(mass)
    else:
      particle = ParticlePhaseCoordinates(mass,**kwargs)
      self.particles.append(particle)

  def translate(self,translation_vector):
    """
    Translates the phase volume to a new coordinate systems
    and returns the translated phase volume.
    """
    new_phase_volume = Phase6DVolume()
    for particle in self:
      new_phase_volume.addParticle(particle.translate(translation_vector))
    return new_phase_volume

  def __iter__(self):
    """
    Allows the ability to iterate over the particles in
    phase space without referencing the particels themselves.
    """
    for particle in self.particles:
      yield particle

  def __str__(self):
    """
    Wraps the str function so that the entire phase space may be printed.
    """
    output = []
    for particle in self:
      output.append(str(particle))
    return "\n".join(output)

  def injectFile(self,filepath,mass=1,header=False,fieldnames=["x","y","z","px","py","pz"],delimiter=" "):
    """
    Reads in the file from filepath and extracts the phasespace according to
    the split on the delimiter.  Fieldnames are only used if header is false, otherwise
    they are derived from the header.  If no fieldname called "mass" is used, the mass of
    the particle is set to mass.
    """
    with open(filepath,'r') as f:
      for line in f:
        line = line.rstrip()
        line = line.lstrip(" ")
        if header == True:
          header = False
          fieldnames = line.split(delimiter)
          continue
        pieces = line.split(delimiter)
        pieces = [p for p in pieces if p != '']
        if len(pieces) != len(fieldnames):
          print pieces
          print fieldnames
          raise Exception("The format of " + filepath + " is inconsistent.")
        row = dict(zip(fieldnames,pieces))
        position = Cartesian3DVector(row["x"],row["y"],row["z"])
        momentum = Cartesian3DVector(row["px"],row["py"],row["pz"])
        self.addParticle(mass,x=position,p=momentum)
        

  def getMean(self,fieldnames):
    """
    Checks to see if the mean has already been caclulated.  If so, return
    the previously identified mean.  Otherwise, calculate the mean and store
    if for later retrieval and returning.
    """
    if not hasattr(self,"means"):
      self.means = {}
    sortedfieldnames = sorted(fieldnames)
    fieldname_key = ",".join(sortedfieldnames)
    if fieldname_key not in self.means:
      self.means[fieldname_key] = self.calcMean(fieldnames)
    return self.means[fieldname_key]

  def calcMean(self,fieldnames):
    """
    This finds the mean of the product of the fieldnames across the ensemble
    of particles.  The mean is calculated with 1/N where N is the number of particles.
    I thought about this, but I think correcting to 1/(N-len(fieldnames)) is not
    appropriate in this instance.
    """
    total = 0
    for particle in self:
      product = 1.
      for fieldname in fieldnames:
        product *= particle.getValueFromFieldname(fieldname)
      total += product
    return total/len(self)

  def boostCoordinatesToCOM(self,direction):
    """ Boosts the coordinates and momentum to the reference frame of the average momentum/position of the provided fieldname.
      Notice that these coordinates are then not at the same time.
      Args:
          fieldname : the position fieldname, i.e. either x, y, or z
      Returns:
          simulation_output : with the boosted coordinates
    """
    speed_light = constants.physical_constants["speed of light in vacuum"][0]
    position = [getattr(particle.x,direction) for particle in self]
    momentum = [getattr(particle.p,direction) for particle in self]
    mean_position = numpy.mean(position)
    mean_momentum = numpy.mean(momentum)
    mean_energy_divided_by_massxc = numpy.sqrt((mean_momentum/self.mass)**2+speed_light*2) 
    lorentz_gamma =  np.sqrt(1 + (mean_momentum/(self.mass*speed_light))**2)
    new_phase_volume = Phase6DVolume()
    for particle in self:
        #The move to spatial COM
        new_position = copy.deepcopy(particle.x)
        direction_position = getattr(new_position,direction)
        new_direction_position = (direction_positon-mean_position)*lorentz_gamma
        setattr(new_position,direction,new_direction_position)

        #The boost to the momentum COM 
        new_momentum = copy.deepcopy(particle.p)
        direction_momentum = getattr(particle.p,direction)
        energy_divided_by_massxc = numpy.sqrt((direction_momentum/self.mass)**2+speed_light*2) 
        new_direction_momentum = lorentz_gamma * (direction_momentum - mean_momentum*energy_divided_by_massxc/mean_energy_divided_by_massxc)
        setattr(new_momentum,direction,new_direction_momentum)
        
        new_phase_volume.addParticle(particle.mass,x=new_position,p=new_momentum)
    return new_phase_volume
  
  def getCovarianceMatrix(self,recalculate=False):
    """
    Wraps the calculation of the covariance matrix.  If it has previously been
    calculated, then skip the calculation.  Returns the 7 x 7 matrix.
    """
    if recalculate or not hasattr(self,"cov_matrix"):
      self.cov_matrix = MyCovarianceMatrix(self)
    return self.cov_matrix

class MeanPhase6DVolume(Phase6DVolume):
  """
  An ensemble of 3D position and 3D momentum vectors taken
  from the means of phase 6D volumes.
  """

  def addParticle(self,phase_volume):
    """
    Adds the 6D mean of the phase volume as a particle
    to the ensemble of mean particles.
    """
    if not isinstance(phase_volume,Phase6DVolume):
      raise Exception("Points in the mean space must be generated from Phase6DVolumes only.")
    x_mean = Cartesian3DVector(phase_volume.getMean("x"),phase_volume.getMean("y"),phase_volume.getMean("z"))
    p_mean = Cartesian3DVector(phase_volume.getMean("px"),phase_volume.getMean("py"),phase_volume.getMean("pz"))
    particle = ParticlePhaseCoordinates(phase_volume.mass,x=x_mean,p=p_mean)
    self.particles.append(particle)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package and use simple commands to get some of the straightforward methods.')
  parser.add_argument('coordinates_files', nargs='+', type=str, help='The path(s) to the phase volume files containing the phase space and therefore specifying that the injection function should be tested.')
  parser.add_argument('-c','--covariance_matrix', dest="covariance_matrix", action="store_true", help='Prints out the upper triangular form of the covariance matrix.', default=False)
  parser.add_argument('-r','--correlation_matrix', dest="correlation_matrix", action="store_true", help='Prints out the upper triangular form of the correlation matrix.', default=False)
  parser.add_argument('-e','--emittances', dest="emittances", action="store_true", help='Prints out the emittances for the phase volume.', default=False)
  parser.add_argument('--mixed_matrix', dest="mixed_matrix", action="store_true", help='Prints out the upper triangular form of the correlation matrix with variance along diagonal.', default=False)
  parser.add_argument('--sub_determinant', dest="sub_determinant", action="store_true", help='Prints the determinant of the phase space portion of the covariance matrix.', default=False)
  parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
  args = parser.parse_args()

  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

  phase_volume = Phase6DVolume()
  for path in args.coordinates_files:
    phase_volume.injectFile(path,mass=mass_of_macroparticle)

  if args.covariance_matrix:
    cov_matrix = phase_volume.getCovarianceMatrix()
    cov_matrix.printCovarianceMatrix()

  if args.correlation_matrix:
    cov_matrix = phase_volume.getCovarianceMatrix()
    cov_matrix.printCorrelationMatrix()

  if args.mixed_matrix:
    cov_matrix = phase_volume.getCovarianceMatrix()
    cov_matrix.printMixedMatrix()

  if args.sub_determinant:
    cov_matrix = phase_volume.getCovarianceMatrix()
    print cov_matrix.getSubDeterminant()
 
  if args.emittances:
    cov_matrix = phase_volume.getCovarianceMatrix()
    print "ex"
    print np.sqrt(cov_matrix.getSubDeterminant(["x","px"]))/mass_of_macroparticle
    print "ey"
    print np.sqrt(cov_matrix.getSubDeterminant(["y","py"]))/mass_of_macroparticle
    print "ez"
    print np.sqrt(cov_matrix.getSubDeterminant(["z","pz"]))/mass_of_macroparticle


