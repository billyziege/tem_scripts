import sys
import os
import math
import argparse
import numpy as np
from scipy import constants
from coordinate_vector import CoordinateException
from coordinate_vector_3d import Cartesian3DVector
from my_covariance_matrix import MyCovarianceMatrix
from my_covariance_matrix import MyCylindricalCovarianceMatrix

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
    self.setR()
    self.setTheta()

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

  def setR(self):
    """
    The standard set function for r determined from x and y.
    """
    self.r = np.sqrt(self.x.x**2 + self.x.y**2)

  def getR(self):
    """
    The standard get function for r determined from x and y.
    """
    return self.r

  def setTheta(self):
    """
    The standard set function for theta from x and y.
    """
    if self.r == 0:
      self.theta = 0
      return
    self.theta = np.arccos(self.x.x/self.r)
    if self.x.y < 0:
      self.theta += np.pi

  def getTheta(self):
    """
    The standard get function for theta determined from x and y.
    """
    return self.theta

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
        raise CoordinateException("Initializing a particle with the incorrect momentum vector type.")
  
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
        raise CoordinateException("Initializing a particle with the incorrect velocity vector type.")

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
      raise CoordinateException("The particle mass needs to be specified to calculate the particle velocity from momentum.")
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
      raise CoordinateException("The particle mass needs to be specified to calculate the particle momentum from velocity.")
    speed_light = constants.physical_constants["speed of light in vacuum"][0]#m/sec by default
    values = {}
    for direction in self.v.order:
      gamma = self.calcLorentzGammaFromVelocity(direction)
      values[direction] = getattr(self.v,direction)*gamma*self.mass*speed_light
    self.setMomentum(Cartesian3DVector(**values))
    return self.getMomentum()

  def calcLorentzGammaFromMomentum(self,direction):
    """
    Calculates the lorenzt gamma in the provided direction from the momentum and mass of the particle.
    """
    if self.mass is None:
      raise CoordinateException("The particle mass needs to be specified to calculate the lorentz gamma.")
    if direction not in self.x.order:   
      raise CoordinateException("The direction, "+str(direction)+ " needs to be one of " +",".join(self.x.order) + " to calculated the lorentz gamma.")
    speed_light = constants.physical_constants["speed of light in vacuum"][0]#m/sec by default
    return np.sqrt(1 + (getattr(self.p,direction)/(self.mass*speed_light))**2)

  def calcLorentzGammaFromVelocity(self,direction):
    """
    Calculates the lorenzt gamma in the provided direction from the velocity of the particle expressed as a fraction of c.
    """
    if direction not in self.v.order:   
      raise CoordinateException("The direction, "+str(direction)+ " needs to be one of " +",".join(self.x.order) + " to calculated the lorentz gamma.")
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

  def getEnergy(self,**kwargs):
    """
    Wraps the calc energy method skipping it if the energy is already stored.
    """
    if not hasattr(self,"energy") or self.energy is None:
      self.energy = self.calcEnergy(**kwargs)
    return self.energy

  def calcEnergy(self,units="MeV",**kwargs):
    """
    Calculates the kinetic energy of the particle in MeV.
    """
    speed_light = constants.physical_constants["speed of light in vacuum"][0]#m/sec by default
    if self.mass is None:
      raise CoordinateException("The particle mass needs to be specified to calculate the energy.")
    if units=="MeV":
      return np.sqrt(self.p*self.p + self.mass**2)
    elif units=="SI":
      return np.sqrt(self.p*self.p + self.mass**2*speed_light**2)*speed_light
    else:
      raise CoordinateException("The units need to be specified.")

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
    if isinstance(mass,ParticlePhaseCoordinates):
      self.particles.append(mass)
      if self.mass is None:
        self.mass = mass.mass
      return
    if self.mass is None:
      self.mass = mass
    if self.mass != mass:
      raise Exception("Right now phase volumes are only supported for particles of the same mass: "+ str(self.mass) + "," + str(mass))
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

  def injectFile(self,filepath,mass=1,momentum_weight=1,velocity_weight=1.,header=False,fieldnames=["x","y","z","px","py","pz"],delimiter=" "):
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
        try:
          momentum = Cartesian3DVector(row["px"],row["py"],row["pz"])
          momentum = (1.0/momentum_weight)*momentum #Rescales momentum
          self.addParticle(mass,x=position,p=momentum)
        except KeyError:
          velocity = Cartesian3DVector(row["vx"],row["vy"],row["vz"])
          velocity = velocity_weight*velocity
          self.addParticle(mass,x=position,v=velocity)
        

  def filterByRadius(self,min_r=None,max_r=None):
    """
    Returns a new phase volume from the provided phase volume but only with
    particles that lie between min_r and max_r.  If min_r or max_r are None,
    the min_r and max_r are determined from the min and max, respectively, of
    all particles in the phase space.  That is, if both are None, then 
    the phase volume returned is identical...
    """
    if min_r is None or max_r is None:
      current_min = None
      current_max = None
      for particle in self:
        if min_r is None:
          try: 
            if particle.getR() < current_min:
              current_min = particle.getR()
          except TypeError:
            current_min = particle.getR()
        if max_r is None:
          try: 
            if particle.getR() > current_max:
              current_max = particle.getR()
          except TypeError:
            current_max = particle.getR()
      if min_r is None:
        min_r = current_min
      if max_r is None:
        max_r = current_max
    output_phase_volume = Phase6DVolume()
    for particle in self:
      if particle.getR() >= min_r and particle.getR() <= max_r:
        output_phase_volume.addParticle(particle)
    return output_phase_volume

  def filterByTheta(self,min_theta=None,max_theta=None):
    """
    Returns a new phase volume from the provided phase volume but only with
    particles that lie between min_theta and max_theta.  If min_theta or max_theta are None,
    the min_theta and max_theta are determined from the min and max, respectively, of
    all particles in the phase space.  That is, if both are None, then 
    the phase volume returned is identical...
    """
    if min_theta is None or max_theta is None:
      current_min = None
      current_max = None
      for particle in self:
        if min_theta is None:
          try: 
            if particle.getTheta() < current_min:
              current_min = particle.getTheta()
          except TypeError:
            current_min = particle.getTheta()
        if max_theta is None:
          try: 
            if particle.getTheta() > current_max:
              current_max = particle.getTheta()
          except TypeError:
            current_max = particle.getTheta()
      if min_theta is None:
        min_theta = current_min
      if max_theta is None:
        max_theta = current_max
    output_phase_volume = Phase6DVolume()
    for particle in self:
      if particle.getTheta() >= min_r and particle.getTheta() <= max_r:
        output_phase_volume.addParticle(particle)
    return output_phase_volume

  def getSelfPotentialEnergy(self,particle_charge,**kwargs):
    """
    Checks to see if the electric potential energy of the
    ensemble has already been calculated.  If so, returns
    this value.  Otherwise, calculates and stores this 
    value for later retrieval.
    """
    if not hasattr(self,"potential_energy"):
      self.potential_energy = self.calcSelfPotentialEnergy(particle_charge,**kwargs)
    return self.potential_energy

  def calcSelfPotentialEnergy(self,particle_charge,particle_2=None):
    """
    Uses the k q^2/r electric potential energy to calculate
    the potential energy stored in the entire ensemble of particles.
    If particle_2 is not none, the potential energy between particle_2
    and all particles in the phase volume are returned.
    """
    k = 8.99e9 #N m^2 / C^2
    numerator = k*particle_charge**2
    number_particles = len(self)
    pe = 0
    if particle_2 is None:
      for i in range(0,number_particles-1):
        particle_i = self.particles[i]
        for j in range(i+1,number_particles):
          particle_j = self.particles[j]
          r = abs(particle_j.x - particle_i.x)
          pe += 1/r
    else:
      for i in range(0,number_particles-1):
        particle_i = self.particles[i]
        r = abs(particle_2.x - particle_i.x)
        pe += 1/r
    return numerator*pe

  def getExtractionFieldPotentialEnergy(self,particle_charge,field_strength):
    """
    Checks to see if the extraction field potential energy of the
    ensemble has already been calculated.  If so, returns
    this value.  Otherwise, calculates and stores this 
    value for later retrieval.
    """
    if not hasattr(self,"extraction_energy"):
      self.extraction_energy = self.calcExtractionFieldPotentialEnergy(particle_charge,field_strength)
    return self.extraction_energy

  def calcExtractionFieldPotentialEnergy(self,particle_charge,field_strength):
    """
    Calculates the potential energy stored in the extraction field 
    (assumed to be in the z direction only) for the entire ensemble of particles.
    """
    pe = 0
    numerator = field_strength*particle_charge
    for particle in self:
      pe += particle.x.z
    return pe*numerator

  def getKineticEnergy(self,**kwargs):
    """
    Checks to see if the kinetic energy of the
    ensemble has already been calculated.  If so, returns
    this value.  Otherwise, calculates and stores this 
    value for later retrieval.
    """
    if not hasattr(self,"kinetic_energy"):
      self.kinetic_energy = self.calcKineticEnergy(**kwargs)
    return self.kinetic_energy

  def calcKineticEnergy(self,unit_conversion="MeV",**kwargs):
    """
    Sums up the energies (kinetic) of all of the particles 
    and returns it with SI units.
    """
    ke = 0.0
    for particle in self:
      if unit_conversion == "SI":
        options = {"units":"SI"}
      elif unit_conversion == "MeV":
        options = {"units":"MeV"}
      else:
        raise CoordinateExceptionError("This function only supports SI and MeV unit choices.")
      ke += particle.getEnergy(**options)
    if unit_conversion == "MeV":
      ke = ke/(constants.physical_constants["electron volt-joule relationship"][0]*1e6) #Converts to Joules.
    return ke

  def getTotalEnergy(self,particle_charge,field_strength,**kwargs):
    """
    Gets the kinetic, self electric potential, and extraction
    field potential energies and returns their sum.
    """
    ke = self.calcKineticEnergy(**kwargs)
    fpe = self.calcExtractionFieldPotentialEnergy(particle_charge,field_strength)
    spe = self.calcSelfPotentialEnergy(particle_charge)
    return ke + fpe + spe

  def getEnergySpread(self,**kwargs):
    """
    Checks to see if the energy spread has already been calculated, and if
    so returns it.  Otherwise, calculates and stores this value for
    later retrieval.
    """
    if not hasattr(self,"energy_spread"):
      self.energy_spread = self.calcEnergySpread(**kwargs)
    return self.energy_spread

  def calcEnergySpread(self,**kwargs):
    """
    Returns the standard deviation of the energy distribution.
    """
    energies = np.empty(len(self))
    index = 0
    for particle in self:
      energies[index] = particle.getEnergy(**kwargs)
      index += 1
    return np.std(energies)

  def getCoherenceLength(self,**kwargs):
    """
    Checks to see if the coherence length has already been calculated, and if
    so returns it.  Otherwise, calculates and stores this value for
    later retrieval.
    """
    if not hasattr(self,"coherence_length"):
      self.coherence_length = self.calcCoherenceLength(**kwargs)
    return self.coherence_length

  def calcCoherenceLength(self,**kwargs):
    """
    Calculates the coherence length v*h/energy_spread
    in m.
    """
    speed_light = constants.physical_constants["speed of light in vacuum"][0] #m/sec by default
    h = constants.physical_constants["Planck constant in eV s"][0]/1e6 #Planck's constant in MeV/s

    mean_pz = self.getMean(["pz"])
    gamma = np.sqrt(1 + (mean_pz/(self.mass*speed_light))**2)
    mean_vz = mean_pz / (gamma*self.mass)*speed_light
    energy_spread = self.getEnergySpread()

    return mean_vz*h/energy_spread

  def getMax(self,fieldnames):
    """
    Checks to see if the max has already been caclulated.  If so, return
    the previously identified max.  Otherwise, calculate the max and store
    if for later retrieval and r`eturning.
    """
    if not hasattr(self,"maxs"):
      self.maxs = {}
    sortedfieldnames = sorted(fieldnames)
    fieldname_key = ",".join(sortedfieldnames)
    if fieldname_key not in self.maxs:
      self.maxs[fieldname_key] = self.calcMax(fieldnames)
    return self.maxs[fieldname_key]

  def calcMax(self,fieldnames):
    """
    This finds the max of the product of the fieldnames across the ensemble
    of particles.  
    """
    current_max = None
    for particle in self:
      product = 1.
      for fieldname in fieldnames:
        product *= particle.getValueFromFieldname(fieldname)
      if current_max is None:
        current_max = product
        continue
      elif current_max < product:
        current_max = product
    return current_max

  def getMin(self,fieldnames):
    """
    Checks to see if the min has already been caclulated.  If so, return
    the previously identified min.  Otherwise, calculate the min and store
    if for later retrieval and returning.
    """
    if not hasattr(self,"mins"):
      self.mins = {}
    sortedfieldnames = sorted(fieldnames)
    fieldname_key = ",".join(sortedfieldnames)
    if fieldname_key not in self.mins:
      self.mins[fieldname_key] = self.calcMin(fieldnames)
    return self.mins[fieldname_key]

  def calcMin(self,fieldnames):
    """
    This finds the min of the product of the fieldnames across the ensemble
    of particles.  
    """
    current_min = None
    for particle in self:
      product = 1.
      for fieldname in fieldnames:
        product *= particle.getValueFromFieldname(fieldname)
      if current_min is None:
        current_min = product
        continue
      elif current_min > product:
        current_min = product
    return current_min

  def getMaxR(self):
    """
    Checks to see if the rmax has already been calculated.  If so, return
    the previously identified rmax.  Otherwise, calculate the rmax and store
    if for later retrieval and returning.
    """
    if "rmax" not in self.__dict__:
      self.rmax = self.calcMaxR()
    return self.rmax

  def calcMaxR(self):
    """
    This finds the max radius in xy across the ensemble
    of particles.  
    """
    current_max = None
    for particle in self:
      r = np.sqrt(particle.x.x*particle.x.x + particle.x.y*particle.x.y)
      if current_max is None:
        current_max = r
        continue
      elif current_max < r:
        current_max = r
    return current_max

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

  def getAverageDensity(self):
    """
    Checks to see if the average density has already been caclulated.  If so, return
    the previously identified vale.  Otherwise, calculate it.
    """
    if not hasattr(self,"average_density"):
      self.average_density = self.calcAverageDensity()
    return self.average_density

  def calcAverageDensity(self,**kwargs):
    """
    Calculates the average density.
    """
    pi = 3.14159
    r_perp_sq = self.calcRperpSq()
    r_z_sq = self.calcRzSq()
    
    return 3*len(self)/(4*pi*r_perp_sq*np.sqrt(r_z_sq))

  def getDose(self):
    """
    Checks to see if the dose has already been caclulated.  If so, return
    the previously identified vale.  Otherwise, calculate it.
    """
    if not hasattr(self,"dose"):
      self.dose = self.calcDose()
    return self.dose

  def calcDose(self,**kwargs):
    """
    Calculates the average density.
    """
    pi = 3.14159
    r_perp_sq = self.calcRperpSq()
    
    return len(self)/(pi*r_perp_sq)

  def calcRperpSq(self,scale=2.0):
    """
    Calculates <x^2 + y^2> * scale for estimating
    volume elements of rotationally symmetric distributions.
    """
    x_mean = self.getMean(["x"])
    y_mean = self.getMean(["y"])
    r_perp_sq = 0
    for particle in self:
      r_perp_sq += 0.5*( (particle.x.x-x_mean)**2 + (particle.x.y-y_mean)**2 )
    return scale*r_perp_sq/len(self)
    
  def calcRzSq(self,scale=2.0):
    """
    Calculates <z^2> * scale for estimating
    volume elements of rotationally symmetric distributions.
    """
    z_mean = self.getMean(["z"])
    r_z_sq= 0
    for particle in self:
      r_z_sq +=(particle.x.z-z_mean)**2
    return scale*r_z_sq/len(self)
    
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

  def getCylindricalCovarianceMatrix(self,recalculate=False):
    """
    Wraps the calculation of the covariance matrix.  If it has previously been
    calculated, then skip the calculation.  Returns the 7 x 7 matrix.
    """
    if recalculate or not hasattr(self,"cov_matrix"):
      self.cyl_cov_matrix = MyCylindricalCovarianceMatrix(self)
    return self.cyl_cov_matrix

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
  parser.add_argument('--cylindrical', dest="cylindrical", action="store_true", help='Tells the program we are working in cyclindrical coordinates.', default=False)
  parser.add_argument('--cylindrical_emittances', dest="emittances", action="store_true", help='Prints out the cylindrical emittances for the phase volume.', default=False)
  parser.add_argument('--mixed_matrix', dest="mixed_matrix", action="store_true", help='Prints out the upper triangular form of the correlation matrix with variance along diagonal.', default=False)
  parser.add_argument('--sub_determinant', dest="sub_determinant", action="store_true", help='Prints the determinant of the phase space portion of the covariance matrix.', default=False)
  parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
  parser.add_argument('--filter_r_min', dest="r_min", type=float, help='Tells the script to filter the phase volume returning particles with r greater than or equal to r_min.', default=None)
  parser.add_argument('--filter_r_max', dest="r_max", type=float, help='Tells the script to filter the phase volume returning particles with r smaller than or eaual to r_max.', default=None)
  parser.add_argument('--filter_theta_min', dest="theta_min", type=float, help='Tells the script to filter the phase volume returning particles with theta greater than or equal to theta_min.', default=None)
  parser.add_argument('--filter_theta_max', dest="theta_max", type=float, help='Tells the script to filter the phase volume returning particles with theta greater than or equal to theta_max.', default=None)
  args = parser.parse_args()

  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

  phase_volume = Phase6DVolume()
  for path in args.coordinates_files:
    phase_volume.injectFile(path,mass=mass_of_electron,momentum_weight=args.number_of_electrons_per_macroparticle)

  if args.r_min is not None or args.r_max is not None:
    phase_volume = phase_volume.filterByRadius(min_r=args.r_min,max_r=args.r_max)

  if args.cylindrical:
    cov_matrix = phase_volume.getCylindricalCovarianceMatrix()
  else:
    cov_matrix = phase_volume.getCovarianceMatrix()

  if args.covariance_matrix:
    cov_matrix.printCovarianceMatrix()

  if args.correlation_matrix:
    cov_matrix.printCorrelationMatrix()

  if args.mixed_matrix:
    cov_matrix.printMixedMatrix()

  if args.sub_determinant:
    cov_matrix = phase_volume.getCovarianceMatrix()
    print cov_matrix.getSubDeterminant() 

  if args.emittances:
    output = []
    if args.cylindrical:
      output.append(np.sqrt(cov_matrix.getSubDeterminant(["rho","prho"]))/mass_of_electron)
      output.append(np.sqrt(cov_matrix.getSubDeterminant(["phi","pphi"]))/mass_of_electron)
      output.append(np.sqrt(cov_matrix.getSubDeterminant(["z","pz"]))/mass_of_electron)
    else:
      ex = np.sqrt(cov_matrix.getSubDeterminant(["x","px"]))/mass_of_electron
      output.append(np.sqrt(cov_matrix.getSubDeterminant(["x","px"]))/mass_of_electron)
      ey = np.sqrt(cov_matrix.getSubDeterminant(["y","py"]))/mass_of_electron
      output.append(np.sqrt(cov_matrix.getSubDeterminant(["y","py"]))/mass_of_electron)
      ez = np.sqrt(cov_matrix.getSubDeterminant(["z","pz"]))/mass_of_electron
      output.append(np.sqrt(cov_matrix.getSubDeterminant(["z","pz"]))/mass_of_electron)
      n = len(phase_volume)*args.number_of_electrons_per_macroparticle
      output.append(n)
      output.append(n/(ex*ey))
      output.append(n/(ex*ey*ez))
    print ",".join([str(o) for o in output])


