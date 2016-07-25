import sys
import os
import math
import argparse
import numpy as np
from scipy import constants
from my_covariance_matrix import MyEditableCovarianceMatrix
from coordinate_vector_3d import Cartesian3DVector
from phase_volume_6d import Phase6DVolume

class ConventionalCOMBoost():
  """
  The conventional boost, <p^lab>, built into it.
  """

  def __init__(self,phase_volume,direction="z",**kwargs):
    if not isinstance(phase_volume,Phase6DVolume):
      raise Exception("A Phase6DVolume needs to be used to set up the conventional COM boost.")
    self.mass = phase_volume.particles[0].mass #Just the mass of the first particle, but all particles should have same mass
    if direction not in ["x","y","z"]:
      raise Exception("The only acceptable boosts are in the x, y, or z directions.")
    self.direction = direction
    self.lorentz_gamma = None
    self.boost_momentum = None
    self.boost_velocity = None
    self.boosted_cov_matrix = None
    self.translation_vector = self.getTranslationVector(phase_volume)
    self.phase_volume = phase_volume.translate(self.translation_vector)
    self.getCOMLorentzGamma(**kwargs)

  def getTranslationVector(self,phase_volume):
    """
    Calculates the appropriate translation vector to
    take the boosted coordinates average to 0,0,0. 
    """
    return Cartesian3DVector(x=phase_volume.getMean(["x"]),y=phase_volume.getMean(["y"]),z=phase_volume.getMean(["z"]))

  def setLorentzGamma(self,value):
    """
    This overwrites the value for the lorentz gamma
    with the provided value.
    """
    self.lorentz_gamma = value

  def setBoostMomentum(self,value):
    """
    This overwrites the value for the boost momentum
    with the provided value.
    """
    self.boost_momentum = value

  def setBoostVelocity(self,value):
    """
    This overwrites the value for the boost velocity
    with the provided value.
    """
    self.boost_velocity = value

  def getCOMLorentzGamma(self,recalculate=False,**kwargs):
    """
    Wraps the calculation of the lorentz gamma
    and returns it. Unitless always.
    """
    if self.lorentz_gamma is not None and not recalculate:
      return self.lorentz_gamma
    if self.boost_momentum is not None and not recalculate:
      self.setLorentzGamma(np.sqrt(1 + (self.boost_momentum/self.mass)**2))
    elif self.boost_velocity is not None and not recalculate:
      self.setLorentzGamma(1./np.sqrt(1. - self.boost_velocity**2))
    else:
      p_boost = self.getCOMBoostMomentum(recalculate=recalculate,**kwargs)
      self.setLorentzGamma(np.sqrt(1 + (p_boost/self.mass)**2))
    return self.lorentz_gamma

  def getCOMBoostMomentum(self,recalculate=False,**kwargs):
    """
    Wraps the calculation of the average momentum
    and returns it.  I assume that the units are Energy/c.
    """
    if self.boost_velocity is not None and not recalculate:
      lorentz_gamma = self.getCOMLorentzGamma(recalculate=recalculate,**kwargs)
      self.setBoostMomentum(lorentz_gamma*self.mass*self.boost_velocity)
    if self.boost_momentum is None or recalculate:
      self.setBoostMomentum(self.phase_volume.getMean(["p"+self.direction]))
    return self.boost_momentum

  def getCOMBoostVelocity(self,recalculate=False,**kwargs):
    """
    Returns v/c.  Note this is unitless.
    """
    if recalculate or self.boost_velocity is None:
      if self.lorentz_boost is None or recalculate:
        self.getLorentzBoost(recalculate=recalculate,**kwargs)
      if self.boost_momentum is None or recalculate:
        self.getCOMBoostMomentum(recalculate=recalculate,**kwargs)
      self.setBoostVelocity(self.boost_momentum/(self.lorentz_gamma*self.mass))
    return self.boost_velocity

  def getBoostedCovarianceMatrix(self,recalculate=False):
    """
    Wraps the calculation so that the boosted covariance matrix can be stored.
    """
    if recalculate or self.boosted_cov_matrix is None:
      self.boosted_cov_matrix = self.boostCovarianceMatrix()
    return self.boosted_cov_matrix

  def boostCovarianceMatrix(self):
    """
    Boosts the covariance matrix adjusting the compoents of the covariance matrix
    according to the rules of the boost.
    """   
    boosted_cov_matrix = MyEditableCovarianceMatrix()
    order = ["x", "y", "z", "px", "py", "pz", "E"]
    for e1 in order:
      for e2 in order:
        boosted_cov_matrix.setCovarianceElement(e1,e2,self.boostCovElement(e1,e2))
    return boosted_cov_matrix

  def boostCovElement(self,e1,e2):
    """
    Provides the boosting rules for each of the elements of the
    covariance matrix.  Again, these are the conventional rules.
    """
    multiplicative = 1
    lorentz_gamma = self.getCOMLorentzGamma()
    if e1 == self.direction or e1 == "p"+self.direction:
      multiplicative *= lorentz_gamma
    if e2 == self.direction or e2 == "p"+self.direction:
      multiplicative *= lorentz_gamma
    if e1 == "E" or e2 == "E": #E is assumed constant in the conventional boost.
      multiplicative = 0
    return multiplicative*self.phase_volume.getCovarianceMatrix().getCovarianceElement(e1,e2)
    
class EVarCOMBoost(ConventionalCOMBoost):
  """
  The COM boost without the assumption that var(E) = 0
  """

  def getCOMLorentzGamma(self,recalculate=False,**kwargs): #Redefined because velocity is used to calculate the other values.
    """
    Wraps the calculation of the lorentz gamma
    and returns it. Unitless always.
    """
    if self.lorentz_gamma is not None and not recalculate:
      return self.lorentz_gamma
    if self.boost_momentum is not None and not recalculate:
      self.setLorentzGamma(np.sqrt(1 + (self.boost_momentum/self.mass)**2))
    elif self.boost_velocity is not None and not recalculate:
      self.setLorentzGamma(1./np.sqrt(1. - self.boost_velocity**2))
    else:
      v_boost = self.getCOMBoostVelocity(recalculate=recalculate,**kwargs)
      self.setLorentzGamma(1/np.sqrt(1- v_boost**2))
    return self.lorentz_gamma

  def getCOMBoostMomentum(self,recalculate=False,**kwargs):
    """
    Returns v/c.  Note this is unitless.
    """
    if recalculate or self.boost_momentum is None:
      if self.lorentz_gamma is None or recalculate:
        self.getCOMLorentzGamma(recalculate=recalculate,**kwargs)
      if self.boost_velocity is None or recalculate:
        self.getCOMBoostVelocity(recalculate=recalculate,**kwargs)
      self.setBoostMomentum(self.boost_velocity*self.lorentz_gamma*self.mass)
    return self.boost_velocity

  def getCOMBoostVelocity(self,recalculate=False,**kwargs):
    """
    Returns v/c = c <p_direction>/<E>.
      direction must be one of x, y, or z.
    """
    if recalculate or self.boost_velocity is None:
      if self.boost_momentum is not None and not recalculate:
        lorentz_gamma = self.getCOMLorentzGamma(recalculate=recalculate,**kwargs)
        self.setBoostVelocity(self.boost_momentum/(lorentz_gamma*self.mass))
      else:
        mean_momentum = self.phase_volume.getMean(["p"+self.direction])
        mean_energy = self.phase_volume.getMean(["E"])
        self.boost_velocity = mean_momentum/mean_energy
    return self.boost_velocity

  def boostCovElement(self,e1,e2):
    """
    Provides the boosting rules for each of the elements of the
    covariance matrix.  Again, these are the conventional rules.
    """
    lorentz_gamma = self.getCOMLorentzGamma()
    lab_matrix = self.phase_volume.getCovarianceMatrix()
    v = self.getCOMBoostVelocity()
    evar_adjustment = 0
    if e1 == "E" and e2 == "E": 
      ee = lab_matrix.getCovarianceElement("E","E")
      pze = lab_matrix.getCovarianceElement("p"+self.direction,"E")
      pzpz = lab_matrix.getCovarianceElement("p"+self.direction,"p"+self.direction)
      evar_adjustment = lorentz_gamma**2 * (ee -  2*v*pze + v**2 *pzpz)
    elif e1 == "E" or e2 == "E": 
      if e1 == "p"+self.direction or e2 == "p"+self.direction:
        ee = lab_matrix.getCovarianceElement("E","E")
        pze = lab_matrix.getCovarianceElement("p"+self.direction,"E")
        pzpz = lab_matrix.getCovarianceElement("p"+self.direction,"p"+self.direction)
        evar_adjustment = lorentz_gamma**2 * ( (1 + v**2)*pze - v*(pzpz+ee) )
      elif e1 == self.direction or e2 == self.direction:
        ze = lab_matrix.getCovarianceElement(self.direction,"E")
        zpz = lab_matrix.getCovarianceElement(self.direction,"p"+self.direction)
        evar_adjustment = lorentz_gamma**2 * ( ze - v*zpz )
      else:
        ze = lab_matrix.getCovarianceElement(self.direction,"E")
        zpz = lab_matrix.getCovarianceElement(self.direction,"p"+self.direction)
        evar_adjustment = lorentz_gamma * ( ze - v*zpz )
    elif e1 == "p"+self.direction and e2 == "p"+self.direction: 
      ee = lab_matrix.getCovarianceElement("E","E")
      pze = lab_matrix.getCovarianceElement("p"+self.direction,"E")
      evar_adjustment = lorentz_gamma**2 * (v**2 * ee -  2*v*pze)
    elif e1 == "p"+self.direction or e2 == "p"+self.direction:
      if e1 == "p"+self.direction:
        other = e2
        p = e1
      else:
        other = e1
        p = e2
      multiplicative = lorentz_gamma
      if other == self.direction:
        multiplicative *= lorentz_gamma
      if other in ["x","y","z","px","py","pz"]: #Remeber that p + self.direction is already excluded, so this is only 5 not 6.
        ze = lab_matrix.getCovarianceElement(other,"E")
        evar_adjustment = multiplicative * (-v) * ze
    return ConventionalCOMBoost.boostCovElement(self,e1,e2)+evar_adjustment

class PTCOMBoost(EVarCOMBoost):
  """
  The boost both without the assumption Var(E) = 0 and the addition
  of (x,y,z)-(px,py,py)/m * t to the transformation to tcom = 0.
  """

  def getTranslationVector(self,phase_volume):
    """
    Calculates the appropriate translation vector to
    take the boosted coordinates after moving time to 0 to average to 0,0,0. 
    """
    #The boost needs to be calculated before the translation... but this is still exact since the translation does nothing to the pu four vector.
    mean_momentum = phase_volume.getMean(["p"+self.direction])
    mean_energy = phase_volume.getMean(["E"])
    boost_velocity = mean_momentum/mean_energy
    lorentz_gamma = 1./np.sqrt(1-boost_velocity**2)
    p_boost = lorentz_gamma*self.mass*boost_velocity
    translation = {}
    translation[self.direction] = phase_volume.getMean([self.direction]) + lorentz_gamma*p_boost/self.mass * (phase_volume.getMean([self.direction,"p"+self.direction]) - p_boost*phase_volume.getMean([self.direction,"E"])) 
    for direction in ["x","y","z"]:
      if direction == self.direction:
        continue
      translation[direction] = phase_volume.getMean([self.direction]) + lorentz_gamma*p_boost/self.mass * (translation[self.direction] * phase_volume.getMean(["p"+direction]) - phase_volume.getMean([self.direction,"p"+direction]))
    return Cartesian3DVector(x=translation["x"],y=translation["y"],z=translation["z"])

  def boostCovElement(self,e1,e2):
    """
    Provides the boosting rules for each of the elements of the
    covariance matrix.  These are the terms I derived.  I apologize for the 
    naming of the variables... but at least they are used right after they 
    are named.  BTW, zzpzpz = <z*z*pz*pz> and zE = <z*E> if self.direction = z.  Etc.
    """
    lorentz_gamma = self.getCOMLorentzGamma()
    lab_matrix = self.phase_volume.getCovarianceMatrix()
    v = self.getCOMBoostVelocity()
    speed_light = constants.physical_constants["speed of light in vacuum"][0]#m/sec by default
    coeff = lorentz_gamma*v/(self.mass*speed_light**2) #A common coefficient that appears everywhere.  This simplifies code.  Has units inverse energy.
    pt_adjustment = 0
    if e1 == self.direction and e2 == self.direction: 
      zzpzpz = self.phase_volume.getMean([self.direction,self.direction,"p"+self.direction,"p"+self.direction])
      zzpz = self.phase_volume.getMean([self.direction,self.direction,"p"+self.direction])
      zpz = self.phase_volume.getMean([self.direction,"p"+self.direction])
      z = self.phase_volume.getMean([self.direction])
      zzEE = self.phase_volume.getMean([self.direction,self.direction,"E","E"])
      zzE = self.phase_volume.getMean([self.direction,self.direction,"E"])
      zE = self.phase_volume.getMean([self.direction,"E"])
      zzpzE = self.phase_volume.getMean([self.direction,self.direction,"p"+self.direction,"E"])
      pt_adjustment = 2*lorentz_gamma**2 * coeff * ( (zzpz-zpz*z-v*(zzE-zE*z)) + coeff*(zzpzpz-zpz**2-2*v*(zzpzE-zpz*zE)+v**2*(zzEE-zE**2)) )
    elif e1 == self.direction or e2 == self.direction: 
      if e1 == self.direction:
        other = e2
        current = e1
      else:
        other = e1
        current = e2
      if other == "p"+self.direction:
        zpzpz = self.phase_volume.getMean([self.direction,"p"+self.direction,"p"+self.direction])
        zpz = self.phase_volume.getMean([self.direction,"p"+self.direction])
        pz = self.phase_volume.getMean(["p"+self.direction])
        zpzE = self.phase_volume.getMean([self.direction,"p"+self.direction,"E"])
        zE = self.phase_volume.getMean([self.direction,"E"])
        zEE = self.phase_volume.getMean([self.direction,"E","E"])
        z = self.phase_volume.getMean([self.direction])
        E = self.phase_volume.getMean(["E"])
        pt_adjustment = lorentz_gamma**2 * coeff * (zpzpz-zpz*pz -v*(2*zpzE-zpz*E-zE*pz)+v**2*(zEE-zE*E))
      elif other == "E":
        zpzE = self.phase_volume.getMean([self.direction,"p"+self.direction,"E"])
        zpzpz = self.phase_volume.getMean([self.direction,"p"+self.direction,"p"+self.direction])
        zpz = self.phase_volume.getMean([self.direction,"p"+self.direction])
        pz = self.phase_volume.getMean(["p"+self.direction])
        zEE = self.phase_volume.getMean([self.direction,"E","E"])
        zE = self.phase_volume.getMean([self.direction,"E"])
        E = self.phase_volume.getMean(["E"])
        pt_adjustment = lorentz_gamma**2 * coeff * (zpzE-zpz*E-v*(zpzpz-zpz*pz+zEE-zE*E)+v**2*(zpzE-pz*zE)) 
      elif other in ["x","y","z"]: #Remeber that self.direction is already excluded, so this is only 2 not 3.
        zzpapz = self.phase_volume.getMean([self.direction,self.direction,"p"+other,"p"+self.direction])
        zzpaE = self.phase_volume.getMean([self.direction,self.direction,"p"+other,"E"])
        zzpa = self.phase_volume.getMean([self.direction,self.direction,"p"+other])
        zpa = self.phase_volume.getMean([self.direction,"p"+other])
        z = self.phase_volume.getMean([self.direction])
        azpz =  self.phase_volume.getMean([other,self.direction,"p"+self.direction])
        zpz = self.phase_volume.getMean([self.direction,"p"+self.direction])
        a = self.phase_volume.getMean([other])
        azE =  self.phase_volume.getMean([other,self.direction,"E"])
        zE = self.phase_volume.getMean([self.direction,"E"])
        pt_adjustment = lorentz_gamma*coeff*( zzpa-zpa*z+azpz-zpz*a-v*(azE-zE*a)+coeff*(zzpapz-zpa*zpz-v*(zzpaE-zpa*zE)) )
      elif other in ["px","py","pz"]:#Remeber that p+self.direction is already excluded, so this is only 2 not 3.
        zpapz = self.phase_volume.getMean([self.direction,other,"p"+self.direction])
        zpz = self.phase_volume.getMean([self.direction,"p"+self.direction])
        pa = self.phase_volume.getMean([other])
        zpaE = self.phase_volume.getMean([self.direction,other,"E"])
        zE = self.phase_volume.getMean([self.direction,"E"])
        pt_adjustment = lorentz_gamma*coeff*( zpapz-zpz*pa - v*(zpaE-zE*pa) )
    elif e1 in ["x","y","z"] or e2 in ["x","y","z"]: #Remeber that self.direction is already excluded, so this is only 2 not 3.
      if e1 in ["x","y","z"]:
        other = e2
        current = e1
      else:
        other = e1
        current = e2
      if other=="p"+self.direction:
        zpapz = self.phase_volume.getMean([self.direction,"p"+current,"p"+self.direction])
        zpz = self.phase_volume.getMean([self.direction,"p"+self.direction])
        pa = self.phase_volume.getMean(["p"+current])
        zpaE = self.phase_volume.getMean([self.direction,"p"+current,"E"])
        zE = self.phase_volume.getMean([self.direction,"E"])
        pt_adjustment = lorentz_gamma*coeff*( zpapz-zpz*pa - v*(zpaE-zE*pa) )
      if other=="E":
        zpapz = self.phase_volume.getMean([self.direction,"p"+current,"p"+self.direction])
        zpaE = self.phase_volume.getMean([self.direction,"p"+current,"E"])
        zpa = self.phase_volume.getMean([self.direction,"p"+current])
        E = self.phase_volume.getMean(["E"])
        pz = self.phase_volume.getMean(["p"+self.direction])
        pt_adjustment = coeff*( zpaE-zpa*E - v*(zpapz-zpa*pz) )
      if other in ["x","y","z"]:#Remeber that self.direction is already excluded, so this is only 2 not 3.
        zzpapb = self.phase_volume.getMean([self.direction,self.direction,"p"+current,"p"+other])
        azpb = self.phase_volume.getMean([current,self.direction,"p"+other])
        bzpa = self.phase_volume.getMean([other,self.direction,"p"+current])
        zpa = self.phase_volume.getMean([self.direction,"p"+current])
        zpb = self.phase_volume.getMean([self.direction,"p"+other])
        a = self.phase_volume.getMean([current])
        b = self.phase_volume.getMean([other])
        pt_adjustment = coeff*(azpb-a*zpb+bzpa-b*zpa + coeff*(zzpapb-zpa*zpb))
      if other in ["px","py","pz"]:#Remeber that p+self.direction is already excluded, so this is only 2 not 3.
        zpapb = self.phase_volume.getMean([self.direction,"p"+current,other])
        zpa = self.phase_volume.getMean([self.direction,"p"+current])
        pb = self.phase_volume.getMean([other])
        pt_adjustment = coeff*(zpapb-zpa*pb)
    return EVarCOMBoost.boostCovElement(self,e1,e2)+pt_adjustment


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package and use simple commands to get some of the straightforward methods.')
  parser.add_argument('coordinates_files', nargs='+', type=str, help='The path(s) to the phase volume files containing the phase space and therefore specifying that the injection function should be tested.')
  parser.add_argument('-c','--conventional_boosted_covariance_matrix', dest="conventional", action="store_true", help='Prints out the upper triangular form of the conventionally boosted covariance matrix.', default=False)
  parser.add_argument('-e','--evar_boosted_covariance_matrix', dest="evar", action="store_true", help='Prints out the upper triangular form of the Energy varianted boosted covariance matrix.', default=False)
  parser.add_argument('-p','--pt_boosted_covariance_matrix', dest="pt", action="store_true", help='Prints out the upper triangular form of the full energy varianted boosted covariance matrix while taking account of the velocity of the particles in the COM frame.', default=False)
  parser.add_argument('--sub_determinant', dest="sub_determinant", action="store_true", help='Prints the determinant of the phase space portion of the covariance matrix.', default=False)
  parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
  args = parser.parse_args()

  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

  phase_volume = Phase6DVolume()
  for path in args.coordinates_files:
    phase_volume.injectFile(path,mass=mass_of_macroparticle)

  if args.conventional:
    con_boost = ConventionalCOMBoost(phase_volume)
    print con_boost.getCOMLorentzGamma()
    boosted_cov_matrix = con_boost.getBoostedCovarianceMatrix()
    boosted_cov_matrix.printCovarianceMatrix()
    if args.sub_determinant:
      print boosted_cov_matrix.getSubDeterminant()

  if args.evar:
    evar_boost = EVarCOMBoost(phase_volume)
    print evar_boost.getCOMLorentzGamma()
    boosted_cov_matrix = evar_boost.getBoostedCovarianceMatrix()
    boosted_cov_matrix.printCovarianceMatrix()
    print  boosted_cov_matrix.getCovarianceElement("z","pz")
    if args.sub_determinant:
      print boosted_cov_matrix.getSubDeterminant()

  if args.pt:
    pt_boost = PTCOMBoost(phase_volume)
    print pt_boost.getCOMLorentzGamma()
    boosted_cov_matrix = pt_boost.getBoostedCovarianceMatrix()
    boosted_cov_matrix.printCovarianceMatrix()
    if args.sub_determinant:
      print boosted_cov_matrix.getSubDeterminant()

