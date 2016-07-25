import sys
import os
import shutil
import csv
import argparse
import re
import numpy
import math
from scipy import constants
from template import fill_template
from handle_simulation_parameters import ParameterSetList
from my_data_table import MyDataTable, MyDataField, MyDataRow

class SimulationOutput(MyDataTable):
  """
  Creates a class to store the phase space table in memory.
  """

  def __init__(self):
    fields = []
    fields.append(MyDataField("x",unit="m"))
    fields.append(MyDataField("y",unit="m"))
    fields.append(MyDataField("z",unit="m"))
    fields.append(MyDataField("p_x",unit="MV/m"))
    fields.append(MyDataField("p_y",unit="MV/m"))
    fields.append(MyDataField("p_z",unit="MV/m"))
    MyDataTable.__init__(self,fields)

  def injectFile(self,coordinates_file):
    """
    Reads in a coordinates file as a list of dictionary with 
    the following keys (stuff in parenthesis are comments and not part of the key):
     x
     y
     z
     p_x 
     p_y 
     p_z 
    """
    MyDataTable.injectFile(self,coordinates_file,delimiter="\s+")

  def convertCoordinatesToBetterUnits(self,mass,position_conversion=10**6,time_conversion=10**12):
    """
    Converts the units of the output.  By default, from m to um, sec to psec, and kg to units of mass of particles. NEED TO CONVERT UNITS STILL.
      Args:
          mass (real) : mass of particles in kg
      Returns:
          simulation_output : with the converted coordinates
    """
    mass_conversion = 1/mass
    speed_light = constants.physical_constants["speed of light in vacuum"][0]*position_conversion/time_conversion#in um/psec by default
    eV_to_J_conversion = constants.physical_constants["electron volt-joule relationship"][0]*mass_conversion*position_conversion**2/time_conversion**2 #From eV to m_part*um^2/psec^2
    momentum_conversion = eV_to_J_conversion*10**6/speed_light #From MeV/c to m_part um/psec by default
    new_simulation_output = SimulationOutput()
    for row in self.rows:
      fieldnames = []
      for field in self.fields.list:
        fieldnames.append(field.name)
      values = []
      for name in fieldnames:
        value = row.getCellWithFieldname(name).getValue()
        if name.startswith("p"):
          value = value*momentum_conversion
        else:
          value = value*position_conversion
        values.append(value)
      new_row = MyDataRow(self.fields,values)
      new_simulation_output.rows.append(new_row)
    return new_simulation_output
    

  def boostCoordinates(self,fieldname,position_conversion=10**6,time_conversion=10**12):
    """ Boosts the coordinates and momentum to the reference frame of the average momentum/position of the provided fieldname.
      Args:
          fieldname : the position fieldname, i.e. either x, y, or z
      Returns:
          simulation_output : with the boosted coordinates
    """
    speed_light = constants.physical_constants["speed of light in vacuum"][0]*position_conversion/time_conversion#in um/psec by default
    position = self.returnAllColumnValues(fieldname)
    momentum = self.returnAllColumnValues("p_"+fieldname)
    mean_position = numpy.mean(position)
    mean_momentum = numpy.mean(momentum)
    #mean_momentum = self.calcWeightedMeanMomentum(fieldname)
    mean_energy_divided_by_massxc = numpy.sqrt(mean_momentum**2+speed_light*2) 
    lorentz_gamma = self.calcLorentzGamma(fieldname)
    new_simulation_output = SimulationOutput()
    for row in self.rows:
      fieldnames = []
      for field in self.fields.list:
        fieldnames.append(field.name)
      values = []
      for name in fieldnames:
        value = row.getCellWithFieldname(name).getValue()
        if name == fieldname:
          value = (value-mean_position)*lorentz_gamma
        if name == "p_" + fieldname:
          energy_divided_by_massxc = numpy.sqrt(value**2+speed_light*2) 
          value = lorentz_gamma * (value - mean_momentum*energy_divided_by_massxc/mean_energy_divided_by_massxc)
        values.append(value)
      new_row = MyDataRow(self.fields,values)
      new_simulation_output.rows.append(new_row)
    return new_simulation_output
    

  def calcWeightedMeanMomentum(self,fieldname,position_conversion=10**6,time_conversion=10**12):
    """Funciton to calculate the average of momentum divided by energy time m_part *speel_light.  Assumes momentum has mass units of m_part
      Args:
          fieldname : the position fieldname, i.e. either x, y, or z
      Returns:
          velocity (real): value of average velocity
    """
    speed_light = constants.physical_constants["speed of light in vacuum"][0]*position_conversion/time_conversion#in um/psec by default
    momentum = self.returnAllColumnValues("p_"+fieldname)
    return numpy.mean([p/numpy.sqrt(1+(p/speed_light)**2) for p in momentum])
    

  def calcEmittance(self,fieldname):
      """ Function to calculate normalized emittance
      Args:
          fieldname : the position fieldname, i.e. either x, y, or z
      Returns:
          emittance (real): value of emittance
      """
      position = self.returnAllColumnValues(fieldname)
      momentum = self.returnAllColumnValues("p_"+fieldname)
      mean_positionxmomentum = numpy.mean([x*p for x,p in zip(position,momentum)])
      mean_position_squared = numpy.mean([x*x for x in position])
      mean_momentum_squared = numpy.mean([p*p for p in momentum])
      emittance = numpy.sqrt(mean_position_squared * mean_momentum_squared - (mean_positionxmomentum)**2)
      return emittance
  
  def calcLorentzGamma(self,fieldname,position_conversion=10**6,time_conversion=10**12):
      """ Function to calculate lorenzt boost in fieldname direction.  Assuming momentum has mass units expressed in m_part
      and we are boosting to the mean velocity frame.
      Args:
          fieldname : the position fieldname, i.e. either x, y, or z
      Returns:
          gamma (real): value of Lorentz gamma
      """
      speed_light = constants.physical_constants["speed of light in vacuum"][0]*position_conversion/time_conversion#in um/psec by default
      #mean_momentum = self.calcWeightedMeanMomentum(fieldname)
      momentum = self.returnAllColumnValues("p_"+fieldname)
      mean_momentum = numpy.mean(momentum)
      return numpy.sqrt(1 + (mean_momentum/speed_light)**2)

  def calcEnergySpread(self, fieldname):
      """ Function to calculate energy spread
      Args:
          fieldname : the position fieldname, i.e. either x, y, or z
      Returns:
          gamma (real): value of Lorentz gamma
          E_spread (real): value of energy spread
      """
      #position = self.returnAllColumnValues(fieldname)
      #momentum = self.returnAllColumnValues("p_"+fieldname)
      #mean_position = numpy.mean(position)
      #mean_momentum = numpy.mean(momentum)
      #mean_positionxmomentum = numpy.mean([x*p for x,p in zip(position,momentum)])
      #var_position = numpy.var(position, ddof = 1)
      #var_momentum = numpy.var(momentum, ddof = 1)
      #gamma = numpy.sqrt(1 + (mean_momentum/mass)**2)
      lorentz_gamma = self.calcLorentzGamma(fieldname)
      #delta_p = numpy.sqrt(var_momentum - (mean_positionxmomentum - mean_position * mean_momentum)**2/var_position)
      delta_p = self.calcEta(fieldname)
      delta_gamma = numpy.sqrt(lorentz_gamma**2 -1)/lorentz_gamma * delta_p
      E_spread = delta_gamma/lorentz_gamma
      return lorentz_gamma, E_spread
  
  def calcBunchLength(self,fieldname,position_conversion=10**6,time_conversion=10**12):
      """ Function to calculate electron bunch length 
      Args:
          fieldname : the position fieldname, i.e. either x, y, or z
      Returns:
          numpy.std(t_dist) (real): standard deviation of temporal distribution
      """
      speed_light = constants.physical_constants["speed of light in vacuum"][0]*position_conversion/time_conversion#in um/psec by default
      position = self.returnAllColumnValues(fieldname)
      momentum = self.returnAllColumnValues("p_"+fieldname)
      vel = [p/numpy.sqrt(1 + p**2) * speed_light for p in momentum]
      t_dist = [x/v for x,v in zip(position,vel)]
      return numpy.std(t_dist)
  
  def calcVariance(self,fieldname,ddof=0):
     """ Function to calculate the of a fieldname
     Args:
       fieldname : the fieldname for which the variance should be calculated, i.e. x, y, z, p_x, p_y, or p_z
     Returns:
       sigma of fieldname : the variance of the provided fieldname
     """
     field_value = self.returnAllColumnValues(fieldname)
     return numpy.var(field_value, ddof = ddof)

  def calcGamma(self,fieldname):
    """ Function to calculate the gamma of the distribution. Units are m_part um^2/psec.
    Args:
      fieldname : the position fieldname, i.e. either x, y, or z
    Returns:
      gamma of the fieldname : <xp> - <x><p>
    """
    position = self.returnAllColumnValues(fieldname)
    momentum = self.returnAllColumnValues("p_"+fieldname)
    mean_position = numpy.mean(position)
    mean_momentum = numpy.mean(momentum)
    mean_positionxmomentum = numpy.mean([x*p for x,p in zip(position,momentum)])
    return (mean_positionxmomentum - mean_position*mean_momentum)

  def calcEta(self,fieldname):
    """ Function to calculate the eta of the distribution.  Units are m_part^2 um^2/psec^2.
    Args:
      fieldname : the position fieldname, i.e. either x, y, or z
    Returns:
      gamma of the fieldname : <p^2> - <p>^2 - (<xp> - <x><p>)^2/(<x^2> - <x>^2)
    """
    position = self.returnAllColumnValues(fieldname)
    momentum = self.returnAllColumnValues("p_"+fieldname)
    var_position = numpy.var(position,ddof=0)
    var_momentum = numpy.var(momentum,ddof=0)
    gamma = self.calcGamma(fieldname)
    return (var_momentum - gamma**2/var_position)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package and make reports for the output of the simulation statistics.')
  parser.add_argument('coordinates_files', nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space and therefore specifying that the injection function should be tested.')
  parser.add_argument('-a','--average', dest="average", action="store_true", help='Prints out the mean of all components.', default=False)
  parser.add_argument('-s','--std', dest="std", action="store_true", help='Prints out the standard deviation of all components.', default=False)
  parser.add_argument('-v','--var', dest="var", action="store_true", help='Prints out the variance of all components.', default=False)
  parser.add_argument('-e','--emittance', dest="emittance", action="store_true", help='Prints out the emmitance of the spatial components.', default=False)
  parser.add_argument('-n','--energy_spread', dest="energy_spread", action="store_true", help='Prints out the energy spread in the longitudinal direction and the kinectic energy.', default=False)
  parser.add_argument('-b','--bunch_length', dest="bunch_length", action="store_true", help='Prints out the bunch_length in the longitudinal direction.', default=False)
  parser.add_argument('-c','--coherence_flux', dest="coherence_flux", action="store_true", help='Prints out the coherence flux using x in the transversal direction.', default=False)
  parser.add_argument('-g','--gamma', dest="gamma", action="store_true", help='Prints out the gamma of the spatial component.',default=False)
  parser.add_argument('-z','--eta', dest="eta", type=int, help='Prints out the eta of the spatial component.',default=False)
  parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
  parser.add_argument('-r','--relativistic_boost', dest="boost", action="store_true", help='Tells the program to boost the coordinates.  Prints the boosted and unboosted coorindates side by side.', default=False)
  parser.add_argument('-f','--full_report', dest="full_report", action='store_true', help='Creates a full report of the provided coordinate files with all of the possible statistics.', default=None)
  args = parser.parse_args()

  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron
  simulation_output = SimulationOutput()
  simulation_output.injectFile(args.coordinates_files[0])
  simulation_output.convertCoordinatesToBetterUnits(mass_of_macroparticle)
  if args.boost:
    boosted_simulation_output = simulation_output.boostCoordinates("z")
    for i in range(len(simulation_output.rows)):
      output = []
      output.append(simulation_output.returnTableValue("z",i))
      output.append(boosted_simulation_output.returnTableValue("z",i))
      output.append(simulation_output.returnTableValue("p_z",i))
      output.append(boosted_simulation_output.returnTableValue("p_z",i))
      print ",".join([str(o) for o in output])
    output = []
    output.append(numpy.mean(simulation_output.returnAllColumnValues("z")))
    output.append(numpy.mean(boosted_simulation_output.returnAllColumnValues("z")))
    output.append(numpy.mean(simulation_output.returnAllColumnValues("p_z")))
    output.append(numpy.mean(boosted_simulation_output.returnAllColumnValues("p_z")))
    output.append(numpy.var(simulation_output.returnAllColumnValues("z")))
    output.append(numpy.var(boosted_simulation_output.returnAllColumnValues("z")))
    output.append(numpy.var(simulation_output.returnAllColumnValues("p_z")))
    output.append(numpy.var(boosted_simulation_output.returnAllColumnValues("p_z")))
    print ",".join([str(o) for o in output])
    exit()
  simulation_output = simulation_output.boostCoordinates("z")
  if args.average or args.full_report:
    for field in simulation_output.fields.list:
      fieldname = field.name
      print( "average %s: %.3e %s" % (str(field.name),simulation_output.returnMean(fieldname),str(simulation_output.getFieldAttrByName(fieldname,"unit"))) )
  if args.std or args.full_report:
    for field in simulation_output.fields.list:
      fieldname = field.name
      print( "std %s: %.3e %s" % (str(field.name),simulation_output.returnStandardDeviation(fieldname),str(simulation_output.getFieldAttrByName(fieldname,"unit"))) )
  if args.var or args.full_report:
    for field in simulation_output.fields.list:
      fieldname = field.name
      print( "var %s: %.3e %s" % (str(field.name),simulation_output.calcVariance(fieldname),str(simulation_output.getFieldAttrByName(fieldname,"unit")) + "^2") )
  if args.emittance or args.full_report:
    for fieldname in ["x","y","z"]:
      print( "emittance in %s: %.3e m" % (str(fieldname),simulation_output.calcEmittance(fieldname)))
  if args.energy_spread or args.full_report:
    lorentz_gamma, energy_spread = simulation_output.calcEnergySpread("z")
    K = mass_of_electron * (lorentz_gamma - 1)
    print( "longitudinal energy spread: %.3e MeV" % (energy_spread) )
    print( "Kinetic energy (calculated only in longitudinal direction): %.3e MeV" % (K) )
  if args.gamma or args.full_report:
    for fieldname in ["x","y","z"]:
      print( "gamma in %s: %.3e " % (str(fieldname),simulation_output.calcGamma(fieldname)))
  if args.eta or args.full_report:
    for fieldname in ["x","y","z"]:
      print( "eta in %s: %.3e " % (str(fieldname),simulation_output.calcEta(fieldname)))
  if args.bunch_length or args.full_report:
    bunch_length = simulation_output.calcBunchLength("z")
    FWHM_bunch_length = bunch_length * 2 * numpy.sqrt(2.0 * numpy.log(2.0))
    print( "Bunch length = %.3e sec" % (bunch_length) )
    print( "FWHM Bunch length = %.3e sec" % (FWHM_bunch_length) )
  if args.full_report:
    h = constants.physical_constants["Planck constant in eV s"][0]
    speed_light = constants.physical_constants["speed of light in vacuum"][0]
    constant = h*speed_light/(2.0*mass_of_electron*1e6)
    Ne = len(simulation_output.rows) * args.number_of_electrons_per_macroparticle
    print("Ne %e" %Ne)
    coh_flux = constant**2 * Ne/(simulation_output.calcEmittance("x"))**2
    print("Coh. flux = %e " % coh_flux)

