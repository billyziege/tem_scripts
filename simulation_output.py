import sys
import os
import shutil
import csv
import argparse
import re
import numpy
import math
from template import fill_template
from handle_simulation_parameters import ParameterSetList
from my_data_table import MyDataTable, MyDataField

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

  def returnEmittance(self,fieldname):
    """
    Calculates the emittance for x, y or z.
    """
    if fieldname not in ["x", "y", "z"]:
      raise Exception("Emmitance is only defined for x, y, and z")
    m0c = 0.511 #mass of electron in MeV/c, same as momentum
    p_key = "p_" + fieldname
    sum_key_squared = 0
    sum_p_key_squared = 0
    sum_key_times_p_key = 0
    n = len(self.rows)
    if n == 0:
      return 0
    for row in self.rows:
     cell1 = row.getCellWithFieldname(fieldname)
     cell2 = row.getCellWithFieldname(p_key)
     sum_key_squared += pow(cell1.getValue(),2) 
     sum_p_key_squared += pow(cell2.getValue(),2)
     sum_key_times_p_key += cell1.getValue()*cell2.getValue()
    sum_key_squared = sum_key_squared/float(n)
    sum_p_key_squared = sum_p_key_squared/float(n)
    sum_key_times_p_key = sum_key_times_p_key/float(n)
    mean_key = self.returnMean(fieldname)
    mean_p_key = self.returnMean(p_key)
    second_moment_key_squared = sum_key_squared - mean_key*mean_key
    second_moment_p_key_squared = sum_p_key_squared - mean_p_key*mean_p_key
    second_moment_key_times_p_key = sum_key_times_p_key - mean_key*mean_p_key

    difference = second_moment_key_squared*second_moment_p_key_squared - second_moment_key_times_p_key*second_moment_key_times_p_key
    return math.sqrt(difference)/m0c

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package')
  parser.add_argument('-c','--coordinates_files', dest="coordinates_files", nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space and therefore specifying that the injection function should be tested.  Prints out the last row of the first phase space file as an object.', default=None)
  parser.add_argument('-a','--average_coordinates_files', dest="average_coordinates_files", nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space and therefore specifying that the averaging function should be tested.  Prints out the standard deviation of all components.', default=None)
  parser.add_argument('-s','--std_coordinates_files', dest="std_coordinates_files", nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space and therefore specifying that the standard deviation function should be tested.  Prints out the standard deviation of all components.', default=None)
  parser.add_argument('-e','--emittance_coordinates_files', dest="emittance_coordinates_files", nargs='+', type=str, help='The path(s) to the simultion output files containing the phase space and therefore specifying that the emmitance function should be tested.  Prints out the emmitance of all components.', default=None)
  args = parser.parse_args()

  simulation_output = SimulationOutput()
  if args.coordinates_files is not None:
    simulation_output.injectFile(args.coordinates_files[0])
    print str(simulation_output.rows[-1])
  if args.average_coordinates_files is not None:
    for coordinates_file in args.average_coordinates_files:
      simulation_output.injectFile(coordinates_file)
    print len(simulation_output.rows)
    for field in simulation_output.fields:
      fieldname = field.name
      print "average " + str(field.name) + ":" + str(simulation_output.returnMean(fieldname))
  if args.std_coordinates_files is not None:
    for coordinates_file in args.std_coordinates_files:
      simulation_output.injectFile(coordinates_file)
    for field in simulation_output.fields:
      fieldname = field.name
      print "std " + str(fieldname) + ":" + str(simulation_output.returnStandardDeviation(fieldname))
  if args.emittance_coordinates_files is not None:
    for coordinates_file in args.emittance_coordinates_files:
      simulation_output.injectFile(coordinates_file)
    for fieldname in ["x","y","z"]:
      print "emittance " + str(fieldname) + ":" + str(simulation_output.returnEmittance(fieldname))

