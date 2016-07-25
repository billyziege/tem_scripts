import sys
import os
import shutil
import csv
import argparse
import re
import numpy
from scipy import constants
from template import fill_template
from handle_simulation_parameters import ParameterSetList
from my_data_table import MyDataTable, MyDataField

class CosyOutput(MyDataTable):
  """
  Creates a class to store the screen table in memory.
  """

  def __init__(self):
    fields = []
    fields.append(MyDataField("step number",fieldtype="int"))
    fields.append(MyDataField("number macroparticles",fieldtype="int"))
    fields.append(MyDataField("time",unit="sec"))
    fields.append(MyDataField("x",unit="m"))
    fields.append(MyDataField("y",unit="m"))
    fields.append(MyDataField("z",unit="m"))
    fields.append(MyDataField("std_x",unit="m"))
    fields.append(MyDataField("std_y",unit="m"))
    fields.append(MyDataField("std_z",unit="m"))
    fields.append(MyDataField("p_x",unit="MV/m"))
    fields.append(MyDataField("p_y",unit="MV/m"))
    fields.append(MyDataField("p_z",unit="MV/m"))
    fields.append(MyDataField("std_px",unit="MV/m"))
    fields.append(MyDataField("std_py",unit="MV/m"))
    fields.append(MyDataField("std_pz",unit="MV/m"))
    MyDataTable.__init__(self,fields)

  def injectFile(self,filepath,adjust_step_number=0):
    step_number_field = self.getFieldByName("step number")
    step_number_field.adjust = adjust_step_number
    MyDataTable.injectFile(self,filepath,delimiter="\s+")

  def writeToFile(self,**kwargs):
    MyDataTable.writeToFile(self,delimiter=" ",**kwargs)

  def addGammaField(self,mass):
    gamma_field = MyDataField("gamma")
    self.fields.addField(gamma_field)
    for row in self.rows:
      row.addCell(gamma_field,numpy.sqrt(1 + (row.getCellWithFieldname("p_z").getValue()/mass)**2))

  def addVelocityField(self,mass):
    if self.getFieldByName("gamma") is None:
      self.addGammaField(mass)
    velocity_field = MyDataField("v_z",unit="c")
    self.fields.addField(velocity_field)
    for row in self.rows:
      row.addCell(velocity_field,numpy.sqrt(1 - (1/row.getCellWithFieldname("gamma").getValue())**2))

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package')
  parser.add_argument('-s','--screen_file', dest="screen_file", type=str, help='The path to the cosy output file, sometimes called screen.txt, specifying that the read_screen_file function should be tested.  Prints out the last row of the screen file as an object.', default=None)
  parser.add_argument('-m','--number_of_electron_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons per macroparticle for the simulation.  This defaults to 100 unless specified.', default=100)
  args = parser.parse_args()

  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  mass_of_macroparticle = args.number_of_electrons_per_macroparticle*mass_of_electron

  cosy_output = CosyOutput()
  if args.screen_file is not None:
    cosy_output.injectFile(args.screen_file)
    cosy_output.addVelocityField(mass_of_macroparticle)
    print str(cosy_output.rows[-1])
