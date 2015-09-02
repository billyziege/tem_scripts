import sys
import os
import shutil
import csv
import argparse
import re
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

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package')
  parser.add_argument('-s','--screen_file', dest="screen_file", type=str, help='The path to the cosy output file, sometimes called screen.txt, specifying that the read_screen_file function should be tested.  Prints out the last row of the screen file as an object.', default=None)
  args = parser.parse_args()

  cosy_output = CosyOutput()
  if args.screen_file is not None:
    cosy_output.injectFile(args.screen_file)
    print str(cosy_output.rows[-1])
