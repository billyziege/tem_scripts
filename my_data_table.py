import sys
import os
import shutil
import csv
import argparse
import re
import numpy
import math

class MyDataException(Exception):
  pass

class MyDataField():
  """
  A class to store the fields for a data table.
  """

  def __init__(self,fieldname,fieldtype="float",unit=None):
    """
    Initializes the field with the default type of float.  
    """
    self.setName(fieldname)
    self.setType(fieldtype) 
    self.setUnit(unit)

  def setName(self,fieldname):
    self.name = fieldname

  def getName(self):
    return self.name

  def hasName(self,fieldname):
    return self.name == fieldname

  def setType(self,fieldtype):
    self.type = fieldtype

  def getType(self):
    return self.type

  def hasType(self,fieldtype):
    return self.type == fieldtype

  def setUnit(self,unit):
    self.unit = unit

  def getUnit(self):
    return self.unit

class MyDataCell():
  """
  A class to store the data in a specific column/row of a table.
  Right now it only stores the value and field, but this class will allow a
  hook to manipulate cells independentally.  Field connects the
  cell back to the columns of the table.
  """

  def __init__(self,field,value):
    self.setField(field)
    self.setValue(value,field)

  def setField(self,field):
    self.field = field

  def getField(self):
    return self.field

  def hasFieldName(self,fieldname):
    return self.field.hasName(fieldname)

  def setValue(self,value,field):
    if field.type == "float":
      if hasattr(field,"adjust"):
        self.value = float(value) + field.adjust
      else:
        self.value = float(value)
    elif field.type == "int":
      if hasattr(field,"adjust"):
        self.value = int(value) + field.adjust
      else:
        self.value = int(value)
    else:
      self.value = str(value)

  def getValue(self):
     return self.value

  def __str__(self):
    return str(self.getValue())

class MyDataRow():
  """
  A class to store the data in a spefic row of a table.  Right
  now it is pretty basic, but it provides a hook for further 
  development.
  """

  def __init__(self,fields=[],values=[],comment=None):
    """
    Fields is a list of MyDataFields
    and values is a corresponding list of
    values for those fields.
    """
    self.cells = []
    self.addData(fields,values)
    self.setComment(comment)

  def __str__(self):
    output = []
    for cell in self.cells:
      output.append(str(cell))
    return ",".join(output)

  def addData(self,fields,values):
    """
    First checks to make sure that the fields and
    values are of the same length, then creates 
    and fills the cells of the row.
    """
    if len(fields) != len(values):
      message = []
      message.append("The number of elements identified: " + str(len(values)))
      message.append("The expected number of elements identfied:" + str(len(fields)))
      message.append("The elements: " + str(values))
      raise MyDataException("\n".join(message))
    for i in range(len(fields)):
      self.cells.append(MyDataCell(fields[i],values[i]))
        

  def setComment(self,comment):
    self.comment = comment

  def getComment(self):
    return self.comment

  def getCellWithFieldname(self,fieldname):
    """
    Returns the cell of the provided row with that has the provided fieldname.
    """
    for cell in self.cells:
      if cell.hasFieldName(fieldname):
        return cell
    return None

class MyDataTable():
  """
  A class to store  a generalization of data tables.
  """

  def __init__(self,fields=[]):
    """
    Fields need to be defined before the table can be created.
    """
    self.rows = []
    self.fields = fields

  def addField(self,fieldname,fieldtype="float"):
    """
    Checks to make sure that the field does not violate any expectations
    and then adds the field if it passes.
    """
    if self.getFieldByName is not None:
      raise MyDataException("There are multiple fields with the fieldname " + fieldname + ".")
    field = MyDataField(fieldname,fieldtype)
    self.fields.append(field)
    return field

  def getFieldByName(self,fieldname):
    """
    Returns the first field with the provided fieldname or None.  There should
    at most be one field with this fieldname unless there is an error.
    """
    for field in self.fields:
      if field.hasName(fieldname):
        return field
    return None

  def injectFile(self,filepath,header=False,ignore_first_line=False,delimiter=','):
    """
    Reads in the data from filepath, splits it by the delimiter 
    (can be regular expression),
    and puts it into the table in the SAME ORDER as the tables 
    fields. Ignores the first line if ignore_first_line is True
    and skips any lines that start with "#" and any pieces of lines
    after a "#", i.e. # indicates a comment.
    """
    if header is True:
      orig_field_order = self.fields    
    with open(filepath) as f:
      if header is True:
        ignore_first_line = False
        line = f.readline()
        while line.startswith("#"):
          line = f.readline
        line = line.rstrip()
        line = line.lstrip()
        line, comment = re.split("#", line, 1) #Separates comments.
        header_fieldname_order = re.split(delimiter,line)
        self.reorganizeFields(header_fieldname_order)
      if ignore_first_line is True:
        line = f.readline
      for line in f:
        if line.startswith("#"):
          continue
        elif ignore_first_line is True:
          ignore_first_line = False #So that it only ignores 1 line.
          continue
        line = line.rstrip()
        line = line.lstrip()
        comment = None
        pieces = re.split("#", line, 1) #Separates comments.
        line = pieces[0]
        if len(pieces) == 2:
          comment = pieces[1]
        values = re.split(delimiter,line)
        row = MyDataRow(self.fields,values,comment=comment)
        self.rows.append(row)
    if header is True: #Reset the order of the fields back to the original order.
      self.fields = orig_field_order

  def reorganizeFields(self,fieldname_order):
    """
    Reorders the list of fields by the provided fieldname order.
    """
    fieldname_order = list(set(fieldname_order))
    if len(fieldname_order) != len(self.fields):
      raise MyDataException("The provided fieldname order, " + str(fieldname_order) + ", does not have the expected number of elements.")
    new_field_order = []
    for fieldname in fieldname_order:
      for field in self.fields:
        if field.hasName(fieldname):
          new_field_order.append(field)
          continue
    if len(new_field_order) != len(self.fields):
      raise MyDataException("The provided fieldname order, " + str(fieldname_order) + ", does not have the expected number of elements with the correct fieldnames.")
    self.fields = new_field_order

  def getRow(self,row_number):
    """
    Returns the row with the provided row_number.
    """
    return self.rows[row_number]

  def returnTableValue(self,fieldname,row_number=-1):
    """
    Returns the value of the column keyed by fieldname of row row_number.  By default, the
    last value of the table will be returned.
    """
    row = self.getRow(row_number)
    cell = row.getCellWithFieldname(fieldname)
    if cell is None:
      raise MyDataException("No cell for row, " + row_number + ", with the field " + fieldname + ".")
    return cell.getValue()

  def getRowNumbersWithKeyValuePair(self,fieldname,filter_value,comparison = 'eq'):
    """
    Returns a list with the row numbers that have the cells of the given fieldname that has value that matches the comparison.
    """
    row_number_list = []
    for i in range(len(self.rows)):
      value = self.returnTableValue(fieldname,i)
      if comparison == 'eq' and filter_value == value:
        row_number_list.append(i)
      elif comparison == 'gt' and value >  filter_value:
        row_number_list.append(i)
      elif comparison == 'gteq' and value >=  filter_value:
        row_number_list.append(i)
      elif comparison == 'lt' and value <  filter_value:
        row_number_list.append(i)
      elif comparison == 'lteq' and value <=  filter_value:
        row_number_list.append(i)
    return row_number_list

  def returnAllColumnValues(self,fieldname):
    """
    Returns a list of all of the cells values for the given fieldname.
    """
    values= []
    for i in range(len(self.rows)):
      values.append(self.returnTableValue(fieldname,i))
    return values

  def returnMean(self,fieldname):
    """
    Gets the mean of all of the table elements from the provided fieldname
    provided the fieldtype is float or int.
    """
    field = self.getFieldByName(fieldname)
    permitted_types = ["float","int"]
    if field.type not in permitted_types:
      raise MyDataException("Trying to calculate a mean of a field of the incorrect type.  The type of the field is ," + field.type + ", whereas a mean only works on with the following types: " + " ".join(permitted_types) + "")
    return numpy.mean(self.returnAllColumnValues(fieldname))
    
  def returnStandardDeviation(self,fieldname):
    """
    Gets the standard deviation of all of the table elements from the provided key.
    """
    field = self.getFieldByName(fieldname)
    permitted_types = ["float","int"]
    if field.type not in permitted_types:
      raise MyDataException("Trying to calculate a mean of a field of the incorrect type.  The type of the field is ," + field.type + ", whereas a std only works on with the following types: " + " ".join(permitted_types) + "")
    return numpy.std(self.returnAllColumnValues(fieldname))

  def writeToFile(self,filepath,header=False,include_units=False,delimiter=","):
    """
    Writes the table to the given filepath.
    """
    with open(filepath,'w') as f:
      if header is True:
        output = []
        for field in self.fields:
          column = fieldname
          if include_units is True and field.unit is not None:
            column += " (" + field.unit + ")"
          output.append(column)
        f.write(delimiter.join(output))
      for row in self.table:
        output= []
        for field in self.fields:
          cell = row.returnCellWithFieldname(field.name)
          output.append(str(cell))
        f.write(delimiter.join(output))
