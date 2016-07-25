import sys
import os
import math
import argparse
import numpy as np
from scipy import constants

class MyCovarianceMatrix():
  """
  An object to provide how I like to interact with the x,y,z,px,py,pz,E
  covariance matrix.
  """

  def __init__(self,phase_volume):
    self.cov_matrix = self.calcCovarianceMatrix(phase_volume)

  def __len__(self):
    return self.cov_matrix.size
  
  def getCovarianceMatrix(self):
    """
    Wraps the calculation of the covariance matrix.  If it has previously been
    calculated, then skip the calculation.  Returns the 7 x 7 matrix.
    """
    return self.cov_matrix

  def calcCovarianceMatrix(self,phase_volume):
    """
    Uses numpy to calculate the 7 x 7 covariance matrix.
    """
    input_array = []
    for particle in phase_volume:
      row = []
      row.append(particle.x.x)
      row.append(particle.x.y)
      row.append(particle.x.z)
      row.append(particle.p.x)
      row.append(particle.p.y)
      row.append(particle.p.z)
      row.append(particle.getEnergy())
      input_array.append(row)
    input_array = np.array(input_array).T #Transpose
    return np.cov(input_array)

  def getCovarianceElement(self,first_variable,second_variable):
    """
    Returns the appropriate element of the covariance matrix that
    corresponds to Cov(first_variable,second_variable) where both of these
    variables need to be in ["x", "y", "z", "px", "py", "pz", "E"].
    """
    order = ["x", "y", "z", "px", "py", "pz", "E"]
    if first_variable not in order:
      raise CoordinateException("Covariance is only specified between these propetries: " + ", ".join(order))
    if second_variable not in order:
      raise CoordinateException("Covariance is only specified between these propetries: " + ", ".join(order))
    first_index = order.index(first_variable)
    second_index = order.index(second_variable)
    cov_matrix = self.getCovarianceMatrix()
    return cov_matrix[first_index][second_index]

  def printCovarianceMatrix(self,order=["x", "y", "z", "px", "py", "pz", "E"]):
    """
    Returns the covariance matrix with columns/rows in the provided order in upper
    triangular form.
    """
    for o1 in range(len(order)):
      string = ""
      for o2 in range(len(order)):
        if o2 < o1:
          empty_string = "            "
          string += empty_string
        else:
          element1 = order[o1] 
          element2 = order[o2] 
          current_element = self.getCovarianceElement(element1,element2)
          formatted_string = "{:>12.2e}".format(current_element)
          string += formatted_string
      print string
    return

  def printCorrelationMatrix(self,order=["x", "y", "z", "px", "py", "pz", "E"]):
    """
    Returns the correlation matrix with columns/rows in the provided order in upper
    triangular form.
    """
    for o1 in range(len(order)):
      string = ""
      for o2 in range(len(order)):
        if o2 < o1:
          empty_string = "            "
          string += empty_string
        else:
          element1 = order[o1] 
          element2 = order[o2] 
          current_element = self.getCovarianceElement(element1,element2)
          std_1 = np.sqrt(self.getCovarianceElement(element1,element1))
          std_2 = np.sqrt(self.getCovarianceElement(element2,element2))
          formatted_string = "{:>12.3f}".format(current_element/(std_1*std_2))
          string += formatted_string
      print string
    return

  def printMixedMatrix(self,order=["x", "y", "z", "px", "py", "pz", "E"]):
    """
    Returns the correlation matrix with columns/rows in the provided order in upper
    triangular form for the off diagonal values.  The diagonal values are the
    diagonal values of the covariance matrix (i.e. the variances).
    """
    for o1 in range(len(order)):
      string = ""
      for o2 in range(len(order)):
        if o2 < o1:
          empty_string = "                "
          string += empty_string
        else:
          element1 = order[o1] 
          element2 = order[o2] 
          current_element = self.getCovarianceElement(element1,element2)
          if o1 != o2:
            std_1 = np.sqrt(self.getCovarianceElement(element1,element1))
            std_2 = np.sqrt(self.getCovarianceElement(element2,element2))
            formatted_string = "{:>16.3f}".format(current_element/(std_1*std_2))
          elif o1 == o2:
            formatted_string = "{:>16.4e}".format(current_element)
          string += formatted_string
      print string
    return

  def getSubDeterminant(self,subelements=["x","y","z","px","py","pz"]):
    """
    Gets the determinant of the sub matrix defined by the sub elements.
    """
    subkey = ",".join(sorted(subelements))
    if not hasattr(self,"sub_determinant"):
      self.sub_determinant = {}
    if subkey not in self.sub_determinant:
      order = ["x","y","z","px","py","pz","E"]
      subindices = []
      for element in subelements:
        subindices.append(order.index(element))
      rows = np.array(subindices, dtype=np.intp)
      rows = rows[:, np.newaxis]
      columns = np.array(subindices, dtype=np.intp)
      (sign, logdet) = np.linalg.slogdet(self.cov_matrix[rows,columns])
      self.sub_determinant[subkey] = sign * np.exp(logdet)
    return self.sub_determinant[subkey]

  def writeCovarianceMatrix(self,filepath,sep=" "):
    """
    Writes the covariance matrix to a text file.
    """
    cov_matrix = self.getCovarianceMatrix()
    with open(filepath,'w') as f:
      for i in range(7):
        output = []
        for j in range(7):
          if j<i:
            continue
          output.append(cov_matrix[i][j])
        f.write(sep.join([str(o) for o in output])+"\n")
    
class MyEditableCovarianceMatrix(MyCovarianceMatrix):
  """
  An object to provide how I like to interact with the x,y,z,px,py,pz,E
  covariance matrix initially setting all elements to 0 and allowing writing
  to individual elements 
  """

  def __init__(self):
    self.cov_matrix = np.zeros(shape=(7,7))

  def setCovarianceElement(self,first_variable,second_variable,value):
    """
    Sets the appropriate element of the covariance matrix that
    corresponds to Cov(first_variable,second_variable) where both of these
    variables need to be in ["x", "y", "z", "px", "py", "pz", "E"].
    """
    order = ["x", "y", "z", "px", "py", "pz", "E"]
    if first_variable not in order:
      raise CoordinateException("Covariance is only specified between these propetries: " + ", ".join(order))
    if second_variable not in order:
      raise CoordinateException("Covariance is only specified between these propetries: " + ", ".join(order))
    first_index = order.index(first_variable)
    second_index = order.index(second_variable)
    self.cov_matrix[first_index,second_index] = value

  def readCovarianceMatrix(self,filepath,sep=" "):
    """
    Read covariance matrix from a text file into object.
    """
    order = ["x", "y", "z", "px", "py", "pz", "E"]
    cov_matrix = self.getCovarianceMatrix()
    with open(filepath,'r') as f:
      i = 0
      for line in f:
        line.rstrip()
        elements = line.split(sep)
        for j in range(i,i+len(elements)):
          self.setCovarianceElement(order[i],order[j],elements[j-i])
          self.setCovarianceElement(order[j],order[i],elements[j-i])
        i += 1
